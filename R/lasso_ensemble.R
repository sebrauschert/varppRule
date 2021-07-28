#' LASSO cross validation of rules
#'
#' This function performs nested cross validation on the generated rule data set. It is the final step in the rulefit algorithm and returns the predictions
#'
#' @param data is a list of data with the rules added. The benign and the pathogenic variants files are necessary for the sampling.
#' @param rules is the list of rules that were generated in the varpp function. This is necessary for the annotation of the final results.
#' @param bootstrap.rounds number of bootstrap rounds for the outer loop of the LASSO cross-validation, defaults to 100.
#' @param cores number of cores for parallel, defaults to 4
#' @return A list of predictions for the CADD raw rankscore and the tissue/cell specific expression added. Further, a variable importance list for all rules and variables tested.
#'
#' @import caret
#' @import glmnet
#' @import tidyverse
#' @import magrittr
#' @importFrom caret trainControl train nearZeroVar
#' @import doMC
#' @import doParallel
#' @import parallel
#' @import grid
#' @export
lasso_ensemble <- function(data,
                             rules,
                             bootstrap.rounds,
                             cores){


  # Specify benign and pathogenic gene names for sampling
  cls_pathogenic_genes <- unique(data$dat$Gene[data$dat$Pathogenic %in% 1])
  cls_benign_genes     <- unique(data$dat$Gene[data$dat$Pathogenic %in% 0])


  #=====================================================================================
  # SET UP RESULTS TABLE
  #=====================================================================================

  # These will be progressively filled in loop for random forest votes
  rf_trees <- list()

  rf_trees$CADD_expression <- data.frame(
    bind_rows(data$disease_variants[ , c("GeneVariant","Pathogenic")], data$benign_variants[ , c("GeneVariant","Pathogenic")]),
    predictNum=vector(length=nrow(data$disease_variants) + nrow(data$benign_variants), mode="numeric"),
    predictDenom=vector(length=nrow(data$disease_variants) + nrow(data$benign_variants), mode="numeric"))

  rf_trees$CADD_expression_varimp <- data.frame(
    Variable=c("CADD_raw_rankscore", colnames(data$dat)[!colnames(data$dat) %in% "Gene"]),
    sumVarimp=vector(length=length(colnames(data$dat)[!colnames(data$dat) %in%  "Gene"]) + 1, mode="numeric"),
    ntree=vector(length=length(colnames(data$dat)[!colnames(data$dat) %in% "Gene"]) + 1, mode="numeric"),
    stringsAsFactors=FALSE)



  #==========================================================================================
  # THE SAMPLING LOOP
  #==========================================================================================

  # Create an empty list to be filled in the loop
  inbag = NULL

  # Sampling
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  message(paste0("Cross validating LASSO with k: ",bootstrap.rounds))

  inbag <- foreach (j = 1:bootstrap.rounds) %dopar% {

    cls_patho <- sample(cls_pathogenic_genes, replace=TRUE)

    sub_patho <- lapply(cls_patho, FUN=function(x) data$disease_variants[which(data$disease_variants$Gene == x), "GeneVariant"])
    sub_patho <- lapply(sub_patho, FUN=function(x) sample(x, size=1))
    sub_patho <- dplyr::combine(sub_patho)
    sub_patho <- sub_patho[match(gsub("_.*$", "", sub_patho), gsub("_.*$", "", sub_patho))]

    # Sampling the benign variants
    cls_benign  <- sample(cls_benign_genes, replace=TRUE)

    benign_variants_sub <- data$benign_variants[data$benign_variants$Gene %in% unique(cls_benign),
                                                c("Gene","GeneVariant","Pathogenic","CADD_raw_rankscore")]

    # The new sampling step based on the function in utilities
    sub_benign <- .sample_benign_variants(benign_variants_sub, cls_benign)

    # This is still necessary to make sure that only one variant is ever selected per gene
    sub_benign <- sub_benign[match(gsub("_.*$", "", sub_benign), gsub("_.*$", "", sub_benign))]
    WEIGHTS    <- rbind(data.frame(plyr::count(sub_benign)),data.frame(plyr::count(sub_patho)))

    # The merge function does not preserve the order of the original 'dat' data.frame.
    # That means, the sampling inbag order I use to subset later is not in the right order
    # and the wrong Gene Variants are selected. Solution: create index id and sort the merged data after

    data$dat$id   <- seq_len(nrow(data$dat))
    weighted_data <- merge(data$dat, WEIGHTS, by.x="GeneVariant", by.y="x", all=T)
    data$dat$id   <- NULL
    weighted_data <- weighted_data[order(weighted_data$id), ] %>% select(-id)
    inbag         <-
      as.numeric(unlist(weighted_data %>%
                          #select(-x) %>%
                          mutate(freq = ifelse(freq %in% NA, 0, freq)) %>%
                          rename(weight = freq) %>%
                          select(weight)))


    # Remove redundancies
    rm(sub_patho, sub_benign, benign_variants_sub, WEIGHTS)
    gc()

    # This step is crucial: the index so far only collects the position in the original data
    # set, but some Variants are sampled multiple times, hence we need to duplciate the selection of
    # those variants n times, as extracted by the count in hte previous loop

    indices   <- rep(which(inbag >= 1 ), inbag[which(inbag >= 1 )])

    # Specify the in bag samples: based on the indices
    dat_in        <- data$dat[indices,]
    dat_in$weight <- 1
    dat_in$Gene   <- NULL

    # Remove those variants that are in the same genes as the ones selected in this round of the loop
    dat_out        <- data$dat[!(gsub("_.*$", "",data$dat$GeneVariant) %in% gsub("_.*$", "",dat_in$GeneVariant)), ]
    dat_out$weight <- 0
    dat_out$Gene   <- NULL


    dat_boot <- list(training=data.frame(bind_rows(dat_in, dat_out), check.names=FALSE))

    rm(dat_in, dat_out)
    gc()

    dat_boot$training$Pathogenic <- factor(dat_boot$training$Pathogenic)

    # Random forest does not handle missing data
    dat_boot$training_CADD <- na.omit(dat_boot$training[ , colnames(dat_boot$training) != "MetaSVM_rankscore"])
    dat_boot$training_CADD <- dat_boot$training_CADD %>% arrange(desc = GeneVariant)
    dat_boot$training      <- NULL

    # Split into features and class: training set
    # The long expression that gets assigned to X_train is the training data set:
    # training = dat_boot$training_CADD %>% filter(weight %in% 1) %>% select(-c(GeneVariant, weight)).
    # It is just more memory efficient than creating a "training" object first.

    X_train  <- (dat_boot$training_CADD %>% filter(weight %in% 1) %>% select(-c(GeneVariant, weight)))[,c(2:length(names(dat_boot$training_CADD %>% filter(weight %in% 1) %>% select(-c(GeneVariant, weight)))))]
    Y_train  <- as.factor((dat_boot$training_CADD %>% filter(weight %in% 1))[,"Pathogenic"])
    X_train  <- Matrix::as.matrix(X_train)


    modelling_start_time <- Sys.time()

    model_glmnet <- cv.glmnet(X_train,
                              Y_train,
                              alpha = 1,
                              family = "binomial",
                              parallel = FALSE,
                              standardize= TRUE)



    modelling_end_time   = Sys.time()
    modelling_time_taken = modelling_end_time - modelling_start_time
    message(paste0("Modelling finished after ", paste(round(modelling_time_taken,2), units(modelling_time_taken), sep=" ")))


    # Predict on the data
    predictions <- as.factor(predict(model_glmnet,
                                     Matrix::as.matrix(dat_boot$training_CADD[,4:(length(names(dat_boot$training_CADD)))-1]),
                                     s="lambda.min",
                                     type="class"))


    Coefs <- coef(model_glmnet, s="lambda.min")

    Results <- data.frame(
      features = Coefs@Dimnames[[1]][ which(Coefs[,1] != 0 ) ], #intercept included
      coefs    = Coefs              [ which(Coefs[,1] != 0 ) ]  #intercept included
    )

    Results$rules <- rules[as.character(Results$features)]
    Results %>%
      arrange(desc(abs(coefs)))  %>% filter(!features %in% "(Intercept)") -> varimp

    rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictNum"] <- rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictNum"] + as.numeric(as.character(predictions[dat_boot$training_CADD$weight == 0]))

    # Here we keep track of how and if a GeneVariant was "predicted" by teh classifier or not. Later, if this is 0, the Variant will be excluded.
    rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictDenom"] <- rf_trees$CADD_expression[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$CADD_expression$GeneVariant), "predictDenom"] + 1
    rf_trees$CADD_expression_varimp[match(as.character(varimp$features), rf_trees$CADD_expression_varimp$Variable), "sumVarimp"] <- rf_trees$CADD_expression_varimp[match(as.character(varimp$features), rf_trees$CADD_expression_varimp$Variable), "sumVarimp"] + as.numeric((varimp$coefs))
    rf_trees$CADD_expression_varimp[match((varimp$features), rf_trees$CADD_expression_varimp$Variable), "ntree"] <- rf_trees$CADD_expression_varimp[match((varimp$features), rf_trees$CADD_expression_varimp$Variable), "ntree"] + ifelse(varimp$coefs == 0, 0, 1)

     # Remove redundancies from memory
    rm(dat_boot, Coefs, Results, varimp)
    gc()


    rf_trees
  }

  stopCluster(cl)

  CADD   = list()
  varIMP = list()

  for (i in 1:bootstrap.rounds){

    CADD[[i]]   <- inbag[[i]]$CADD_expression[,c(3:4)]
    varIMP[[i]] <- inbag[[i]]$CADD_expression_varimp[,c(2:3)]

  }

  CADD             <- Reduce('+', CADD)
  varIMP           <- Reduce('+', varIMP)
  CADD_expression  <- cbind(inbag[[i]]$CADD_expression[,c(1,2)], CADD)


  rf_trees                        <- list()
  rf_trees$CADD_expression        <- CADD_expression
  rf_trees$CADD_expression_varimp <- cbind(Variable= inbag[[i]]$CADD_expression_varimp[,1], varIMP)

  # Drop variants not seen by the classifier
  rf_trees$CADD_expression <- rf_trees$CADD_expression[rf_trees$CADD_expression$predictDenom != 0, ]

  accuracy_CADD_expression <- data.frame(GeneVariant=rf_trees$CADD_expression$GeneVariant,
                                         Pathogenic=rf_trees$CADD_expression$Pathogenic,
                                         CADD_expression=rf_trees$CADD_expression$predictNum/rf_trees$CADD_expression$predictDenom)


  # Add CADD and MetaSVM scores for comparison to VARPP
  accuracy <- merge(accuracy_CADD_expression, bind_rows(data$disease_variants[ , c("GeneVariant", "CADD_raw_rankscore")],
                                                        data$benign_variants[ , c("GeneVariant", "CADD_raw_rankscore")]),
                    by.x="GeneVariant", by.y="GeneVariant", all.x=TRUE, all.y=FALSE)

  varimp <- data.frame(Variable = rf_trees$CADD_expression_varimp[,1],
                       CADD_expression = rf_trees$CADD_expression_varimp[,2]/rf_trees$CADD_expression_varimp$ntree)



  # Filter out the variables that did not get picked for variable importance, add the rules and reorder the columns
  varimp       <- varimp %>% filter(!CADD_expression %in% 0)
  varimp$Rules <- rules[as.character(varimp$Variable)]
  varimp       <- varimp[,c("Variable", "Rules", "CADD_expression")]
  varimp       <- varimp %>% arrange(desc(abs(as.numeric(as.character(CADD_expression)))))

  list(accuracy = accuracy,
       varimp = varimp)
}




#===========#
#           #
#  /(_M_)\  #
# |       | #
#  \/~V~\/  #
#           #
#===========#
