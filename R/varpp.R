#'varpp: VARiant Prioritisation by Phenotype
#'
#' @param HPO_genes HPO term associated list of genes
#' @param type the prediction data; either hcl (single cell) or gtex (tissue specific)
#' @param ntree is the number of trees that should be built for ranger. It defaults to 1000
#' @param max.depth is the maximum tree depth for the ranger trees. IT defaults to 3.
#' @param cores number of cores for parallel, defaults to 4
#' @import progress
#' @import doMC
#' @import tidyverse
#' @import ranger
#' @import magrittr
#' @importFrom magrittr "%>%"
#' @import foreach
#' @importFrom iterators icount
#' @export
varpp <- function(HPO_genes,
                  type="gtex",
                  ntree=500,
                  max.depth=NULL,
                  cores=4){

  #=====================================================================================
  # Prepare the data
  #=====================================================================================

  if(type %in% "gtex"){

    patho  <- patho_gtex
    benign <- benign_gtex
  } else {
    patho  = patho_hcl
    benign = benign_hcl
  }


  message(paste0("VARPP initiated with ", ntree," trees, maximum depth of ", max.depth," and ",type," data."))
  # Get the gene names
  hpo_gene_names <- HPO_genes

  # Filter the genes that we got from phenolyzer
  patho %>%
    filter(Gene %in% hpo_gene_names) %>%
    select(-c(CADD_raw_rankscore)) %>%
    rename(CADD_raw_rankscore = CADD_PHRED_SCORE) -> varpp_patho

  varpp_patho <- varpp_patho[,c(1,2,4,3,7:length(names(varpp_patho)))]


  # Filter out the benign genes that are in the pathogenic gene list
  benign %>%
    filter(!Gene %in% intersect(benign$Gene, patho$Gene)) %>%
    select(-c(CADD_raw_rankscore)) %>%
    rename(CADD_raw_rankscore = CADD_PHRED_SCORE) -> varpp_benign

  varpp_benign <- varpp_benign[,c(1,2,4,3,7:length(names(varpp_benign)))]

  # Create the input data for varppRuleFit
  dat <- list(dat=data.frame(rbind(varpp_patho, varpp_benign)), disease_variants=data.frame(varpp_patho), benign_variants=data.frame(varpp_benign))




  # Spedify the benign and pathogenic gene names
  cls_pathogenic_genes = unique(dat$dat$Gene[dat$dat$Pathogenic %in% 1])
  cls_benign_genes     = unique(dat$dat$Gene[dat$dat$Pathogenic %in% 0])

  #=====================================================================================
  # SET UP RESULTS TABLE
  #=====================================================================================

  # These will be progressively filled in loop for random forest votes
  rf_trees <- list()

  rf_trees$rf_results <- data.frame(
    bind_rows(dat$disease_variants[ , c("GeneVariant","Pathogenic")], dat$benign_variants[ , c("GeneVariant","Pathogenic")]),
    predictNum=vector(length=nrow(dat$disease_variants) + nrow(dat$benign_variants), mode="numeric"),
    predictDenom=vector(length=nrow(dat$disease_variants) + nrow(dat$benign_variants), mode="numeric"))

  rf_trees$rf_results_varimp <- data.frame(
    Variable=c("CADD_raw_rankscore", colnames(dat$dat)[!colnames(dat$dat) %in% c("Gene", "GeneVariant", "Pathogenic", "CADD_raw_rankscore")]),
    sumVarimp=vector(length=length(colnames(dat$dat)[!colnames(dat$dat) %in% c("Gene", "GeneVariant", "Pathogenic", "CADD_raw_rankscore")]) + 1, mode="numeric"),
    ntree=vector(length=length(colnames(dat$dat)[!colnames(dat$dat) %in%  c("Gene", "GeneVariant", "Pathogenic", "CADD_raw_rankscore")]) + 1, mode="numeric"),
    stringsAsFactors=FALSE)

  rf_trees$rules <- list()


  #==========================================================================================
  # THE SAMPLING LOOP
  #==========================================================================================

  # Create an empty list to be filled in the loop
  inbag = NULL

  # Sampling
  message("Starting varpp...")

  cl <- makeCluster(cores)
  doParallel::registerDoParallel(cl)

  inbag <- foreach::foreach (j = 1:ntree) %dopar% {


    # Sampling the pathogenic variants
    cls_patho <- sample(cls_pathogenic_genes, replace=TRUE)

    sub_patho <- lapply(cls_patho, FUN=function(x) dat$disease_variants[which(dat$disease_variants$Gene == x), "GeneVariant"])
    sub_patho <- lapply(sub_patho, FUN=function(x) sample(x, size=1))
    sub_patho <- dplyr::combine(sub_patho)
    sub_patho <- sub_patho[match(gsub("_.*$", "", sub_patho), gsub("_.*$", "", sub_patho))]

    # Sampling the benign variants
    cls_benign  <- sample(cls_benign_genes, replace=TRUE)

    benign_variants_sub <- dat$benign_variants[dat$benign_variants$Gene %in% unique(cls_benign), c("Gene","GeneVariant","Pathogenic","CADD_raw_rankscore")]

    # The new sampling step based on the function in utilities
    sub_benign <- .sample_benign_variants(benign_variants_sub, cls_benign)

    # This is still necessary to make sure that only one variant is ever selected per gene
    sub_benign <- sub_benign[match(gsub("_.*$", "", sub_benign), gsub("_.*$", "", sub_benign))]

    WEIGHTS        = rbind(data.frame(plyr::count(sub_benign)),data.frame(plyr::count(sub_patho)))

    # The merge function does not preserve the order of the original 'dat' data.frame.
    # That means, the sampling inbag order I use to subset later is not in the right order
    # and the wrong Gene Variants are selected. Solution: create index id and sort the merged data after

    dat$dat$id     <- seq_len(nrow(dat$dat))
    weighted_data  <- merge(dat$dat, WEIGHTS, by.x="GeneVariant", by.y="x", all=T)
    dat$dat$id     <- NULL
    weighted_data  <- weighted_data[order(weighted_data$id), ] %>% select(-id)
    inbag <-
      as.numeric(unlist(weighted_data %>%
                          mutate(freq = ifelse(freq %in% NA, 0, freq)) %>%
                          rename(weight = freq) %>%
                          select(weight)))

    # Remove redundancies
    rm(sub_patho, sub_benign, benign_variants_sub, WEIGHTS)


    # This step is crucial: the index so far only collects the position in the original data
    # set, but some Variants are sampled multiple times, hence we need to duplciate the selection of
    # those variants n times, as extracted by the count in hte previous loop

    indices   <- rep(which(inbag >= 1 ), inbag[which(inbag >= 1 )])

    # Specify the in bag samples: based on the indices
    dat_in        <- dat$dat[indices,]
    dat_in$weight <- 1
    dat_in$Gene   <- NULL

    # Remove those variants that are in the same genes as the ones selected in this round of the loop
    dat_out        <- dat$dat[!(gsub("_.*$", "",dat$dat$GeneVariant) %in% gsub("_.*$", "",dat_in$GeneVariant)), ]
    dat_out$weight <- 0
    dat_out$Gene <- NULL

    dat_boot <- list(training=data.frame(bind_rows(dat_in, dat_out), check.names=FALSE))

    rm(dat_in, dat_out)

    dat_boot$training$Pathogenic <- factor(dat_boot$training$Pathogenic)

    # Random forest does not handle missing data
    dat_boot$training_CADD <- na.omit(dat_boot$training[ , colnames(dat_boot$training) != "MetaSVM_rankscore"])
    dat_boot$training_CADD <- dat_boot$training_CADD %>% arrange(desc = GeneVariant)
    dat_boot$training <- NULL


    rf_results <-
      ranger::ranger(
        data = dat_boot$training_CADD[ , !colnames(dat_boot$training_CADD) %in% c("GeneVariant", "weight")],
        dependent.variable.name = "Pathogenic",
        num.trees = 1,
        replace = FALSE,
        sample.fraction = 0.9999,
        importance = "permutation",
        case.weights = dat_boot$training_CADD$weight,
        scale.permutation.importance = FALSE,
        holdout = TRUE,
        max.depth = max.depth
      )


    rf_trees$rf_results[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$rf_results$GeneVariant), "predictNum"] <- rf_trees$rf_results[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$rf_results$GeneVariant), "predictNum"] + as.numeric(as.character(rf_results$predictions[dat_boot$training_CADD$weight == 0]))
    rf_trees$rf_results[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$rf_results$GeneVariant), "predictDenom"] <- rf_trees$rf_results[match(dat_boot$training_CADD$GeneVariant[dat_boot$training_CADD$weight == 0], rf_trees$rf_results$GeneVariant), "predictDenom"] + 1
    rf_trees$rf_results_varimp[match(names(rf_results$variable.importance), rf_trees$rf_results_varimp$Variable), "sumVarimp"] <- rf_trees$rf_results_varimp$sumVarimp + rf_results$variable.importance
    rf_trees$rf_results_varimp[match(names(rf_results$variable.importance), rf_trees$rf_results_varimp$Variable), "ntree"] <- rf_trees$rf_results_varimp$ntree + ifelse(rf_results$variable.importance == 0, 0, 1)
    rf_trees$rules <- .extract_ranger_rules(rf_results)

    # Remove redundancies from memory
    rm(dat_boot)

    rf_trees

  }
  stopCluster(cl)

  # Collect the results for ranger at every step and put it together

  CADD   = list()
  varIMP = list()
  RULES  = list()

  for (i in 1:ntree){

    CADD[[i]]   <- inbag[[i]]$rf_results[,c(3:4)]
    varIMP[[i]] <- inbag[[i]]$rf_results_varimp[,c(2:3)]
    RULES[[i]]  <- inbag[[i]]$rules
  }

  CADD   <- Reduce('+', CADD)
  varIMP <- Reduce('+', varIMP)

  rf_trees                   <- list()
  rf_trees$rf_results        <- cbind(inbag[[1]]$rf_results[,c(1,2)], CADD)
  rf_trees$rf_results_varimp <- cbind(Variable= inbag[[1]]$rf_results_varimp[,1], varIMP)

  rm(CADD,varIMP)

  # Drop variants not seen by the classifier
  rf_trees$rf_results <- rf_trees$rf_results[rf_trees$rf_results$predictDenom != 0, ]


  accuracy <- merge(data.frame(GeneVariant=rf_trees$rf_results$GeneVariant,
                               Pathogenic=rf_trees$rf_results$Pathogenic,
                               rf_results=rf_trees$rf_results$predictNum/rf_trees$rf_results$predictDenom),
                    bind_rows(dat$disease_variants[ , c("GeneVariant", "CADD_raw_rankscore")],
                                        dat$benign_variants[ , c("GeneVariant", "CADD_raw_rankscore")]),
                    by.x="GeneVariant", by.y="GeneVariant", all.x=TRUE, all.y=FALSE)

  varimp <- data.frame(Variable=rf_trees$rf_results_varimp$Variable,
                       rf_results=rf_trees$rf_results_varimp$sumVarimp/rf_trees$rf_results_varimp$ntree)#,

  # Combine the rules
  RULES        <- do.call(c, RULES) #Reduce(c,inbag)
  names(RULES) <- paste0("rule_",1:length(RULES))


  results <- list(accuracy=accuracy, varimp=varimp, rules=RULES)
  class(results)<- "varpp"
  message("VARPP finished!")
  return(results)

}
