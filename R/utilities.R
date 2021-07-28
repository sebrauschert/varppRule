#' Class specific functions
#'
#' @param x, an object of class varppRuleFit
#' @export
print.varppRule <- function(x){

  if(!class(x) =="varppRule"){
    stop("Argument 'x' should be of class 'varppRuleFit'")
  }

  cat("=========================================\n")
  cat("\n")

  cat(c("", "  _____       _      ______ _ _", " |  __ \\     | |    |  ____(_) |",
    " | |__) |   _| | ___| |__   _| |_", " |  _  / | | | |/ _ \\  __| | | __|",
    " | | \\ \\ |_| | |  __/ |    | | |_", " |_|  \\_\\__,_|_|\\___|_|    |_|\\__|",
    "", ""), sep="\n")

  cat("\n")
  cat("========== Final RuleFit model ==========\n")


  cat("\n  Cohen's Kappa: ", as.numeric(metrics(x)$RuleFit$overall["Kappa"]),
      "\n  F1-score:", as.numeric(metrics(x)$RuleFit$byClass["F1"]),
      "\n  Number of rules :", length(grep("rule",x$RuleFit$varimp$Variable)),
      "\n  Top Rule :", as.character(x$RuleFit$varimp$Variable[grep("rule",x$RuleFit$varimp$Variable)][1]),
      "\n  Prediction finished after :", x$duration[[1]], units(x$duration), "\n", sep = " ")

  cat("\n")

  cat("=========================================\n")

  cat("\n")
  cat("Top predictors with performance:\n")

  cat("\n")
  print(head(x$RuleFit$varimp[grep("rule",x$RuleFit$varimp$Variable),], n=10), print.gap = 2, quote = FALSE, row.names = FALSE)
  invisible(x$RuleFit$varimp)

  cat("\n")

  cat("=========================================\n")

  cat("\n")
  cat("Model call: \n")
  print(x$rulefit_call)
  cat("\n")
}


#' Class specific functions
#'
#' @param x, an object of class varpp
#' @export
print.varpp <- function(x){

  if(!class(x) =="varpp"){
    stop("Argument 'x' should be of class 'varpp'")
  }

  cat("=========================================\n")

  cat(c("", " __   ___   ___ ___ ___ ", " \\ \\ / /_\\ | _ \\ _ \\ _ \\",
        "  \\ V / _ \\|   /  _/  _/", "   \\_/_/ \\_\\_|_\\_| |_|  ",
        "                        "), sep="\n")

  cat("\n")

  cat("Top 10 predicted Variants:\n")

  cat("\n")

  print(x$accuracy %>% arrange(desc(rf_results)) %>% head(n=10) %>% rename(VARPP_Score = rf_results), print.gap = 2, quote = FALSE, row.names = FALSE)
  cat("\n")

  cat("========== Final VARPP model ==========\n")


  cat("\n  Cohen's Kappa: ", as.numeric(.class_by_threshold(x$accuracy[,2], x$accuracy[,3])$overall["Kappa"]),
      "\n  F1-score:", as.numeric(.class_by_threshold(x$accuracy[,2], x$accuracy[,3])$byClass["F1"]))

  cat("\n")

  cat("\n")

  cat("======== Model performance ============\n")

  cat("\n")

  print(performance_varpp(x), print.gap = 2, quote = FALSE, row.names = FALSE)

  cat("\n")

  cat("=========================================\n")
  cat("\n")
  cat("Top predictors with importance score:\n")

  cat("\n")
  print(head(x$varimp, n=10)%>%arrange(desc(rf_results))%>%rename(Importance=rf_results), print.gap = 2, quote = FALSE, row.names = FALSE)
  invisible(x$varimp)

  cat("\n")

  cat("=========================================\n")
}



#' Predict function for varppRule
#'
#' @param patient_data based on a patient .vcf file,
#'                     a preprocessed inpout file that is annotated with GTEx and CADD scores
#' @param hpo_term patient hpo terms
#' @export
predict.varppRule <- function(patient_data, model_results, predict=c("probability", "class")){

  results = model_results
  attr(results, "class") <- "varppRule" # This is legacy, as all results have been created with an old attribute name

  # Prepare coefficients
  rules <- results$RuleFit$varimp[grep("rule",results$RuleFit$varimp$Variable),] %>%
    select(Description) %>%
    unlist() %>%
    as.character()

  coefficients <- results$RuleFit$varimp[grep("rule",results$RuleFit$varimp$Variable),] %>%
    select(Coefficient) %>%
    unlist() %>%
    as.numeric()



  # Create Rule variables based on the patient data
  for( i in 1:length(rules)){

    if(grepl('CADD_raw_rankscore >', rule)) {
      patient_data <-  cbind(patient_data, patient_data %>%
                               mutate(!!rules[i] := ifelse(eval(parse(text=rules[i])), 1, 0 )) %>%
                               select(!!rules[i]) )

    }else{

      patient_data <- cbind(patient_data,patient_data %>%
                              mutate(!!rules[i] := ifelse(eval(parse(text=rules[i])), 0, 1 )) %>%
                              select(!!rules[i]))
    }

  }

  # Prediction step
  features <- as.character(results$RuleFit$varimp$Variable)

  # Make sure the rule variables are called the same in the feature set as they were called above in the matrix
  features[grep("rule", results$RuleFit$varimp$Variable)] <-
    as.character(results$RuleFit$varimp$Description[grep("rule", results$RuleFit$varimp$Variable)])

  feature_vector <- results$RuleFit$varimp$Coefficient
  names(feature_vector) <- features

  # Here we return the probabilities for every variant to be pathogenic
  patient_data$Prediction <- logit2prob(rowSums(as.matrix(patient_data[,features]) %*%diag(feature_vector)))

  if(predict=="probability"){

    patient_data[,c("Gene", "CADD_raw_rankscore", "Prediction")]

  }else{

    patient_data %>%
      mutate(Prediction = ifelse(Prediction > 0.5, "Pathogenic", "Benign")) %>%
      select(Gene, CADD_raw_rankscore, Prediction)
  }


}


#====================================================================================================================


#' Sampling and sub-setting for the benign variant data
#'
#' @param benign_data the subset of benign variants
#' @param sampled_genes the sampled genes (with replacement) from the sampling step
#'
#' @return a list of randomly sampled gene variants, without replacement
#' @importFrom plyr count
#' @importFrom  magrittr "%>%"
#'
#' @export
.sample_benign_variants <- function(benign_data, sampled_genes){

  # Create a table that has the gene name and the number of variants as a second column
  # we want this as this will allow the sampling of a random number between 1 and the number
  # of gene variants.
  count_table   <- plyr::count(benign_data$Gene)

  # Now we merge the count data with the genes that were sampled in the step before this one
  benign_subset <- merge(count_table, data.frame(Gene = sampled_genes), by.x="x", by.y="Gene")
  gene_variants <- benign_subset %>%
    mutate(variant = paste0(x, "_", unlist(lapply(benign_subset$freq, function(x) sample(1:x, 1))))) %>%
    select(variant) %>%
    unlist() %>%
    as.character()

  gene_variants
}

#====================================================================================================================

#' Area under the precision recall curve plot for Class \code{varppRuleFit}
#' @param x an object of class \code{varppRuleFit}
#' @import ggplot2
#' @import precrec
#' @export
auPRC <- function(x, model = c("rf", "rulefit")){

  if(model == "rf"){
    randomforest <- evalmod(scores = x$RandomForest$accuracy[,3], labels = x$RandomForest$accuracy[,2])
    autoplot(randomforest, "PRC")
  }else{
    rulefit <- evalmod(scores = x$RuleFit$accuracy[,3], labels = x$RuleFit$accuracy[,2])
    autoplot(rulefit,"PRC")

  }

}

#' Area under the receiver operator curve plot for Class \code{varppRuleFit}
#' @param x an object of class \code{varppRuleFit}
#' @import ggplot2
#' @import precrec
#' @export
auROC <- function(x, model = c("rf", "rulefit")){

  if(model == "rf"){
    randomforest <- evalmod(scores = x$RandomForest$accuracy[,3], labels = x$RandomForest$accuracy[,2])
    autoplot(randomforest, "ROC")
  }else{
    rulefit <- evalmod(scores = x$RuleFit$accuracy[,3], labels = x$RuleFit$accuracy[,2])
    autoplot(rulefit, "ROC")

  }

}

#====================================================================================================================

#' Function to return all rules ranked by variable importance
#' @param x an object of class varppRuleFit
#' @export
ruleVarImp <- function(x){

  x$RuleFit$varimp

}

#====================================================================================================================

#' Return model metrics for both Random Forest and RuleFit
#' @param actual values and predicted values of the model
#' @import pROC
#' @export
.threshold <- function(actual, predicted){
  suppressMessages(suppressWarnings(
    roc_obj <- pROC::roc(actual, predicted)))

  #get the "best" "threshold"
  # there are lots of other options for other metrics as well
  suppressMessages(suppressWarnings(
    cut_off <- pROC::coords(roc_obj, "best", "threshold")$threshold))
  cut_off

}

#====================================================================================================================

#' Return model metrics for both Random Forest and RuleFit
#' @param actual values and predicted values of the model
#' @import pROC
#' @export
.class_by_threshold <- function(actual, predicted) {

  suppressMessages(suppressWarnings(
    roc_obj <- pROC::roc(actual, predicted)))

  #get the "best" "threshold"
  # there are lots of other options for other metrics as well
  suppressMessages(suppressWarnings(
    cut_off <- pROC::coords(roc_obj, "best", "threshold")$threshold))

  tab <- table(predicted = ifelse(predicted > cut_off , 1, 0), original = actual)
  if(nrow(tab)!=ncol(tab)){

    missings      <- setdiff(colnames(tab),rownames(tab))
    missing_mat   <- mat.or.vec(nr = length(missings), nc = ncol(tab))
    tab           <- as.table(rbind(as.matrix(tab), missing_mat))
    rownames(tab) <- colnames(tab)
  }

  caret::confusionMatrix(tab)
}

#====================================================================================================================

#' Return model metrics for both Random Forest and RuleFit
#' @param x an object of class \code{varppRuleFit}
#' @import pROC
#' @importFrom caret confusionMatrix
#' @export
metrics <- function(x) {
  #========================================================================================
  # Model quality measures
  # 1. Random Forest predictions
  #========================================================================================

  if(!is.null(x$RandomForest$accuracy[,3])){

    RandomForestPrediction <- .class_by_threshold(x$RandomForest$accuracy[,2], x$RandomForest$accuracy[,3])

    #========================================================================================
    # 2. LASSO
    #========================================================================================

    RuleFitPrediction <-  .class_by_threshold(x$RuleFit$accuracy[,2], x$RuleFit$accuracy[,3])

    list(RandomForest = RandomForestPrediction, RuleFit = RuleFitPrediction)
  }else{

    RuleFitPrediction <-  .class_by_threshold(x$RuleFit$accuracy[,2], x$RuleFit$accuracy[,3])

    list(RuleFit=RuleFitPrediction)


  }
}


#====================================================================================================================

#' Return a table with model names, auPRC, PP100 and ntree of the model: only for two level bootstrap model
#' @param x an object of class \code{varppRuleFit}
#' @import precrec
#' @export
performance <- function(x, ntree=x$ntree){
  results_varpp <-
    evalmod(
      mmdata(scores=join_scores(x$RuleFit$accuracy$CADD_expression[!is.na(x$RuleFit$accuracy$CADD_expression)],
                                x$RuleFit$accuracy$CADD_raw_rankscore[!is.na(x$RuleFit$accuracy$CADD_raw_rankscore)],
                                chklen=FALSE),

             join_labels(x$RuleFit$accuracy$Pathogenic[!is.na(x$RuleFit$accuracy$CADD_expression)],
                         x$RuleFit$accuracy$Pathogenic[!is.na(x$RuleFit$accuracy$CADD_raw_rankscore)],
                         chklen=FALSE),
             modnames=c("RuleFit","CADD_raw_rankscore"),
             dsids=1:2)
    )


  model_results <-     auc(results_varpp) %>%
    filter(curvetypes %in% "PRC") %>%
    select(-c(dsids,curvetypes)) %>%
    mutate(PP100 = c(sum(x$RuleFit$accuracy[order(-(x$RuleFit$accuracy$CADD_expression)), ][1:100, "Pathogenic"])/100,
                     sum(x$RuleFit$accuracy[order(-(x$RuleFit$accuracy$CADD_raw_rankscore)), ][1:100, "Pathogenic"])/100)) %>%
    mutate(ntree = ntree) %>%
    rename(auPRC = aucs)

  model_results
}

#====================================================================================================================

#' Return a table with model names, auPRC, PP100 and ntree of the model: only for two level bootstrap model
#' @param x an object of class \code{varppRuleFit}
#' @import precrec
#' @export
performance_varpp <- function(x){
  results_varpp <-
    evalmod(
      mmdata(scores=join_scores(x$accuracy$rf_results[!is.na(x$accuracy$rf_results)],
                                x$accuracy$CADD_raw_rankscore[!is.na(x$accuracy$CADD_raw_rankscore)],
                                chklen=FALSE),

             join_labels(x$accuracy$Pathogenic[!is.na(x$accuracy$rf_results)],
                         x$accuracy$Pathogenic[!is.na(x$accuracy$CADD_raw_rankscore)],
                         chklen=FALSE),
             modnames=c("VARPP","CADD_raw_rankscore"),
             dsids=1:2)
    )


  model_results <-     auc(results_varpp) %>%
    filter(curvetypes %in% "PRC") %>%
    select(-c(dsids,curvetypes)) %>%
    mutate(PP100 = c(sum(x$accuracy[order(-(x$accuracy$rf_results)), ][1:100, "Pathogenic"])/100,
                     sum(x$accuracy[order(-(x$accuracy$CADD_raw_rankscore)), ][1:100, "Pathogenic"])/100)) %>%
    rename(auPRC = aucs)

  model_results
}

#====================================================================================================================

#' Calculate kappa statistic
#' @param cross_table the confusion Matrix of predictions and actual data
#' @import precrec
#' @export
kappa_stats <- function(cross_table){

  diagonal.counts    <- diag(cross_table)
  N                  <- sum(cross_table)
  row.marginal.props <- rowSums(cross_table)/N
  col.marginal.props <- colSums(cross_table)/N

  # Compute kappa (k)
  Po <- sum(diagonal.counts)/N
  Pe <- sum(row.marginal.props*col.marginal.props)
  k <- (Po - Pe)/(1 - Pe)
  k
}

#====================================================================================================================


#' Prediction of single rules
#'
#' @param rulename is the name of one of the rules as returned by the varppRuleFit model
#' @param rulefit_results_object, a varppRuleFit object
#'
#' @export
selected_rule_performance <- function(rulename,
                                      rulefit_results_object
){

  data           = rulefit_results_object$RuleData
  predictions    = rulefit_results_object$RuleFit$accuracy[,1]
  predicted_data = subset(data, data$GeneVariant %in% predictions)

  tab <- table(rule = predicted_data[,rulename], Pathogenic=predicted_data[,"Pathogenic"])

  if(nrow(tab)!=ncol(tab)){

    missings      <- setdiff(colnames(tab),rownames(tab))
    missing_mat   <- mat.or.vec(nr = length(missings), nc = ncol(tab))
    tab           <- as.table(rbind(as.matrix(tab), missing_mat))
    rownames(tab) <- colnames(tab)
  }

  caret::confusionMatrix(tab)

}

#====================================================================================================================


#' Extract rules from ranger trees
#'
#' This function returns rules based on the decision trees built in ranger. It depends on the function varpp
#'
#' @param rf_results the results fro mthe ranger tree generation within the varpp function
#'
#' @return a named vector of rules
#' @import tidypredict
#' @import stringr
#' @export
.extract_ranger_rules <- function(rf_results){

  RULES <- gsub('\n', '',gsub('\"', '', tidypredict::tidypredict_fit(rf_results), ','))
  RULES <- gsub('case_when\\(','',RULES)
  RULES <- gsub('\\)','',RULES)


  rules <- unlist(strsplit(as.character(unlist(RULES)), ','))
  rules <- gsub("\\s+", " ", rules) # This means one or more spaces
  rules <- gsub('~ 1', '~ "1"', rules)
  rules <- gsub('~ 0', '~ "0"', rules)
  rules <- stringr::str_trim(rules)
  rules
}

#====================================================================================================================


#' Simulate data to test the rulefit function
#'
#' @param n the number of instances
#' @param p the number of variables
#' @param seed the seed, to create reproducible results
#'
#' @return a simulated data set with a binary outcome variable and 35 predictors that are associated with the outcome
#
#' @export
simulate_data <- function(n=1000,
                          p = 100,
                          seed = 777){

  n = n
  p = p
  #no influence : 65

  #xij ∼ U(0, 1)
  set.seed(seed)
  x <- NULL
  for( i in 1:p){

    x <- cbind(x,runif(n, min = 0, max = 1))

  }

  x <- data.frame(x)
  names(x) <- paste0("x", 1:dim(x)[2])

  # Create noise variable
  e <- rnorm(n, mean = 0, sd =2)

  y <-  10*(exp(x[,1]^2) * exp(x[,2]^2) * exp(x[,3]^2) * exp(x[,4]^2) * exp(x[,5]^2)) + rowSums(x[,6:35] + e)

  y_bin <- ifelse(y > median(y), 1, 0)


  x$y <- y_bin
  x$index <- c(1:dim(x)[1])
  x <- x[,c("index", "y", names(x)[grep("x", names(x))])]
  x$index.1 <- NULL

  x
}


#====================================================================================================================


#' Density plot for the rule predictions
#'
#' @param rulefit_results an object of class varppRuleFit
#'
#' @return a density plot for the predictions based on the final rules
#
#' @export
density_plot <- function(rulefit_results){

  # Create a data set to show the ratio of predicted classes in the rules versus actual classes
  data.all  <- rulefit_results$RuleData[,c(1,2,3,grep("rule",names(rulefit_results$RuleData)))]
  plot.data <- data.frame(pathogenic = data.all[,3], predictions= rowSums(data.all[,grep("rule", names(data.all))]))

  # Plot the density
  plot.data %>%
    mutate(pathogenic = ifelse(pathogenic %in% 1, "Pathogenic", "Benign")) %>%
    ggplot(aes(x=predictions)) +
    geom_density(aes(x=predictions, y=..scaled.., fill=as.factor(pathogenic)), alpha=1/2) +
    theme_minimal() +
    scale_fill_discrete(name = names(rulefit_results$RuleData)[2], labels = c("negative", "positive")) +
    facet_wrap(~pathogenic)
}


#====================================================================================================================


#' Scatterplot of # of rules that predict correctly versus # of variants per gene
#'
#' @param rulefit_results an object of class varppRuleFit
#' @param y the outcome variable
#'
#' @return a scatterplot of # of rules and # of variants per gene
#
#' @export
ruleVariantPlot <- function(rulefit_results) {

  DATA          <- rulefit_results$RuleData
  TEST          <- data.frame(Gene=DATA$Gene , GeneVariant= DATA$GeneVariant, pathogenic = DATA[,names(rulefit_results$RuleData)[[3]]], predictions= rowSums(DATA[,grep("rule", names(DATA))]))
  variant_count <- data.frame(plyr::count((DATA$Gene)))

  # Merge the variant counts with the count of rules that predict a variant correctly
  dat <- merge(TEST, variant_count, by.x="Gene", by.y="x")


  # Filter out the pathogenic variants only
  dat_patho <- dat %>%
    filter(pathogenic %in% 1) %>%
    rename(Number_of_variants = freq)

  misclassified <- plyr::count(dat_patho[which(dat_patho$predictions %in% 0), "Gene"])
  misclassified <- misclassified %>% rename(Number_of_misclassifications = freq)
  # Misclassified versus number of variants
  test_dat <- merge(dat_patho[,c("Gene", "GeneVariant", "Number_of_variants")], misclassified, by.x="Gene", by.y="x")

  # Plot the #Rules with the number of variants, one variant at a time

  # In this plot, 100% on the y-axis means 100% misclassification (0 instead of 1), larger % on the x-axis
  # means there is a larger percantage of variants in this gene relative to the others
  test_dat %>%
    ggplot(aes(x=(Number_of_variants/sum(unique(Number_of_variants))) *100, y= (Number_of_misclassifications/Number_of_variants)*100, size=Number_of_variants, color=Gene)) +
    geom_point() +
    xlab("% of pathogenic variants per gene (as a percent of all variants)") +
    ylab("% of misclassification ('0 rules get it right')") +
    theme_minimal() +
    guides(color=FALSE) +
    scale_size_area(name= "# of variants per gene") +
    coord_flip()
  #) %>% hide_legend()
}


#' varIMP: Function to extract the variable importance of the expression data variables
#'
#' This function is provided on top of ruleVarImp. It Re-weights the rule kappas by the variabels selected per rule
#' and returns a 0 to 1 scaled importance value per tissue. The most important variabel will have a value of 1.
#' This is based on the variabel importance described in the RuleFit publication by XXX
#'
#' @param rule_model This is the RuleFit model object
#' @param HPOterm  Add the HPO term name for the model
#' @export

varIMP <- function(rule_model=NULL,
                   HPOterm){

  if(is.null(rule_model)){
    stop("Please provide a RuleFit object")
  }

  range01 <- function(x){
    (x-min(x))/(max(x)-min(x))
  }

  rule_kappas_to_merge <- data.frame(rules = names(rule_model$rule_kappas), kappa= as.numeric(rule_model$rule_kappas))
  rules                <- merge(rule_model$RuleFit$varimp, rule_kappas_to_merge, by.x="Variable", by.y="rules")

  DAT <- NULL

  for(i in 1:dim(rules)[1]){
    tissues <- unlist(str_extract_all(as.vector(unlist(rules[grep("rule",rules$Variable), "Description"][i])), "[A-Z_a-z_0-9.]+"))
    DAT <- rbind(DAT, data.frame(tissues, kappa=rules[i,"kappa"], count=length(unique(tissues))))
  }

  DAT %>%
    group_by(tissues) %>%
    summarise(kappa=mean(kappa), count=mean(count)) -> results
  tiss <- plyr::count(unlist(str_extract_all(as.vector(unlist(rules[grep("rule",rules$Variable), "Description"])), "[A-Z_a-z_0-9.]+")))
  tissue_count <- tiss[-grep("^\\d", tiss$x),] %>%
    arrange(desc(freq))

  all_info <- merge(results, tissue_count, by.x="tissues", by.y="x")

  all_info %>%
    mutate(heat_value=kappa/count*freq) %>%
    arrange(desc(heat_value)) %>%
    filter(!tissues %in% "CADD_raw_rankscore") %>%
    rename(!!HPOterm:=heat_value) %>%
    rename(Variables = tissues) -> all_info

  all_info[,5] <- range01(all_info[,5])
  RESULTS <- all_info[,c(1,5)]

  RESULTS <- RESULTS %>%
    replace(is.na(.),0)
  RESULTS
}



#' Function to remove the CADD score variable including '>' , '<' and '=' from the rules
#'
#' @param varppRuleObject the results from varppRule
#'
#' @importFrom stringr str_squish str_remove
#' @export
removeCADD <- function(varppRuleObject){
  topRule <- names(varppRuleObject$rule_kappas[which(varppRuleObject$rule_kappas %in% max(varppRuleObject$rule_kappas))])

  rule <- varppRuleObject$RuleFit$varimp[grep("rule",varppRuleObject$RuleFit$varimp$Variable),] %>%
    filter(Variable %in% topRule) %>%
    select(Description) %>%
    unlist() %>%
    as.character()

  rule_trimmed <- str_remove(str_squish(str_remove(rule, "[&]?\\s?CADD_raw_rankscore [<|>|<=|>=]* [0-9]*\\.?[0-9]*\\s?")),
                             "[&]?\\s?CADD_raw_rankscore [<|>|<=|>=]* [0-9]*\\.?[0-9]*\\s?")
  str_squish(str_remove(rule_trimmed, "^&\\s?"))
}



#' Function to extract CADD score including cut-off
#'
#' @param varppRuleObject the results from varppRule
#'
#' @importFrom stringr str_squish str_extract
#' @export
getCADDcutOff <- function(varppRuleObject) {
  topRule <- names(varppRuleObject$rule_kappas[which(varppRuleObject$rule_kappas %in% max(varppRuleObject$rule_kappas))])

  rule <- varppRuleObject$RuleFit$varimp[grep("rule",varppRuleObject$RuleFit$varimp$Variable),] %>%
    filter(Variable %in% topRule) %>%
    select(Description) %>%
    unlist() %>%
    as.character()

  rule_trimmed <- paste(str_extract_all(rule, "[&]?\\s?CADD_raw_rankscore [<|>|<=|>=]* [0-9]*\\.?[0-9]*\\s?") %>% unlist() %>% str_replace("& ", ""), collapse=", ")
  str_squish(str_remove(rule_trimmed, "^&\\s?"))
}



#' Return a gene panel based on the tissues in the top rule
#'
#' @param varppRuleFitObject the results from varppRule
#'
#' @importFrom stringr str_squish str_remove
#' @importFrom magrittr '%>%'
#' @import dplyr
#' @export

genePanelTop <- function(varppRuleObject) {

  PANEL <- varppRuleObject$RuleData[,c("Gene", "Pathogenic")]
  topRule <- names(varppRuleObject$rule_kappas[which(varppRuleObject$rule_kappas %in% max(varppRuleObject$rule_kappas))])

  rule <- varppRuleObject$RuleFit$varimp[grep("rule",varppRuleObject$RuleFit$varimp$Variable),] %>%
    filter(Variable %in% topRule) %>%
    select(Description) %>%
    unlist() %>%
    as.character()


    if(grepl('CADD_raw_rankscore >', rule)) {
      PANEL <-  cbind(PANEL,varppRuleObject$RuleData %>%
                        mutate(!!rule := ifelse(eval(parse(text=removeCADD(varppRuleObject))), 1, 0 )) %>%
                        select(!!rule) )

    }else{

      PANEL <- cbind(PANEL,varppRuleObject$RuleData %>%
                       mutate(!!rule := ifelse(eval(parse(text=removeCADD(varppRuleObject))), 0, 1 )) %>%
                       select(!!rule))


    }

  varppRuleObject$RuleData %>%
    mutate(panel_selection = ifelse(PANEL[,3] %in% 1, TRUE, FALSE)) %>%
    filter(panel_selection %in% TRUE) %>%
    select(Gene) %>%
    unlist() %>%
    as.character() %>%
    unique()

}


#===========#
#           #
#  /(_M_)\  #
# |       | #
#  \/~V~\/  #
#           #
#===========#
