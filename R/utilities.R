#' Class specific functions
#'
#' @param x, an object of class varppRuleFit
#' @export
print.varppRuleFit <- function(x){

  if(!class(x) =="varppRuleFit"){
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
  print(head(x$RuleFit$varimp,n = 10), print.gap = 2, quote = FALSE, row.names = FALSE)
  invisible(x$RuleFit$varimp)

  cat("\n")

  cat("=========================================\n")

  cat("\n")
  cat("Model call: \n")
  print(x$rulefit_call)
  cat("\n")
}


#====================================================================================================================


#' Sampling and subsettign for the benign variant data
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

#' Return a table with model names, auPRC, PP100 and ntree of the model: only for regular rulefit, not two way
#' @param x an object of class \code{varppRuleFit}
#' @import precrec
#' @export
auPRC_table <- function(x, ntree=x$ntree){
  results_varpp <-
    evalmod(
      mmdata(scores=join_scores(x$RuleFit$accuracy[,3][!is.na(x$RuleFit$accuracy[,3])],
                                x$RandomForest$accuracy[,3][!is.na(x$RandomForest$accuracy[,3])],
                                chklen=FALSE),
             join_labels(x$RuleFit$accuracy[,2][!is.na(x$RuleFit$accuracy[,3])],
                         x$RandomForest$accuracy[,2][!is.na(x$RandomForest$accuracy[,3])],
                         chklen=FALSE),
             modnames=c("RuleFit", "RandomForest"),
             dsids=1:2)
    )


  model_results <-     auc(results_varpp) %>%
    filter(curvetypes %in% "PRC") %>%
    select(-c(dsids,curvetypes)) %>%
    mutate(PP100 = c(sum(as.numeric(as.character(x$RuleFit$accuracy[order(-(x$RuleFit$accuracy[,3])), ][1:100, 2])))/100,
                     sum(as.numeric(as.character(x$RandomForest$accuracy[order(-(x$RandomForest$accuracy[,3])), ][1:100,2])))/100)) %>%
    mutate(ntree = ntree) %>%
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

  data           = rulefit_results_object$RuleData$dat
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


#' Density plot for the rule predictions
#'
#' @param rulefit_results an object of class varppRuleFit
#'
#' @return a density plot for the predictions based on the final rules
#
#' @export
density_plot <- function(rulefit_results){

  # Create a data set to show the ratio of predicted classes in the rules versus actual classes
  data.all  <- rulefit_results$RuleData$dat[,c(1,2,grep("rule",names(rulefit_results$RuleData$dat)))]
  plot.data <- data.frame(pathogenic = data.all[,2], predictions= rowSums(data.all[,grep("rule", names(data.all))]))

  # Plot the density
  plot.data %>%
    ggplot(aes(x=predictions)) +
    geom_density(aes(x=predictions, y=..scaled.., fill=as.factor(pathogenic)), alpha=1/2) +
    theme_minimal() +
    scale_fill_discrete(name = names(rulefit_results$RuleData$dat)[2], labels = c("negative", "positive")) +
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
