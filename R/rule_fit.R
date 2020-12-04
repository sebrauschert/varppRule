#' The RuleFit function
#'
#' This is rulefit in general. When not using the VARPP framework, this function allows to perform rulefit on any data set with a
#' binary outcome.
#'
#' @param data a data set with features, an outcome variable and an index variable (probably created within the function in the end)
#' @param y the outcome variable name
#' @param two_level_bootstrap this is TRUE by default and is necessar when sampling has to be performed on agene and gene variant level
#' @param ntree number of trees to be built, defaults to 200
#' @param max.depth maximum tree depth, defaults to 3
#' @param rule.filter filter the top n rules based on kappa statistic. If NULL, the rules are filter above a kappa of 0.05
#' @param bootstrap.rounds number of bootstrap rounds for the outer loop of the LASSO cross-validation, defaults to 100
#' @param rule.extract.cores number of cores for parallel, defaults to 4. This is specifically for the varpp rule extract step (less memory hungry than the cv.glmnet step)
#' @param kappa.cores number of cores used for the rule filtering by kappa. This needs to be seperate, as it is quite memory intensive when the input + rule data is very large. Defaults to 2.
#' @param lasso.cores number of cores for the cv.glmnet step, as this is quite memory hungry, it is seperated
#' @param ranger.and.rulefit in some cases, if not all, we want to be able to report the Random forest and the RuleFit accuracy together, to see how much better one performs over the other.
#' This defaults to TRUE
#' @return A list of predictions for the outcome. Further, a variable importance list for all rules and variables tested.
#'
#' @import doMC
#' @import doParallel
#' @import parallel
#' @import ranger
#' @import foreach
#' @import caret
#' @import DMwR
#' @import glmnet
#' @import tidyverse
#' @import magrittr
#' @importFrom parallel detectCores
#' @importFrom caret trainControl train nearZeroVar
#' @import DMwR
#' @import lattice
#' @import ggplot2
#' @import grid
#' @import precrec
#' @import DT
#' @import rmarkdown
#' @import knitr
#' @import pander
#' @import tidyr
#' @import tidyselect
#' @export
rule_fit <- function(HPO,
                     ntree = 200,
                     max.depth = 3,
                     rule.filter = 10,
                     bootstrap.rounds = 100,
                     rule.extract.cores = 4,
                     kappa.cores = 2,
                     lasso.cores = 4,
                     ranger.and.rulefit = FALSE
                    ){


  #========================================================================================
  # Setup
  #========================================================================================
  # Record the call and return it with the results later
  CALL = match.call()

  # Time measurment
  start_time <- Sys.time()

  # # Function checklist and stop commands for missing values
  if(is.null(HPO))
    stop("HPO genes list not provided")

  #========================================================================================
  # Prepare the data
  #========================================================================================
  #patho  <- varppRule:::patho
  #benign <- varppRule:::benign

  # Get the gene names
  hpo_gene_names <- HPO$Gene

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
  data <- list(dat=data.frame(rbind(varpp_patho, varpp_benign)), disease_variants=data.frame(varpp_patho), benign_variants=data.frame(varpp_benign))



  #========================================================================================
  # Step I: Random Forest and rule extraction
  #========================================================================================

  if(ranger.and.rulefit == TRUE){

    message("Initiate varpp with ranger and Rulefit results...")
    message(paste0("Extracting rules from ",ntree, " trees, with a max.depth of ", max.depth, ", using ", rule.extract.cores, " cores"))


      # Bootstrap rounds is a new addition as I try to save the dat_boot oject for
      # the lasso step; the lasso step currently does the exact same sampling as we do here, which should not be necessary.
      # There is, however, an additional element of randomness in the model when doing the sampling again later. Maybe revise?
      varpp_plus_rules <- varpp(data,
                                     ntree = ntree,
                                     max.depth = max.depth,
                                     cores = rule.extract.cores,
                                     bootstrap.rounds = bootstrap.rounds)


      rules            <- varpp_plus_rules$rules
      rf_varpp_results <- list(accuracy = varpp_plus_rules$accuracy, varimp =  varpp_plus_rules$varimp)

      # remove that as it might stick in memory
      rm(varpp_plus_rules)
      gc()



  }else{
    message(paste0("Extracting rules from ",ntree, " trees, with a max.depth of ", max.depth, ", using ", rule.extract.cores, " cores"))


      # Bootstrap rounds is a new addition as I try to save the dat_boot oject for
      # the lasso step; the lasso step currently does the exact same sampling as we do here, which should not be necessary.
      # There is, however, an additional element of randomness in the model when doing the sampling again later. Maybe revise?
      # VARPP and extract the rules
      varpp_start_time   = Sys.time()
      rules <- varpp_for_rulefit(dat=data,
                     ntree = ntree,
                     max.depth = max.depth,
                     cores = rule.extract.cores)

      varpp_end_time   = Sys.time()
      varpp_time_taken = varpp_end_time - varpp_start_time
      message(paste0("varpp finished after ", paste(round(varpp_time_taken,2), units(varpp_time_taken), sep=" ")))

    }


  #========================================================================================
  # CREATING RULE VARIABLES
  #========================================================================================

  rules_data = list()

  message("Creating rule variables...")
  rule_start_time <- Sys.time()

  # Create the values that the rules take (0 or 1, this is part of the extracted "rule sentence")
  positive    = as.vector(unlist(lapply(rules, function(x) if(grepl('\"1\"', x) %in% TRUE) return(1) else(return(0)) )))
  alternative = as.vector(unlist(lapply(rules, function(x) if(grepl('\"1\"', x) %in% TRUE) return(0) else(return(1)) )))

  # remove the ~ 0 and  ~1 parts of the rules to make them evaluatable
  rules = unlist(lapply(rules, function(x) str_split(x, " ~")[[1]][1]))


    cl <- makeCluster((rule.extract.cores/2))
    registerDoParallel(cl)

    rules_data <- foreach( i = 1:length(rules)) %dopar% {
      data$dat %>%
        mutate(!!rules[i] := ifelse(eval(parse(text=rules[i])), positive[i], alternative[i])) %>%
        select(!!rules[i])

    }
    stopCluster(cl)


  rule_end_time   = Sys.time()
  rule_time_taken = rule_end_time - rule_start_time


  message(paste0("Rules created after ", paste(round(rule_time_taken,2), units(rule_time_taken), sep=" ")))

  # This will combine the rules created in the above loop
  rule_dat = do.call(cbind, rules_data)
  rm(rules_data, positive, alternative)
  gc()

  # Remove duplicated columns
  rule_dat  <- rule_dat[!duplicated(as.list(rule_dat))]

  # Save the total number of Rules
  total_rules_without_duplicates = dim(rule_dat)[2]

  #===================================================================================
  # FILTERING RULES: ONLY 0, ONLY 1 and KAPPA filter
  #===================================================================================
  message(paste0("Rules before filtering: ", total_rules_without_duplicates))
  message("Filtering rules...")

  # Maybe do the conditional zero variance thing here?
  # Remove Rules that are only 0 or only 1 (based on the condition if all tested rules return FALSE, don't remove as ther is nothing to remove )
  if(length(which((as.numeric(unlist(lapply(rule_dat, function(x) var(x, na.rm=T)))) %in% 0) %in% FALSE)) < 1){

  rule_dat <- rule_dat[,-which(as.numeric(unlist(lapply(rule_dat, function(x) var(x)))) %in% 0)]

  }

    patho = data$dat[,y]


  if(is.null(rule.filter)){

    rule_kappas <- unlist(parallel::mclapply(rule_dat, function(x){
      kappa_stats(table(x, patho))
    }, mc.cores = kappa.cores, mc.allow.recursive = TRUE))

    # Filtering based on KAPPA
    rule_kappas <- sort(rule_kappas)

  message(paste0("All ",dim(rule_dat)[2], " rules are used for LASSO..."))

  }else{
  if(rule.filter%%1 == 0){
      # Remove all rules with a Kappa value lower than 0.1
      rule_kappas <- unlist(parallel::mclapply(rule_dat, function(x){
        kappa_stats(table(x, patho))
        }, mc.cores = kappa.cores, mc.allow.recursive = TRUE))

      rule_kappas <- sort(rule_kappas)[(length(rule_kappas)-rule.filter+1):length(rule_kappas)]
      rule_dat <- rule_dat[,names(rule_kappas)]
      message(paste0(dim(rule_dat)[2], " rules are used for LASSO..."))
  }

      if(!rule.filter%%1 == 0){

        # Remove all rules with a Kappa value lower than rule.filter
        rule_kappas <- unlist(parallel::mclapply(rule_dat, function(x){
          kappa_stats(table(x, patho))
          }, mc.cores = kappa.cores, mc.allow.recursive = TRUE))

        rule_dat <- rule_dat[,-which(rule_kappas < rule.filter)]
        message(paste0(dim(rule_dat)[2], " rules are used for LASSO..."))

      }

    }

  rm(patho)
  gc()

  # Save the total number of Rules
  total_rules_after_var_removal = dim(rule_dat)[2]


  if(total_rules_after_var_removal == 0)
    stop("No rules had a kappa value > 0.05. Modelling stopped")

  data$dat        <- cbind(data$dat, rule_dat)
  data_dimensions <- dim(data$dat)

  # Remove unnecessary data
  rm(rule_dat)
  gc()

  #========================================================================================
  # Step II: glmnet
  #========================================================================================

  message("Starting LASSO...")

    LASSO <- lasso_ensemble(data,
                              rules,
                              bootstrap.rounds=bootstrap.rounds,
                              cores = lasso.cores)



  LASSO$time          <- NULL
  LASSO$memory        <- NULL
  LASSO$process_stage <- NULL

  # Calculate Kappa statistic for every rule and sort the rule slot by that
  # SELECTED_RULES, the subset of the LASSO results with "rule" in the name:
  # as.character(LASSO$varimp$Variable[grep("rule",LASSO$varimp$Variable)])

  if(!length(as.character(LASSO$varimp$Variable[grep("rule",LASSO$varimp$Variable)])) == 0){


      predicted_data = subset(data, data$GeneVariant %in% LASSO$accuracy$GeneVariant)

      LASSO$varimp <- LASSO$varimp %>%
        filter(!CADD_expression %in% NaN)

      LASSO$varimp <- data.frame(Variable=LASSO$varimp$Variable, Coefficient= LASSO$varimp$CADD_expression, Description = LASSO$varimp$Rules) %>%
      arrange(desc(abs(Coefficient)))
  }else{

    message("No rules identified to keep after LASSO.")

  }

  message("LASSO finished! Cleaning up and preparing results.")

  end_time   = Sys.time()
  time_taken = end_time - start_time


          results        <- list(RuleFit = LASSO,
                                 y=y,
                                   RuleData = data$dat,
                                   duration = time_taken,
                                   ntree = ntree,
                                   maxdepth = max.depth,
                                   ranger.and.rulefit = ranger.and.rulefit,
                                   rulefit_call = CALL,
                                   total_rules_without_duplicates = total_rules_without_duplicates,
                                   total_rules_after_var_removal = total_rules_after_var_removal,
                                   data_dimensions = data_dimensions,
                                   rule_kappas = rule_kappas,
                                   bootstrap_rounds = bootstrap.rounds)


  class(results) <- "varppRule"


  #========================================================================================
  # Return the results
  #========================================================================================

  results

}
