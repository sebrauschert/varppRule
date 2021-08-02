#' The RuleFit function
#'
#' RuleFit creates variant predictions and human interpretable rules
#'
#'
#' @param HPO_genes HPO term associated list of genes, or any list of patient genes.
#' @param HPO_term_name In case the model is for one specific HPO term, this can be provided, otherwise it is assigned as "custom"
#' @param type the prediction data; either hcl (single cell), gtex (tissue specific) or custom (requires the user to provide custom_patho and custom_benign).
#' @param user_patho a user provided file for the pathogenic variants. This needs to have the following first few columns:Gene, GeneVariant, CADD_raw_rankscore, CADD_PHRED_SCORE, Pathogenic, gene_id,gene_biotype
#' @param user_benign a user provided file for the benign variants. This needs to have the following first few columns:Gene, GeneVariant, CADD_raw_rankscore, CADD_PHRED_SCORE, Pathogenic, gene_id,gene_biotype
#' @param ntree number of trees to be built, defaults to 200.
#' @param max.depth maximum tree depth, defaults to 3.
#' @param rule.filter filter the top n rules based on kappa statistic. If NULL, the rules are filter above a kappa of 0.05.
#' @param bootstrap.rounds number of bootstrap rounds for the outer loop of the LASSO cross-validation, defaults to 100.
#' @param rule.extract.cores number of cores for parallel, defaults to 4. This is specifically for the varpp rule extract step (less memory hungry than the cv.glmnet step).
#' @param kappa.cores number of cores used for the rule filtering by kappa. This needs to be separate, as it is quite memory intensive when the input + rule data is very large. Defaults to 2.
#' @param lasso.cores number of cores for the cv.glmnet step, as this is quite memory hungry, it is separated.
#' @return A list of predictions for the outcome. Further, a variable importance list for all rules and variables tested.
#'
#' @import doMC
#' @import doParallel
#' @import parallel
#' @import ranger
#' @import foreach
#' @import caret
#' @import glmnet
#' @import tidyverse
#' @import magrittr
#' @importFrom caret trainControl train nearZeroVar
#' @import lattice
#' @import ggplot2
#' @import grid
#' @import precrec
#' @import DT
#' @import tidyr
#' @import tidyselect
#' @export
rule_fit <- function(HPO_genes,
                     HPO_term_name = "custom",
                     type = c("gtex", "hcl", "custom"),
                     user_patho = NULL,
                     user_benign = NULL,
                     ntree = 200,
                     max.depth = 3,
                     rule.filter = 10,
                     bootstrap.rounds = 100,
                     rule.extract.cores = 4,
                     kappa.cores = 2,
                     lasso.cores = 4
){


  #========================================================================================
  # Setup
  #========================================================================================
  # Record the call and return it with the results later
  CALL = match.call()

  # Time measurement
  start_time <- Sys.time()

  # # Function checklist and stop commands for missing values
  if(is.null(HPO_genes))
    stop("HPO genes list not provided")

  if(is.null(ntree))
    stop("Please specify the number of trees for the model")

  if(is.null(max.depth))
    stop("Please specify the max.depth for the model")

  if(is.null(type))
    stop("Please specify the data type as either 'gtex' or 'hcl'")


  #========================================================================================
  # Prepare the data: based on gtex or hcl, the respective internal data is chosen
  #========================================================================================


  if(type %in% "gtex"){

    patho  <- patho_gtex
    benign <- benign_gtex

  } else if (type %in% "hcl"){
    patho  <- patho_hcl
    benign <- benign_hcl

  } else {
    patho = user_patho
    benign = user_benign
  }



  # Specify Pathogenic as the outcome
  y <- "Pathogenic"

  # Filter the genes based on the input list
  patho %>%
    filter(Gene %in% HPO_genes) %>%
    select(-c(CADD_raw_rankscore)) %>%
    rename(CADD_raw_rankscore = CADD_PHRED_SCORE) -> varpp_patho

  # As the packages contains the hcl and gtex data, the column names are in a standardized
  # order so the below arranges it correctly
  varpp_patho <- varpp_patho[,c(1,2,4,3,7:length(names(varpp_patho)))]


  # Filter out the benign genes that are in the pathogenic gene list provided by the user
  benign %>%
    filter(!Gene %in% intersect(benign$Gene, patho$Gene)) %>%
    select(-c(CADD_raw_rankscore)) %>%
    rename(CADD_raw_rankscore = CADD_PHRED_SCORE) -> varpp_benign

  # As the packages contains the hcl and gtex data, the column names are in a standardized
  # order so the below arranges it correctly
  varpp_benign <- varpp_benign[,c(1,2,4,3,7:length(names(varpp_benign)))]

  # Create the input data for varppRule
  data <- list(dat = data.frame(rbind(varpp_patho, varpp_benign)),
               disease_variants = data.frame(varpp_patho),
               benign_variants = data.frame(varpp_benign))



  #========================================================================================
  # Step I: Random Forest and rule extraction
  #========================================================================================
  message(paste0("Extracting rules from ",
                 ntree,
                 " trees, with a max.depth of ",
                 max.depth,
                 ", using ",
                 rule.extract.cores,
                 " cores and ",
                 type,
                 " data"))


  # VARPP and extract the rules: Run the custom Random Forest via VARPP and record the rules
  varpp_start_time  <- Sys.time()
  rules <- varpp_for_rulefit(dat = data,
                             ntree = ntree,
                             max.depth = max.depth,
                             cores = rule.extract.cores)

  varpp_end_time   = Sys.time()
  varpp_time_taken = varpp_end_time - varpp_start_time

  message(paste0("varpp finished after ",
                 paste(round(varpp_time_taken,2),
                       units(varpp_time_taken), sep=" ")))


  #========================================================================================
  # CREATING RULE VARIABLES
  #========================================================================================

  rules_data = list()

  message("Creating rule variables...")
  rule_start_time <- Sys.time()

  # Create the values that the rules take (0 or 1, this is part of the extracted "rule sentence")
  # We assign a 1 if the rule contains the assignment to 1, and for the alternative (if the rule is FALSE)
  # we record a 0. Hence we have two complimentary vectors for all rules. We save this in the final ruleFit object,
  # so we ca npredict with the ruls and assign the correct outcome.

  # This selects only those rules that predict the outcome
  rules <- rules[which(grepl('\"1\"', rules))]

  positive    <- as.vector(unlist(lapply(rules, function(x) if(grepl('\"1\"', x) %in% TRUE) return(1) else(return(0)) )))
  alternative <- as.vector(unlist(lapply(rules, function(x) if(grepl('\"1\"', x) %in% TRUE) return(0) else(return(1)) )))

  # remove the ~ 0 and  ~1 parts of the rules to make them evaluatable;
  # this is the reason why we also need to record the 'positive' and 'alternative'
  # information above, as otherwise it is not clear what prediction we make.
  #
  # TO DO: consider only taking forward rules that make a positive prediction, e.g. predict
  #        the pathogenic variants?
  #
  rules <- unlist(lapply(rules, function(x) str_split(x, " ~")[[1]][1]))


  # We parallelise the rule creation step to speed it up
  cl <- makeCluster((rule.extract.cores/2))
  registerDoParallel(cl)

  rules_data <- foreach( i = 1:length(rules)) %dopar% {
      data$dat %>%
        mutate(!!rules[i] := ifelse(eval(parse(text=rules[i])), positive[i], alternative[i])) %>%
        select(!!rules[i])

    }
  stopCluster(cl)

  # Name the positive and alternative variables like the rules, so they can be applied in the predict function
  names(positive)    <- names(rules)
  names(alternative) <- names(rules)

  rule_end_time   <- Sys.time()
  rule_time_taken <- rule_end_time - rule_start_time


  message(paste0("Rules created after ",
                 paste(round(rule_time_taken,2),
                       units(rule_time_taken), sep=" ")))

  # This will combine the rules created in the above loop
  rule_dat <- do.call(cbind, rules_data)
  rm(rules_data)
  gc()



  #===================================================================================
  # FILTERING RULES: ONLY 0, ONLY 1 and KAPPA filter
  #===================================================================================

   # Remove duplicated columns first: This means, if rules lead to the exact same prediction,
  # we remove them. This dose, however, not exclude complimentary rules
  #
  # TO DO: Implement a filter that excludes complimentary rules to remove redundant
  #        Information
  #
  # This is the solution to exclude both duplicated and complementary rules, as we exclude
  # correlated rules. I believe there is really no case to be made for keeping completely correlated rules
  # or complementary rules in the model. We want the model to be as "lean" as possible to perform at its
  # best.

  split_1     <- rule_dat[, -grep("rule", names(rule_dat))]
  split_2     <- rule_dat[, grep("rule", names(rule_dat))]
  split_2     <- split_2[, -caret::findCorrelation(cor(split_2))] # maybe just copy the code from caret?
  rule_dat    <- cbind(split_1, split_2)

  # This is the old function to simply remove duplicated rules.
  #rule_dat  <- rule_dat[!duplicated(as.list(rule_dat))]

  # Save the total number of Rules
  total_rules_without_duplicates <- dim(rule_dat)[2]

  message(paste0("Rules before filtering: ", total_rules_without_duplicates))
  message("Filtering rules...")

  # Maybe do the conditional zero variance thing here?
  # Remove Rules that are only 0 or only 1 (based on the condition if all tested rules return FALSE, don't remove as there is nothing to remove )
  if(length(which((as.numeric(unlist(lapply(rule_dat, function(x) var(x, na.rm=T)))) %in% 0) %in% FALSE)) < 1){

  rule_dat <- rule_dat[,-which(as.numeric(unlist(lapply(rule_dat, function(x) var(x)))) %in% 0)]

  }

  patho <- data$dat[,y]


  # The below is the kappa filter for the rules:
  # depending on the user input, this will be either based on the number of rules (10 by default)
  # or a specific kappa cut-off
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
      rule_dat    <- rule_dat[,names(rule_kappas)]
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
  total_rules_after_var_removal <- dim(rule_dat)[2]


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

  # After creatign the rules and filtering them based on the kappa statistic,
  # we can now use the newly created variables together with the original variables
  # and run the LASSO model
  message("Starting LASSO...")

    LASSO <- lasso_ensemble(data,
                              rules,
                              bootstrap.rounds=bootstrap.rounds,
                              cores = lasso.cores)



  LASSO$time          <- NULL
  LASSO$memory        <- NULL
  LASSO$process_stage <- NULL

  # Calculate Kappa statistic for every rule and sort the rule slot by that
  # SELECTED_RULES, the subset of the LASSO results with "rule" in the name

  if(!length(as.character(LASSO$varimp$Variable[grep("rule",LASSO$varimp$Variable)])) == 0){


      predicted_data <- subset(data, data$GeneVariant %in% LASSO$accuracy$GeneVariant)

      LASSO$varimp <- LASSO$varimp %>%
        filter(!CADD_expression %in% NaN)

      LASSO$varimp <- data.frame(Variable=LASSO$varimp$Variable,
                                 Coefficient= LASSO$varimp$CADD_expression,
                                 Description = LASSO$varimp$Rules) %>%
      arrange(desc(abs(Coefficient)))
  }else{

    message("No rules identified to keep after LASSO.")

  }


  # Save the rules for the prediction. This is necessary to keep track of the
  # rules and their 'positive' and 'alternative' prediction

  rules_to_apply <- names(data)[grep("rule", names(data))]
  rules_to_apply <- rules[rules_to_apply]


  message("LASSO finished! Cleaning up and preparing results.")

  end_time   <- Sys.time()
  time_taken <- end_time - start_time

  # Preparing the final ruleFit output object.

                         # The LASSO results
  results        <- list(RuleFit = LASSO,
                         # The outcome variable name, which will be 'Pathogenic'
                         y = y,
                         # The data used as input for the LASSO step, containing the rule variables
                         RuleData = data$dat,
                         # The names list of positive predictions per rule
                         positive = positive,
                         # The named list of alternative predictions per rule
                         alternative = alternative,
                         # The vector of rules that are actually used in the ruleFit model
                         rules_to_apply = rules_to_apply,
                         duration = time_taken,
                         ntree = ntree,
                         maxdepth = max.depth,
                         rulefit_call = CALL,
                         total_rules_without_duplicates = total_rules_without_duplicates,
                         total_rules_after_var_removal = total_rules_after_var_removal,
                         data_dimensions = data_dimensions,
                         rule_kappas = rule_kappas,
                         bootstrap_rounds = bootstrap.rounds,
                         HPO_term_name = HPO_term_name)


  class(results) <- "varppRule"


  #========================================================================================
  # Return the results
  #========================================================================================

  results

}





#=============#
#             #
#  /(__M__)\  #
# |         | #
#  \/~-V~-\/  #
#             #
#=============#

