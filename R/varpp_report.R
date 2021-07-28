#' Report function for varpp results
#'
#' This function will create a report based on the results from the rulef_fit() function.'
#'
#' @param results the results from the rule_fit() function
#' @param report_filename The path, including the filename for the resulting report.
#'
#' @import DT
#' @import ggplot2
#' @import tidyverse
#' @import grid
#' @import rmarkdown
#' @import knitr
#' @import pander
#' @import tidyr
#' @import tidyselect
#' @export
varpp_report <- function(results,
                         report_filename) {

    gene_var       <- names(results$RuleFit$accuracy)[1]

    # Specify the outcome variable
    y = results$y

    table_function <- performance

    # Pathogenic Variants: Specify for the score to be above 0.5, which means more than 50% of the predictions classified the variant as pathogenic.
    all_res <- results$RuleFit$accuracy[,c(gene_var,y, "CADD_expression", "CADD_raw_rankscore")]
    all_res <- all_res %>%
        rename(RuleFitScore = CADD_expression) %>%
        filter(RuleFitScore > .threshold(results$RuleFit$accuracy[,2], results$RuleFit$accuracy[,3]))


    # Create a data set to show the ratio of predicted classes in the rules versus actual classes
    DATA_all <- results$RuleData[,c(1,2,3,grep("rule",names(results$RuleData)))]


    # Create a data table element with all the pathogenic Variants, that neither of the rules was able to classify as pathogenic
    DATA_patho <- DATA_all %>%
      filter(get(y) %in% 1)

    TEST <- data.frame(pathogenic = DATA_all[,y], predictions= rowSums(DATA_all[,grep("rule", names(DATA_all))]))



    #========================================================================================
    # Set up the list of parameters to be forwarded to the report
    #========================================================================================
    params <- list(

      all_res = datatable(all_res, style='bootstrap4', extensions = 'Buttons', options = list(
        dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))),
        auprcpp100 = datatable(table_function(results), style='bootstrap4',extensions = 'Buttons', options = list(
        dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))),

      # Rule Results
      rules      = datatable(results$RuleFit$varimp , style='bootstrap4', extensions = 'Buttons', options = list(
        dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))),

      # # Model Info
      ntree      = results$ntree,
      maxdepth   = results$maxdepth,
      HPO_term_name    = results$HPO_term_name,
      rules_before_filtering = results$total_rules_without_duplicates,
      rules_after_filtering  = results$total_rules_after_var_removal,


      # Result tables and plots
      RuleFit_results = metrics(results),
      RuleFit_roc     = auROC(results, model = "rulefit"),
      RuleFit_prc     = auPRC(results, model = "rulefit"),


      density_plot =
        TEST %>%
        ggplot(aes(x = predictions)) +
        geom_density(aes(x=predictions, y=..scaled.., fill=as.factor(pathogenic)), alpha=1/2) +
        theme_minimal() +
        scale_fill_discrete(name = y, labels = c("negative", "positive")) +
        facet_wrap(~pathogenic),

      patho_ratio =
        DATA_all %>%
        tidyr::gather(key = "rule", value, 4:12) %>%
        mutate(value=as.factor(value)) %>%
        group_by(!!! rlang::syms(y), rule, value) %>%
        summarise(n = n()) %>%
        ggplot(., aes(x = rule, y = n, fill = value)) +
        geom_col( position = "fill") +
        ggtitle("Predictions versus actual: Rules") +
        ylab("Class counts") +
        xlab("") +
        scale_fill_discrete(name = y, labels = c("negative hit", "positive hit")) +
        theme_minimal() +
        facet_wrap(~get(y)) +
        theme(axis.text.x = element_text(angle = 45)),

      missclass_patho = datatable(data.frame(Misclassified_positives = DATA_patho[which(rowSums(DATA_patho[,grep("rule",names(DATA_patho))]) == 0), 1]),
                                  style='bootstrap4', extensions = 'Buttons', options = list(
        dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))),

      # Kappa plot for all rules
      kappa_plot = results$rule_kappas,
      DATA_all = DATA_all

    )

    if(is.null(report_filename)){
      report_filename <- getwd()
    }


    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    tempReportRmd <- system.file("rmd", "rulefit_report.Rmd", package = "varppRule")
    suppressMessages(rmarkdown::render(input         = tempReportRmd,
                      output_format = "html_document",
                      output_file   = paste0(report_filename,"_ruleFit_report.html"),
                      clean         = TRUE,
                      params        = params,
                      envir         = new.env(parent = globalenv()),
                      quiet         = TRUE))
}

#===========#
#           #
#  /(_M_)\  #
# |       | #
#  \/~V~\/  #
#           #
#===========#
