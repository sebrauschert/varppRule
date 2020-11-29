# Package varpp RuleFit

Human interpretable predictive model for variant prioritisation of rare and other genetic diseases

## Installation

To install, type the following into `R` or `RStudio`:

```r
devtools::install_github("Hobbeist/varppRuleFit")
```

## Run

```r
# Load the package
library(varppRuleFit)

# Specify the file locations
benign_variant_file     <- "/data/benign_variants_simplified_51151.txt"
pathogenic_variant_file <- "/data/vp_gene_2424/HP:0000591"
expression_file         <- "/data/GTExV8_specificity_simplified.csv"

#---------------------------------------------------------------------------------

# Run the function

rulefit_results <- rule_fit(expression_file = expression_file, 
                            benign_variant_file = benign_variant_file,
                            pathogenic_variant_file = pathogenic_variant_file,
                            two_level_bootstrap = TRUE,
                            ntree = 500,
                            max.depth = 3,
                            bootstrap.rounds = 100,
                            smote = FALSE,
                            expression.source = "gtex",
                            rule.extract.cores = 32,
                            lasso.cores = 2,
                            ranger.and.rulefit = TRUE,
                            report = TRUE,
                            report_filename="HP0000591",
                            HPO_term_name = "HP:0000591")
                            )

```



<p align="center">

<img src="varpp-logo.png"  width="40%" height="40%">

</p>