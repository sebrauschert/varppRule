# Package varppRule

Human interpretable predictive model for variant prioritisation of rare and other genetic diseases

<p align="center">

<img src="varpp-logo.png"  width="20%" height="20%">

</p>


## Installation

To install, type the following into `R` or `RStudio`:

```r
devtools::install_github("Hobbeist/varppRule")
```

## Human Phenotype Ontology and related genes  
In order to run RuleFit and VARPP and get a prediction for pathogenic variants, one needs to retrieve a list of genes that are associated with the HPO term of interest. Just provide a vector of genes to the model and it will return the prediction.

The list of genes can either be retrieved manually from https://hpo.jax.org/, or by using the `phenolyzer` software.


## Run RuleFit

```r
# Load the package
library(varppRule)

#---------------------------------------------------------------------------------

# Run the function

rulefit_results <- rule_fit(HPO_genes,
                            ntree = 200,
                            max.depth = 3,
                            rule.filter = 10,
                            bootstrap.rounds = 100,
                            rule.extract.cores = 4,
                            kappa.cores = 2,
                            lasso.cores = 4,
                            ranger.and.rulefit = FALSE
                            )

```

## Run VARPP
```r

varpp_results <- varpp(HPO_genes,
                       ntree = 2000,
                       max.depth = NULL,
                       cores = 4
                      )
                  
```

## Results

The below is an exemplary output from ```varppRule``` for HP:0000002.
Using the ```varpp_report``` function allows to generate an html report based on this output,
which includes more detailed information on the model statistics.

```
=========================================


  _____       _      ______ _ _
 |  __ \     | |    |  ____(_) |
 | |__) |   _| | ___| |__   _| |_
 |  _  / | | | |/ _ \  __| | | __|
 | | \ \ |_| | |  __/ |    | | |_
 |_|  \_\__,_|_|\___|_|    |_|\__|



========== Final RuleFit model ==========

  Cohen's Kappa:  0.7696038 
  F1-score: 0.9546579 
  Number of rules : 10 
  Top Rule : rule_4628 
  Prediction finished after : 14.85761 mins 

=========================================

Top predictors with performance:

    Variable  Coefficient   Description
   rule_4628  -2.21524974   CADD_raw_rankscore >= 17.77 & Pituitary >= 0.00078 & Colon_Sigmoid >= 0.002
   rule_4627   1.64784842   CADD_raw_rankscore < 17.77 & Pituitary >= 0.0007 & Colon_Sigmoid >= 0.002 
   rule_4359   1.04252465   Esophagus_Muscularis >= 0.0019 & Heart_Atrial_Appendage >= 0.0006 & CADD_raw_rankscore < 17.765
   rule_3647   0.72136270   Adipose_Subcutaneous >= 0.001 & Pituitary >= 0.003 & CADD_raw_rankscore < 17.96
  rule_14290   0.69403493   Muscle >= 0.0011 & Adrenal_Gland >= 0.0009 & CADD_raw_rankscore < 17.945
   rule_9185   0.56068510   Pituitary >= 0.001 & CADD_raw_rankscore < 13.725 & CADD_raw_rankscore < 21.65
  rule_12740   0.47650050   CADD_raw_rankscore < 13.43 & Pituitary >= 0.0007 & CADD_raw_rankscore < 20.85
   rule_5233   0.42753515   Adrenal_Gland >= 0.0005 & Brain_Cerebellum >= 0.0001 & CADD_raw_rankscore < 18.575
   rule_8405   0.11213065   Thyroid >= 0.0007 & Artery_Tibial >= 0.0004 & CADD_raw_rankscore < 17.665
   rule_4554   0.05239178   Testis < 0.99 & Adipose_Subcutaneous >= 0.0005 & CADD_raw_rankscore < 18.525
                     
=========================================

Model call: 
rule_fit(data = gtex_data, y = "Pathogenic", two_level_bootstrap = TRUE, 
    ntree = ntree, max.depth = maxdepth, rule.filter = rules, 
    bootstrap.rounds = lasso, rule.extract.cores = 16, kappa.cores = 10, 
    lasso.cores = 16, ranger.and.rulefit = FALSE, report = FALSE)
```

## Predict
Once the model is trained on the HPO term related genes, we can use the model to predict on patient data:

```
predict.varppRule(patient_data, 
                  model_results, 
                  predict=c("probability", "class"))
```

## License

The contents of this repository are distributed under the MIT license. See below for details:
```
The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
