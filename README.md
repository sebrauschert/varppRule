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
