#' Read GTEx or HCL data
#'
#' This function is needed to read in the GTEx or HCL data for modelling
#'
#' @param expression_file is the filepath to the expression data (GTEx or HCL)
#' @param benign_variant_file is the file path to the benign variant file
#' @param pathogenic_variant_file is the filepath to the pathogenic variant file'
#'
#' @return a list of expression data, pathogenic variants and benign variants
#'
#' @import data.table
#' @import tidyverse
#' @export
.gtex_preparation <- function(expression_file, benign_variant_file, pathogenic_variant_file){


  message("Loading GTEx data...")


  expression                  <- read.csv(expression_file, stringsAsFactors=FALSE)
  benign_variants             <- read.table(benign_variant_file, header=TRUE, stringsAsFactors=FALSE)
  benign_variants$Pathogenic  <- 0

  disease_variants            <- data.table::fread(input=pathogenic_variant_file, header=TRUE, na.strings=".")
  colnames(disease_variants)[colnames(disease_variants)== "Gene"] <- "genename"
  disease_variants <- as.data.frame(disease_variants, stringsAsFactors=FALSE)

  # Sequentially number variants within the same gene as this will be required for sampling
  disease_variants$variant <- as.numeric(ave(disease_variants$genename, disease_variants$genename, FUN=seq_along))
  colnames(disease_variants)[colnames(disease_variants)== "genename"] <- "Gene"

  disease_variants$GeneVariant <- paste(disease_variants$Gene, disease_variants$variant, sep="_")

  disease_variants$Pathogenic <- 1

  # Drop any genes in the benign variant set that are in the disease variant set. This will only be relevant if using later versions of Phenolyzer and dbNSFP as there will be no overlap when using stated versions.
  benign_variants <- benign_variants[!(benign_variants$GTExV8_Gene %in% disease_variants$Gene), ]

  exclude_cols <- names(expression)[43:length(names(expression))]

  #=====================================================================================
  # MAKE SURE ALL DATA FRAMES HAVE THE SAME COLNAMES FOR GENES AND OTHER SHARED VARIABLE
  #=====================================================================================

  disease_variants <- disease_variants[,c('GTExV8_Gene', 'GeneVariant', 'Pathogenic', 'MetaSVM_rankscore', 'CADD_raw_rankscore')]
  disease_variants  <- disease_variants %>%
    rename(Gene = GTExV8_Gene) # necessary as the Gene name column is not called Gene here

  benign_variants  <- benign_variants[,c('GTExV8_Gene', 'GeneVariant', 'Pathogenic', 'MetaSVM_rankscore', 'CADD_raw_rankscore')]
  benign_variants  <- benign_variants %>%
    rename(Gene = GTExV8_Gene) # necessary as the Gene name column is not called Gene here

  expression <-
    expression %>%
    rename(Gene = gene_name) %>%
    select(-c("chromosome_name", "start_position",  "end_position", "gene_id", "gene_biotype", "Pathogenic" ))

  # Merge all data
  sub_benign  <- bind_cols(benign_variants, expression[match(gsub("_.*$", "", benign_variants$GeneVariant), expression$Gene),!(colnames(expression) %in% 'Gene')])
  sub_patho   <- bind_cols(disease_variants, expression[match(gsub("_.*$", "", disease_variants$GeneVariant), expression$Gene), !(colnames(expression) %in% 'Gene')])

  dat         <- bind_rows(sub_patho, sub_benign)

  list(dat=dat,disease_variants=disease_variants, benign_variants=benign_variants)

}




#' Read GTEx or HCL data
#'
#' This function is needed to read in the GTEx data without two level bootstrap. Here we jsut want to see if we can
#' identify the disease gene
#'
#' @param expression_file is the filepath to the expression data (GTEx or HCL)
#' @param benign_variant_file is the file path to the benign variant file
#' @param pathogenic_variant_file is the filepath to the pathogenic variant file'
#'
#' @return a list of expression data, pathogenic variants and benign variants
#'
#' @import data.table
#' @import tidyverse
#' @export
.gtex_gene_level <- function(expression_file,
                             benign_variant_file,
                             pathogenic_variant_file){

  message("Loading GTEx data on gene level only...")


  expression                  <- read.csv(expression_file, stringsAsFactors=FALSE)
  benign_variants             <- read.table(benign_variant_file, header=TRUE, stringsAsFactors=FALSE)
  benign_variants$Pathogenic  <- 0

  disease_variants            <- data.table::fread(input=pathogenic_variant_file, header=TRUE, na.strings=".")
  colnames(disease_variants)[colnames(disease_variants)== "Gene"] <- "genename"
  disease_variants <- as.data.frame(disease_variants, stringsAsFactors=FALSE)

  # Sequentially number variants within the same gene as this will be required for sampling
  disease_variants$variant <- as.numeric(ave(disease_variants$genename, disease_variants$genename, FUN=seq_along))
  colnames(disease_variants)[colnames(disease_variants)== "genename"] <- "Gene"

  disease_variants$GeneVariant <- paste(disease_variants$Gene, disease_variants$variant, sep="_")

  disease_variants$Pathogenic <- 1

  # Drop any genes in the benign variant set that are in the disease variant set. This will only be relevant if using later versions of Phenolyzer and dbNSFP as there will be no overlap when using stated versions.
  benign_variants <- benign_variants[!(benign_variants$GTExV8_Gene %in% disease_variants$Gene), ]

  exclude_cols <- names(expression)[43:length(names(expression))]

  #=====================================================================================
  # MAKE SURE ALL DATA FRAMES HAVE THE SAME COLNAMES FOR GENES AND OTHER SHARED VARIABLE
  #=====================================================================================

  disease_variants <- disease_variants[,c('GTExV8_Gene', 'GeneVariant', 'Pathogenic', 'MetaSVM_rankscore', 'CADD_raw_rankscore')]
  disease_variants  <- disease_variants %>%
    rename(Gene = GTExV8_Gene) # necessary as the Gene name column is not called Gene here

  benign_variants  <- benign_variants[,c('GTExV8_Gene', 'GeneVariant', 'Pathogenic', 'MetaSVM_rankscore', 'CADD_raw_rankscore')]
  benign_variants  <- benign_variants %>%
    rename(Gene = GTExV8_Gene) # necessary as the Gene name column is not called Gene here

  expression <-
    expression %>%
    rename(Gene = gene_name) %>%
    select(-c("chromosome_name", "start_position",  "end_position", "gene_id", "gene_biotype", "Pathogenic" ))

  # Merge all data
  sub_benign  <- bind_cols(benign_variants, expression[match(gsub("_.*$", "", benign_variants$GeneVariant), expression$Gene),!(colnames(expression) %in% 'Gene')])
  sub_patho   <- bind_cols(disease_variants, expression[match(gsub("_.*$", "", disease_variants$GeneVariant), expression$Gene), !(colnames(expression) %in% 'Gene')])

  dat         <- bind_rows(sub_patho, sub_benign)
  # Retain only
  dat <- dat %>%
    select(-c(GeneVariant, CADD_raw_rankscore, MetaSVM_rankscore)) %>%
    distinct()

  # We shuffle the data set to make sure not all pathogenic variants are stacked on top of the benign variants.
  # This should be handy for the modelling
  set.seed(5) # there is no need for this step to be not reproducible; in fact it will be easier if it is fixed.
  rows <- sample(nrow(dat))
  dat <- na.omit(dat[rows,])
  dat


  #list(dat=dat,disease_variants=disease_variants, benign_variants=benign_variants)
}
