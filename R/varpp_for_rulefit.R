#'varpp: extract rules from ranger trees
#'
#' This function is meant to only be used internally for the rule_fit function
#'
#' @param dat this is a data list returned from the function load_gtex_or_hcl. It is either GTEx tissue specific gene expression or HCL cell specific expression
#' @param ntree is the number of trees that should be built for ranger. It defaults to 1000
#' @param max.depth is the maximum tree depth for the ranger trees. IT defaults to 3.
#' @param cores number of cores for parallel, defaults to 4
#' @import progress
#' @import doMC
#' @import doParallel
#' @import parallel
#' @import tidyverse
#' @import dplyr
#' @import ranger
#' @import magrittr
#' @importFrom magrittr "%>%"
#' @import foreach
#' @importFrom iterators icount
#' @export
varpp_for_rulefit <- function(dat,
                              ntree,
                              max.depth,
                              cores){



  # Spedify the benign and pathogenic gene names
  cls_pathogenic_genes <- unique(dat$dat$Gene[dat$dat$Pathogenic %in% 1])
  cls_benign_genes     <- unique(dat$dat$Gene[dat$dat$Pathogenic %in% 0])

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
    cls_benign <- sample(cls_benign_genes, replace=TRUE)

    benign_variants_sub <- dat$benign_variants[dat$benign_variants$Gene %in% unique(cls_benign),
                                               c("Gene","GeneVariant","Pathogenic","CADD_raw_rankscore")]

    # The new sampling step based on the function in utilities
    sub_benign <- .sample_benign_variants(benign_variants_sub, cls_benign)

    # This is still necessary to make sure that only one variant is ever selected per gene
    sub_benign <- sub_benign[match(gsub("_.*$", "", sub_benign), gsub("_.*$", "", sub_benign))]
    WEIGHTS    <- rbind(data.frame(plyr::count(sub_benign)),data.frame(plyr::count(sub_patho)))

    # The merge function does not preserve the order of the original 'dat' data.frame.
    # That means, the sampling 'inbag' order we use to subset later is not in the right order
    # and the wrong Gene Variants are selected. Solution: create index id and sort the merged data after

    dat$dat$id     <- seq_len(nrow(dat$dat))
    weighted_data  <- merge(dat$dat, WEIGHTS, by.x="GeneVariant", by.y="x", all=T)
    dat$dat$id     <- NULL
    weighted_data  <- weighted_data[order(weighted_data$id), ] %>%
      select(-id)
    inbag          <-
      as.numeric(unlist(weighted_data %>%
                          mutate(freq = ifelse(freq %in% NA, 0, freq)) %>%
                          rename(weight = freq) %>%
                          select(weight)))

    # Remove redundancies
    rm(sub_patho, sub_benign, benign_variants_sub, WEIGHTS)
    gc()

    # This step is crucial: the index so far only collects the position in the original data
    # set, but some Variants are sampled multiple times, hence we need to duplciate the selection of
    # those variants n times, as extracted by the count in hte previous loop

    indices <- rep(which(inbag >= 1 ), inbag[which(inbag >= 1 )])

    # Specify the in bag samples: based on the indices
    dat_in        <- dat$dat[indices,]
    dat_in$weight <- 1
    dat_in$Gene   <- NULL

    # Remove those variants that are in the same genes as the ones selected in this round of the loop
    dat_out        <- dat$dat[!(gsub("_.*$", "",dat$dat$GeneVariant) %in% gsub("_.*$", "",dat_in$GeneVariant)), ]
    dat_out$weight <- 0
    dat_out$Gene   <- NULL

    dat_boot <- list(training=data.frame(bind_rows(dat_in, dat_out), check.names=FALSE))

    rm(dat_in, dat_out)
    gc()

    dat_boot$training$Pathogenic <- factor(dat_boot$training$Pathogenic)

    # Random forest does not handle missing data
    dat_boot$training_CADD <- na.omit(dat_boot$training[ , colnames(dat_boot$training) != "MetaSVM_rankscore"])
    dat_boot$training_CADD <- dat_boot$training_CADD %>% arrange(desc = GeneVariant)
    dat_boot$training      <- NULL


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

    # This custom function extracts the rules from the ranger trees.
    .extract_ranger_rules(rf_results)


  }
  stopCluster(cl)

  # Collect the results for ranger at every step and put it together

  RULES        <- do.call(c, inbag) #Reduce(c,inbag)
  names(RULES) <- paste0("rule_",1:length(RULES))

  RULES

}


#===========#
#           #
#  /(_M_)\  #
# |       | #
#  \/~V~\/  #
#           #
#===========#
