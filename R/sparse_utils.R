computeSparseStats <- function(pfile, ids, configs) {
    keep_f <- paste0(configs[["gcount.full.prefix"]], ".keep")
    gcount_tsv_f <- paste0(configs[["gcount.full.prefix"]], ".gcount.tsv")
    
    dir.create(dirname(configs[["gcount.full.prefix"]]), showWarnings = FALSE, recursive = TRUE)
    if (file.exists(gcount_tsv_f)) {
        gcount_df <- data.table::fread(gcount_tsv_f)
    } else {
        # To run plink2 --geno-counts, we write the list of IDs to a file
        data.frame(ID = ids) %>% tidyr::separate(ID, into = c("FID", "IID"), sep = "_") %>% 
            data.table::fwrite(keep_f, sep = "\t", col.names = F)
        
        # Run plink2 --geno-counts
        cmd_plink2 <- paste(
            configs[["plink2.path"]], "--silent", 
            "--threads", configs[["nCores"]], 
            "--pfile", pfile, ifelse(configs[["vzs"]], "vzs", ""), 
            "--keep", keep_f, 
            "--out", configs[["gcount.full.prefix"]], 
            "--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs"
            )

        if (!is.null(configs[["mem"]])) 
            cmd_plink2 <- paste(cmd_plink2, "--memory", configs[["mem"]])
        
        snpnetLogger(sprintf("Running plink2: %s", cmd_plink2))
        system(cmd_plink2, intern = F, wait = T)
        
        # read the gcount file stats_means/2 is the MAF
        gcount_df <- data.table::fread(paste0(configs[["gcount.full.prefix"]], ".gcount")) %>% 
            dplyr::rename(original_ID = ID) %>% 
            dplyr::mutate(
                ID = paste0(original_ID, "_", ALT), 
                stats_pNAs = MISSING_CT/(MISSING_CT + OBS_CT), 
                NON_REF_CT = HAP_ALT_CTS + HET_REF_ALT_CTS + TWO_ALT_GENO_CTS, 
                miss_over_non_ref = MISSING_CT/(HAP_ALT_CTS + HET_REF_ALT_CTS + TWO_ALT_GENO_CTS), 
                stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS)/OBS_CT
                )

        gcount_df$index <- seq(nrow(gcount_df))
    }
    
    if (configs[["save"]]) {
        gcount_df %>% data.table::fwrite(gcount_tsv_f, sep = "\t")
    }
    gcount_df
}


getCoxResponseObj <- function(surv, offset = NULL) {
    ord <- order(surv[, 1])
    sorted_y <- surv[ord, 1]
    rankmin <- rank(sorted_y, ties.method = "min") - 1L
    rankmax <- rank(sorted_y, ties.method = "max") - 1L
    ord <- ord - 1L
    normalized_status <- surv[, 2]/sum(surv[, 2])
    if ((!is.null(offset)) && (length(offset) != nrow(surv))) {
        stop("offset does not have the same size as the surv object")
    }
    cox_obj <- pgenlibr::NewCoxResponseObj(normalized_status, ord, rankmin, rankmax, offset)
    cox_obj
}

# Very important: the new genetic data MUST share the same pvar as the training
# genetic data. This function does not perform check of this precondition.
#' @export
sparse_predict <- function(fit_obj, genotype.pfile, phenotype.file, split.col = "split", 
    split = "test", status = NULL) {
    family <- fit_obj$covs_fit$family[[1]]
    if (family != "cox") {
        m <- names(fit_obj$covs_fit$model)
    } else {
        m <- fit_obj$covs_fit$model
        if (is.null(status)) {
            stop("Cox model must specify a status column")
        }
    }
    phenotype <- m[1]
    covs <- m[2:length(m)]

    psamid <- readIDsFromPsam(paste0(genotype.pfile, ".psam"))
    phe <- readPheMaster(phenotype.file, psamid, family, covs, phenotype, status, 
        split.col, list(zstdcat.path='zstdcat', zcat.path='zcat'))
    phe <- phe[phe[[split.col]] == split, ]
    if (family != "cox") {
        covs_pred <- as.matrix(stats::predict(fit_obj$covs_fit, (phe %>% dplyr::select(all_of(covs)))))
        response <- phe[[phenotype]]
    } else {
        covs_pred <- as.matrix(phe %>% dplyr::select(all_of(covs))) %*% fit_obj$covs_fit$beta
        response <- survival::Surv(phe[[phenotype]], phe[[status]])
    }
    
    pgen <- pgenlibr::NewPgen(paste0(genotype.pfile, ".pgen"), sample_subset = match(phe$ID, 
        psamid))
    Xtest <- pgenlibr::NewDense(pgen, fit_obj$snps_used$index, fit_obj$snps_used$stats_means)
    
    nlambda <- ncol(fit_obj$beta)
    prediction <- matrix(nrow = nrow(phe), ncol = nlambda)
    for (i in 1:nlambda) {
        prediction[, i] <- pgenlibr::DenseMultv(Xtest, fit_obj$beta[, i])
    }
    prediction <- sweep(prediction, 1, covs_pred, "+")
    if (family == "gaussian") {
        metric.type <- "r2"
    } else if (family == "binomial") {
        metric.type <- "auc"
    } else {
        metric.type <- "C"
    }
    metric.test <- computeMetric(prediction, response, metric.type)
    
    return(list(prediction = prediction, metric.test = metric.test))
}
