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
        cmd_plink2 <- paste(configs[["plink2.path"]], "--silent", "--threads", configs[["nCores"]], 
            "--pfile", pfile, ifelse(configs[["vzs"]], "vzs", ""), "--keep", keep_f, 
            "--out", configs[["gcount.full.prefix"]], "--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs")
        if (!is.null(configs[["mem"]])) 
            cmd_plink2 <- paste(cmd_plink2, "--memory", configs[["mem"]])
        
        snpnetLogger(sprintf("Running plink2: %s", cmd_plink2))
        system(cmd_plink2, intern = F, wait = T)
        
        # read the gcount file
        gcount_df <- data.table::fread(paste0(configs[["gcount.full.prefix"]], ".gcount")) %>% 
            dplyr::rename(original_ID = ID) %>% dplyr::mutate(
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
