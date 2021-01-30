#' @export
sparse_snpnet <- function(genotype.pfile, phenotype.file, phenotype, group_map, family = NULL, 
    covariates = NULL, nlambda = 100, lambda.min.ratio = 0.01, lambda = NULL, split.col = NULL, 
    p.factor = NULL, status.col = NULL, mem = NULL, configs = NULL) {
    time.start <- Sys.time()
    snpnetLogger("Start snpnet", log.time = time.start)
    
    snpnetLogger("Preprocessing start..")
    
    ### --- Read genotype IDs --- ###
    ids <- list()
    phe <- list()
    ids[["psam"]] <- readIDsFromPsam(paste0(genotype.pfile, ".psam"))
    
    ### --- combine the specified configs with the default values --- ###
    if (!is.null(lambda)) 
        nlambda <- length(lambda)
    configs <- setupConfigs(configs, genotype.pfile, phenotype.file, phenotype, covariates, 
        1, nlambda, split.col, p.factor, status.col, mem)
    ### --- Read phenotype file --- ###
    phe[["master"]] <- readPheMaster(phenotype.file, ids[["psam"]], family, covariates, 
        phenotype, status.col, split.col, configs)
    
    ### --- infer family and update the configs --- ###
    if (is.null(family)) 
        family <- inferFamily(phe[["master"]], phenotype, status.col)
    configs <- updateConfigsWithFamily(configs, family)
    
    ### --- Process phenotypes --- ###
    if (family == "binomial") {
        # The input binary phenotype is coded as 2/1 (case/control) For glmnet, we map
        # this to 1/0 (case/control) The following expression will replace -9 (missing)
        # with -10, but the set of individuals with no-missing values are already
        # computed.
        if (min(phe[["master"]][[phenotype]], na.rm = T) >= 1 && max(phe[["master"]][[phenotype]], 
            na.rm = T) <= 2) {
            phe[["master"]][[phenotype]] <- phe[["master"]][[phenotype]] - 1
        }
    }
    
    ### --- Define the set of individual IDs for training (and validation) set(s) ---
    validation <- (!is.null(split.col))
    if (validation) {
        splits <- c("train", "val")
        for (s in splits) {
            ids[[s]] <- phe[["master"]]$ID[phe[["master"]][[split.col]] == s]
        }
    } else {
        splits <- c("train")
        ids[["train"]] <- phe[["master"]]$ID
    }
    
    ### --- Prepare the feature matrix --- ###
    features <- list()
    for (s in splits) {
        phe[[s]] <- phe[["master"]][match(ids[[s]], phe[["master"]]$ID), ]
        rownames(phe[[s]]) <- phe[[s]]$ID
        if (length(covariates) > 0) {
            features[[s]] <- phe[[s]][, covariates, with = F]
        } else {
            features[[s]] <- NULL
        }
        if (configs[["verbose"]]) 
            snpnetLogger(sprintf("The number of individuals in %s set: %d", s, dim(phe[[s]])[1]))
    }
    
    ### --- Prepare the response --- ###
    response <- list()
    status <- list()
    pred <- list()
    for (s in splits) {
        response[[s]] <- phe[[s]][[phenotype]]
        if (family == "cox") {
            status[[s]] <- phe[[s]][[status.col]]
            response[[s]] <- survival::Surv(response[[s]], status[[s]])
        }
    }
    
    ### --- Read genotypes --- ###
    vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd = paste0(configs[["zstdcat.path"]], 
        " ", paste0(genotype.pfile, ".pvar.zst"))), CHROM = "#CHROM"), VAR_ID = paste(ID, 
        ALT, sep = "_"))$VAR_ID
    pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, ".pvar.zst"))
    
    pgen <- list()
    samples_subset <- list()
    for (s in splits) {
        samples_subset[[s]] <- match(ids[[s]], ids[["psam"]])
        pgen[[s]] <- pgenlibr::NewPgen(paste0(genotype.pfile, ".pgen"), pvar = pvar, 
            sample_subset = samples_subset[[s]])
    }
    pgenlibr::ClosePvar(pvar)
    
    snps_to_use <- computeSparseStats(genotype.pfile, phe[["train"]]$ID, configs = configs)
    
    snps_to_use <- snps_to_use %>% dplyr::filter((NON_REF_CT >= 3) & (stats_pNAs < 
        configs[["missing.rate"]]) & (miss_over_non_ref < 10))
    snps_to_use <- snps_to_use %>% dplyr::select(c("#CHROM", "POS", "original_ID", 
        "index", "NON_REF_CT", "stats_means")) %>% tidyr::separate(original_ID, into = c("CHROM", 
        "POS"), sep = ":", remove = FALSE, convert = TRUE, extra = "drop")
    
    ### Temporary solution, the mapping file must have these columns
    gene_map <- data.table::fread(group_map, select = c("#CHROM", "POS", "SYMBOL", 
        "Csq"))
    names(gene_map)[1] <- "CHROM"
    names(gene_map)[3] <- "gene_symbol"
    
    gene_map$CHROM[gene_map$CHROM == "X"] <- "23"
    gene_map$CHROM[gene_map$CHROM == "Y"] <- "24"
    gene_map$POS <- as.integer(gene_map$POS)
    gene_map$CHROM <- as.integer(gene_map$CHROM)
    unique_chrom <- unique(gene_map$CHROM)
    if (max(unique_chrom) != length(unique_chrom)) {
        stop("The chromosome must be stored in the order 1,2,3,..., X, Y")
    }
    
    cumu_chrom <- integer(length(unique_chrom) + 1)
    for (i in (1:length(unique_chrom))) {
        ind <- which(gene_map$CHROM == i)
        cumu_chrom[i + 1] <- cumu_chrom[i] + length(ind)
        
        requirement <- all(diff(gene_map$POS[ind]) >= 0) && (ind[1] == cumu_chrom[i] + 
            1)
        if (!requirement) {
            stop("The gene mapping must have all SNPs on the same chromosome stored contiguously, and the position of the SNPs within the same chromosome must be non-decreasing")
        }
    }
    snps_gene_ind <- pgenlibr::match_sorted_snp(snps_to_use$CHROM, snps_to_use$POS, 
        gene_map$POS, cumu_chrom)
    
    # Just for testing,
    print(all(snps_to_use$CHROM == gene_map$CHROM[snps_gene_ind]))
    print(all(snps_to_use$POS == gene_map$POS[snps_gene_ind]))
    snps_to_use$gene <- gene_map$gene_symbol[snps_gene_ind]
    snps_to_use$Csq <- gene_map$Csq[snps_gene_ind]
    
    # Only use autosome
    snps_to_use <- snps_to_use %>% dplyr::filter(CHROM < 23) %>% dplyr::filter(gene != 
        "")
    # snps_to_use = snps_to_use %>% filter(Csq %in% c('ptv', 'pav'))
    
    genes_with_comma <- grepl(",", snps_to_use$gene, fixed = TRUE)
    if (any(genes_with_comma)) {
        warning("Some SNPs are mapped to multiple genes. The first gene symbols will be used to group these SNPs")
        snps_to_use$gene[genes_with_comma] <- sapply(strsplit(snps_to_use$gene[genes_with_comma], 
            ","), "[", 1)
    }
    
    unique_genes <- unique(snps_to_use$gene)
    snps_to_use$gene_order <- match(snps_to_use$gene, unique_genes)
    snps_to_use <- snps_to_use %>% dplyr::arrange(gene_order)
    
    gene_cumu <- snps_to_use %>% dplyr::count(gene_order)
    gene_cumu <- c(0, cumsum(gene_cumu$n))
    
    snpnetLoggerTimeDiff("Preprocessing end.", time.start, indent = 1)
    
    
    ### --- Fit a model using only the covariates
    if (family == "cox") {
        stop("I will implement this later")
    } else {
        glmmod <- stats::glm(stats::as.formula(paste(phenotype, " ~ ", paste(c(1, 
            covariates), collapse = " + "))), data = phe[["train"]], family = family)
    }

    proxObj = pgenlibr::NewProxObj(nrow(snps_to_use), gene_cumu)
    
    
    gaussian_response_sd <- NULL
    responseObj <- NULL
    offset <- list()
    # Run linear prediction
    if (family == "gaussian") {
        gaussian_response_sd = sd(glmmod$residuals)
        responseObj = pgenlibr::NewResponseObj(glmmod$residuals/gaussian_response_sd, "gaussian")
    } else {
        for (s in split) {
            offset[[s]] <- stats::predict(glmmod, phe[[s]] %>% dplyr::select(all_of(covs)))
        }
       # Figure this out later
    }
    

    time.load.matrix <- Sys.time()
    snpnetLogger("Start loading training genotype matrix", log.time = time.start)
    
    Xtrain <- pgenlibr::NewSparse(pgen[["train"]], snps_to_use$index)
    snpnetLoggerTimeDiff("End loading genotype matrix", time.load.matrix, indent = 1)

    if (is.null(lambda)) {
        lambda.max <- pgenlibr::ComputeLambdaMax(Xtrain, responseObj, gene_cumu)
        full.lams <- exp(seq(from = log(lambda.max), to = log(lambda.max * lambda.min.ratio), 
            length.out = nlambda))
    } else {
        full.lams <- lambda
    }

    outer.idx = 1
    num.lambda.finished = 0

    # while(num.lambda.finished < length(full.lams)){

    # }
    
    start <- Sys.time()
    #result <- pgenlibr::SparseTest123(Xtrain, response[["train"]], gene_cumu, full.lams[1:30])
    print(proxObj)
    result <- pgenlibr::FitGroupLasso(Xtrain, proxObj, responseObj,full.lams[1:16])
    print(Sys.time() - start)
    print(proxObj)
    result2 <- pgenlibr::FitGroupLasso(Xtrain, proxObj, responseObj,full.lams[17:30])
    print(Sys.time() - start)
    return(0)
    
}
