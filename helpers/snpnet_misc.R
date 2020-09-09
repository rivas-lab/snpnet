suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))

config_params_data_type <- function(){
    list(
        double = c(
            'missing.rate',
            'MAF.thresh',
            'glmnet.thresh',
            'lambda.min.ratio',
            'alpha',
            'p.factor'
        ),
        integer = c(
            'nCores',
            'prevIter',
            'niter',
            'mem',
            'nlambda',
            'nlams.init',
            'nlams.delta',
            'num.snps.batch',
            'increase.size',
            'stopping.lag'
        ),
        logical = c(
            'use.glmnetPlus',
            'standardize.variant',
            'validation',
            'early.stopping',
            'vzs',
            'save',
            'save.computeProduct',
            'verbose',
            'KKT.verbose',
            'KKT.check.aggressive.experimental',
            'rank'
        )
    )
}

parse_covariates <- function(covariates_str){
    if(is.null(covariates_str) || covariates_str == 'None'){
        covariates=c()
    }else{
        covariates=strsplit(covariates_str, ',')[[1]]
    }
    covariates
}

read_lambda_sequence_from_RData_file <- function(rdata_f){
    load(rdata_f)
    fit$full.lams[1:which.max(fit$metric.val)]
}

read_config_from_file <- function(config_file){
    null_strs <- c('null', 'Null', 'NULL', 'none', 'None', 'NONE')

    config_df <- config_file %>% fread(header=T, sep='\t') %>%
    setnames(c('key', 'val'))

    config <- as.list(setNames(config_df$val, config_df$key))
    config_dt <- config_params_data_type()
    for(k in intersect(config_dt[['double']], names(config))){
        config[[k]] <- as.double(config[[k]])
    }
    for(k in intersect(config_dt[['integer']], names(config))){
        config[[k]] <- as.integer(config[[k]])
    }
    for(k in intersect(config_dt[['logical']], names(config))){
        config[[k]] <- as.logical(config[[k]])
    }
    for(key in c('status.col', 'covariates', 'split.col', 'keep')){
        if( (! key %in% names(config)) | (config[[ key ]] %in% null_strs) | is.na(config[[ key ]]) ){
            config[[ key ]] <- NULL
        }
    }
    if(!'validation' %in% names(config)) config[['validation']] <- FALSE
    if(!config[['validation']]) config[['split.col']] <- NULL
    config[['covariates']] = parse_covariates(config[['covariates']])

    if(
        # read penalty factor from a RDS file
        ('p.factor.file' %in% names(config)) &&
        (! is.null(config[['p.factor.file']])) &&
        (! config[['p.factor.file']] %in% null_strs)
    ){
        config[['p.factor']] <- readRDS(
            config[['p.factor.file']]
        )
    }else if (! 'p.factor' %in% names(config)) {
        config[['p.factor']] <- NULL
    }

    if(
        # for refit, read the lambda sequence from the specified RData file
        ('refit_from_RData' %in% names(config)) &&
        (! is.null(config[['refit_from_RData']])) &&
        (! config[['refit_from_RData']] %in% null_strs)
    ){
        config[['lambda']] <- read_lambda_sequence_from_RData_file(config[['refit_from_RData']])
        config[['split.col']] <- NULL
    }else if (! 'lambda' %in% names(config)) {
        config[['lambda']] <- NULL
    }

    config
}

snpnet_fit_to_df <- function(beta, lambda_idx, covariates = NULL, verbose=FALSE){
    if(verbose) snpnet::snpnetLogger(sprintf('Extracting the BETAs (lambda idx: %d)..', lambda_idx), funcname='snpnetWrapper')
    # extract BETAs from snpnet(glmnet) fit as a data frame
    df <- beta[lambda_idx] %>% data.frame()
    colnames(df) <- 'BETA'
    df <- df %>% rownames_to_column("ID")

    df$BETA <- format(df$BETA, scientific = T)

    non.zero.BETAs <- union(
        covariates,
        df %>% filter(as.numeric(BETA) != 0) %>% select(ID) %>% pull()
    )
    if(verbose) snpnet::snpnetLogger(sprintf(
        'The BETAs are extracted for %d variables (%d covariates and %d genetic variants).',
        length(non.zero.BETAs), length(intersect(non.zero.BETAs, covariates)), length(setdiff(non.zero.BETAs, covariates))
    ), indent=1, funcname='snpnetWrapper')
    df %>% filter(ID %in% non.zero.BETAs)
}

save_BETA <- function(df, out.file.head, pvar, vzs = TRUE, covariates = NULL, verbose=FALSE){
    file.geno   <- paste0(out.file.head, ".tsv")
    file.covars <- paste0(out.file.head, ".covars.tsv")
    file.timestamp <- format(Sys.time(), format="%Y%m%d-%H%M%S")

    if(! is.null(covariates)){
        if(verbose) snpnet::snpnetLogger(sprintf('Saving results to %s', file.covars), funcname='snpnetWrapper')
        if(file.exists(file.covars)){
            file.rename(file.covars, paste0(out.file.head, ".", file.timestamp, ".covars.tsv"))
        }
        df %>% filter(ID %in% covariates) %>%
        fwrite(file.covars, sep='\t')
    }

    catcmd <- ifelse(vzs, 'zstdcat', 'cat')

    pvar_df <- fread(cmd=paste(catcmd, pvar, '| sed -e "s/^#//g"', sep=' ')) %>%
    mutate(varID = paste(ID, ALT, sep='_'))

    if(verbose) snpnet::snpnetLogger(sprintf('Saving results to %s', file.geno), funcname='snpnetWrapper')

    if(file.exists(file.geno)){
        file.rename(file.geno, paste0(out.file.head, ".", file.timestamp, ".tsv"))
    }

    df %>% filter(! ID %in% covariates) %>%
    rename('varID' = 'ID') %>%
    left_join(pvar_df, by='varID') %>%
    arrange(CHROM, POS) %>%
    select(c(colnames(pvar_df %>% select(-varID)), 'BETA')) %>%
    fwrite(file.geno, sep='\t')
}

compute_covar_score <- function(phe_df, covar_df){
    phe_df %>%
    mutate(ID = paste0(FID, '_', IID)) %>%
    select(covar_df %>% select(ID) %>% pull() %>% c('ID')) %>%
    column_to_rownames('ID') %>% as.matrix() %*%
    as.matrix(covar_df %>% column_to_rownames('ID')) %>%
    as.data.frame() %>%
    rownames_to_column('ID') %>%
    rename(covar_score = BETA) %>%
    separate(ID, into=c('FID', 'IID'), sep='_') %>%
    select(FID, IID, covar_score)
}

compute_score <- function(phe_df, covar_beta_file, geno_score_file){
    covar_df <- fread(covar_beta_file, sep='\t', colClasses=c(ID="character", BETA="numeric"))

    covar_score_df <- compute_covar_score(phe_df, covar_df)

    fread(
        cmd=paste0('cat ', geno_score_file, ' | sed -e "s/^#//g"'),
        colClasses=c(FID="character", IID="character")
    ) %>%
    rename('geno_score' = 'SCORE1_SUM') %>%
    select(FID, IID, geno_score) %>%
    inner_join(covar_score_df, by=c('FID', 'IID')) %>%
    mutate(score = geno_score + covar_score)
}

compute_auc <- function(df, split_str, score_col){
    df_sub <- df %>%
    filter(split == split_str & phe != -9) %>%
    mutate(phe = as.factor(phe)) %>%
    rename(auc_response = score_col)

    auc(
        df_sub %>% select(phe) %>% pull(),
        df_sub %>% select(auc_response) %>% pull()
    )
}

### old functions
#
# filter_by_percentile_and_count_phe <- function(df, p_l, p_u){
#     df %>% filter(p_l < Percentile, Percentile <= p_u) %>%
#     count(phe)
# }

# compute_OR <- function(df, l_bin, u_bin, cnt_middle){
#     cnt_tbl <- df %>%
#     filter_by_percentile_and_count_phe(l_bin, u_bin) %>%
#     inner_join(cnt_middle, by='phe') %>% gather(bin, cnt, -phe) %>%
#     arrange(-phe, bin)

#     cnt_res <- cnt_tbl %>% mutate(cnt = as.numeric(cnt)) %>% select(cnt) %>% pull()
#     names(cnt_res) <- c('n_TP', 'n_FN', 'n_FP', 'n_TN')

#     OR <- (cnt_res[['n_TP']] * cnt_res[['n_TN']]) / (cnt_res[['n_FP']] * cnt_res[['n_FN']])
#     LOR <- log(OR)
#     se_LOR <- cnt_tbl %>% select(cnt) %>% pull() %>%
#     lapply(function(x){1/x}) %>% reduce(function(x, y){x+y}) %>% sqrt()
#     l_OR = exp(LOR - 1.96 * se_LOR)
#     u_OR = exp(LOR + 1.96 * se_LOR)

#     data.frame(
#         l_bin = l_bin,
#         u_bin = u_bin,
#         n_TP = cnt_res[['n_TP']],
#         n_FN = cnt_res[['n_FN']],
#         n_FP = cnt_res[['n_FP']],
#         n_TN = cnt_res[['n_TN']],
#         OR   = OR,
#         SE_LOR = se_LOR,
#         l_OR = l_OR,
#         u_OR = u_OR,
#         OR_str = sprintf('%.3f (%.3f-%.3f)', OR, l_OR, u_OR)
#     ) %>%
#     mutate(OR_str = as.character(OR_str))
# }

# compute_OR_tbl <- function(df, split_str, score_col){
#     all_ranked_df <- df %>% filter(split == split_str) %>%
#     rename(auc_response = score_col) %>%
#     mutate(Percentile = rank(-auc_response) / n()) %>%
#     filter(phe != -9)

#     cnt_middle <- all_ranked_df %>%
#     filter_by_percentile_and_count_phe(0.4, 0.6) %>%
#     rename('n_40_60' = 'n')

#     bind_rows(lapply(1:20, function(x){
#     compute_OR(all_ranked_df, (x-1)/20, x/20, cnt_middle)
#     }),
#     compute_OR(all_ranked_df, 0, .001, cnt_middle),
#     compute_OR(all_ranked_df, 0, .01, cnt_middle),
#     compute_OR(all_ranked_df, .99, 1, cnt_middle),
#     compute_OR(all_ranked_df, .999, 1, cnt_middle)
#     )
# }

# We incoporated the latest changes in the helper functions
# https://github.com/rivas-lab/covid19/blob/a5592ae4f59da4479847c05338e7621dab6788cd/snpnet/functions.R#L89

### functions for evaluation (PRS bin vs. OR and mean)

compute_mean <- function(df, percentile_col, phe_col, l_bin, u_bin){
    # Compute the mean and sd of the trait value (phe_col), based on the
    # binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
    stratified_df <- df %>%
    rename(Percentile = percentile_col, phe = phe_col) %>%
    filter(l_bin < Percentile, Percentile <= u_bin)

    n     <- stratified_df %>% nrow()
    mean  <- stratified_df %>% select(phe) %>% pull() %>% mean()
    sd    <- stratified_df %>% select(phe) %>% pull() %>% sd()
    std_e <- sd / sqrt(n)
    l_err <- mean - std_e
    u_err <- mean + std_e

    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        mean   = mean,
        std_err = std_e,
        l_err = l_err,
        u_err = u_err,
        mean_str = sprintf('%.3f (%.3f,%.3f)', mean, l_err, u_err),
        bin_str = paste0(100 * l_bin, '% - ', 100 * u_bin, '%'),
        stringsAsFactors=F
    ) %>%
    mutate(mean_str = as.character(mean_str))
}

filter_by_percentile_and_count_phe <- function(df, percentile_col, phe_col, l_bin, u_bin){
    # a helper function for compute_OR.
    # This provides the counts of the descrete phenotype value (phe_col)
    # for the specified bin (l_bin, u_bin], based on the percentile of PRS (percentile_col)
    df %>%
    rename(Percentile = percentile_col, phe = phe_col) %>%
    filter(l_bin < Percentile, Percentile <= u_bin) %>%
    count(phe)
}

compute_OR <- function(df, percentile_col, phe_col, l_bin, u_bin, cnt_middle){
    # Compute the OR and sd of the trait value (phe_col), based on the
    # binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
    # The odds ratio is defined against the "middle" of the PRS distribution, and
    # we assume to have the phenotype counts in that bin (cnt_middle)
    cnt_tbl <- df %>%
    filter_by_percentile_and_count_phe(percentile_col, phe_col, l_bin, u_bin) %>%
    inner_join(cnt_middle, by='phe') %>%
    gather(bin, cnt, -phe) %>% arrange(-phe, bin)

    cnt_res <- cnt_tbl %>% mutate(cnt = as.numeric(cnt)) %>% select(cnt) %>% pull()
    names(cnt_res) <- c('n_TP', 'n_FN', 'n_FP', 'n_TN')

    OR <- (cnt_res[['n_TP']] * cnt_res[['n_TN']]) / (cnt_res[['n_FP']] * cnt_res[['n_FN']])
    LOR <- log(OR)
    se_LOR <- cnt_tbl %>% select(cnt) %>% pull() %>%
    lapply(function(x){1/x}) %>% reduce(function(x, y){x+y}) %>% sqrt()
    l_OR = exp(LOR - 1.96 * se_LOR)
    u_OR = exp(LOR + 1.96 * se_LOR)

    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        n_TP = cnt_res[['n_TP']],
        n_FN = cnt_res[['n_FN']],
        n_FP = cnt_res[['n_FP']],
        n_TN = cnt_res[['n_TN']],
        OR   = OR,
        SE_LOR = se_LOR,
        l_OR = l_OR,
        u_OR = u_OR,
        OR_str = sprintf('%.3f (%.3f,%.3f)', OR, l_OR, u_OR),
        bin_str = paste0(100 * l_bin, '% - ', 100 * u_bin, '%'),
        stringsAsFactors=F
    )
}

compute_summary_OR_df <- function(df, percentile_col, phe_col, bins=((0:10)/10)){
    cnt_middle <- df %>%
    filter_by_percentile_and_count_phe(percentile_col, phe_col, 0.4, 0.6) %>%
    rename('n_40_60' = 'n')

    1:(length(bins)-1) %>%
    lapply(function(i){
        compute_OR(df, percentile_col, phe_col, bins[i], bins[i+1], cnt_middle)
    }) %>%
    bind_rows()
}

compute_summary_mean_df <- function(df, percentile_col, phe_col, bins=c(0, .0005, .01, (1:19)/20, .99, .9995, 1)){
    1:(length(bins)-1) %>%
    lapply(function(i){
        compute_mean(df, percentile_col, phe_col, bins[i], bins[i+1])
    }) %>%
    bind_rows()
}

compute_summary_df <- function(df, percentile_col, phe_col, bins=c(0, .0005, .01, (1:19)/20, .99, .9995, 1), metric='mean'){
    if(metric == 'mean'){
        compute_summary_mean_df(df, percentile_col, phe_col, bins)
    }else if(metric == 'OR'){
        compute_summary_OR_df(df, percentile_col, phe_col, bins)
    }else{
        stop(sprintf('metric %s is not supported!', metric))
    }
}

### functions for plotting

plot_PRS_vs_phe <- function(plot_df, plot_bin2d_x=0.05, plot_bin2d_y=NULL){
    if(is.null(plot_bin2d_y)){
        plot_bin2d_y <- diff(quantile(plot_df$pheno, c(.4, .6))) / 4
    }
    plot_df %>%
    filter(
        # 99.9% coverage
        quantile(plot_df$pheno, .0005) < pheno,
        quantile(plot_df$pheno, .9995) > pheno
    ) %>%
    ggplot(aes(x = SCORE_geno_z, y = pheno)) +
    geom_bin2d(binwidth = c(plot_bin2d_x, plot_bin2d_y)) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()
}

plot_PRS_bin_vs_phe <- function(summary_plot_df, horizontal_line){
    summary_plot_df %>%
    mutate(x_ticks_labels = paste0('[', bin_str, ']')) %>%
    ggplot(aes(x=reorder(x_ticks_labels, -u_bin), y=mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = l_err, ymax = u_err)) +
    geom_hline(yintercept = horizontal_line, color='gray')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
}

## functions for additional evaluation

cat_or_zcat <- function(f){
    ifelse(endsWith(f, '.zst'), 'zstdcat', ifelse(endsWith(f, '.gz'), 'zcat', 'cat'))
}

read_PRS <- function(sscore_f){
    fread(
        cmd=paste(cat_or_zcat(sscore_f), sscore_f),
        select=c('#FID', 'IID', 'SCORE1_SUM'),
        colClasses=c('#FID'='character', 'IID'='character'),
        data.table=F
    ) %>%
    rename('FID'='#FID', 'geno_score'='SCORE1_SUM')
}

read_BETAs <- function(snpnet_BETAs_f){
    fread(snpnet_BETAs_f, colClasses=c('ID'='character'), data.table=F)
}

read_predicted_scores <- function(phe_df, sscore_f, snpnet_covar_BETAs_f, covariates){
    as.matrix(
        phe_df %>%
        mutate(ID = paste(FID, IID, sep='_')) %>%
        column_to_rownames('ID') %>%
        select(all_of(covariates))
    ) %*% as.matrix(
        read_BETAs(snpnet_covar_BETAs_f) %>%
        column_to_rownames('ID')
    ) %>%
    as.data.frame() %>%
    rownames_to_column('ID') %>%
    separate(ID, c('FID', 'IID'), sep = '_') %>%
    rename('covar_score'='BETA') %>%
    left_join(read_PRS(sscore_f), by=c('FID', 'IID')) %>%
    select(FID, IID, geno_score, covar_score)
}

fit_covar_model <- function(df, phe, covariates, family){
    # fit a regression model with covariates
    # this is useful when one wants to refit such model for external validation set
    stats::as.formula(sprintf('%s ~ 1 + %s', phe, paste(covariates, collapse =' + '))) %>%
    glm(family=family, data=df %>% mutate(phe = phe - ifelse(family=='binomial', 1, 0)))
}

fit_to_df <- function(fit){
    # convert the fit object to a data frame
    fit_df <- summary(fit)$coeff %>%
    as.data.frame() %>% rownames_to_column('variable')
    colnames(fit_df) <- c('variable', 'estimate', 'SE', 'z_or_t_value', 'P')
    fit_df
}

compute_covar_score <- function(phe_df, phe, covariates, family){
    # fit a regression model with covariates and compute covariate-only score
    as.matrix(
        phe_df %>%
        mutate(ID = paste(FID, IID, sep='_')) %>%
        column_to_rownames('ID') %>%
        select(all_of(covariates)) %>%
        as.matrix()
    ) %*% as.matrix(
        phe_df %>%
        fit_covar_model(phe, covariates, family) %>%
        fit_to_df() %>%
        filter(variable %in% covariates) %>%
        select(estimate) %>%
        as.matrix()
    ) %>%
    as.data.frame() %>%
    rownames_to_column('ID') %>%
    separate(ID, c('FID', 'IID'), sep = '_') %>%
    rename('covar_score'='estimate')
}

perform_eval <- function(response, pred, metric.type){
    if(metric.type == 'r2'){
        summary(lm(response ~ 1 + pred))$r.squared
    }else{
#         pROC::auc(pROC::roc(response, pred))
        pred.obj <- ROCR::prediction(pred, factor(response - 1))
        auc.obj <- ROCR::performance(pred.obj, measure = 'auc')
        auc.obj@y.values[[1]]
    }
}

build_eval_df <- function(phe_score_df, split_strs, metric.type){
    # build a data frame of the preformance metric
    split_strs %>%
    lapply(function(s){
        phe_score_df %>%
        filter(split == s) %>%
        mutate(geno_covar_score = geno_score + covar_score) -> sdf

        data.frame(
            split      = s,
            geno       = perform_eval(sdf$phe, sdf$geno_score, metric.type),
            covar      = perform_eval(sdf$phe, sdf$covar_score, metric.type),
            geno_covar = perform_eval(sdf$phe, sdf$geno_covar_score, metric.type),
            stringsAsFactors = F
        )
    }) %>% bind_rows()
}

compute_phe_score_df <- function(phe_df, phenotype, sscore_f, snpnet_covar_BETAs_f, covariates, family, refit_split_strs=NULL){
    # join the given dataframe of phenotype with PRS and covariate-based scores
    phe_df %>%
    read_predicted_scores(sscore_f, snpnet_covar_BETAs_f, covariates) %>%
    drop_na(geno_score, covar_score) %>%
    inner_join(
        phe_df %>% rename(!!'phe':= all_of(phenotype)) %>% select(FID, IID, phe, split, all_of(covariates)),
        by=c('FID', 'IID')
    ) %>%
    drop_na(split, phe) %>% filter(phe != -9) -> phe_score_before_refit_df

    if(is.null(refit_split_strs)){
        phe_score_before_refit_df %>%
        select(FID, IID, split, phe, geno_score, covar_score) -> phe_score_df
    }else{
        # refit covar models for the specified (typically non-WB) populations
        refit_split_strs %>%
        lapply(function(split_str){
            phe_score_before_refit_df %>%
            filter(split == split_str) %>%
            compute_covar_score('phe', covariates, family)
        }) %>% bind_rows() %>%
        select(FID, IID, covar_score) -> refit_df

        bind_rows(
            # the ones without refit
            phe_score_before_refit_df %>%
            filter(!split %in% refit_split_strs) %>%
            select(FID, IID, split, phe, geno_score, covar_score),

            # the ones from refit
            phe_score_before_refit_df %>%
            select(FID, IID, split, phe, geno_score) %>%
            inner_join(refit_df, by=c('FID', 'IID'))
        )  -> phe_score_df
    }

    phe_score_df
}

eval_performance <- function(phe_df, phenotype, sscore_f, snpnet_BETAs_f, snpnet_covar_BETAs_f, covariates, family, refit_split_strs=NULL){
    compute_phe_score_df(
        phe_df, phenotype, sscore_f, snpnet_covar_BETAs_f,
        covariates, family, refit_split_strs
    ) -> phe_score_df

    if(family=='binomial'){
        # count the number of cases and controls
        phe_score_df %>% count(split, phe) %>%
        mutate(phe = if_else(phe==2, 'case_n', 'control_n')) %>%
        spread(phe, n) %>% filter(control_n>0, case_n>0) %>%
        arrange(-case_n) -> split_cnt_df
    }else{
        # count the number of non-NA individuals
        phe_score_df %>% count(split) %>%
        arrange(-n) -> split_cnt_df
    }

    phe_score_df %>% build_eval_df((split_cnt_df %>% pull(split)), ifelse(family=='binomial', 'auc', 'r2')) %>%
    mutate(
        geno_delta = geno_covar - covar,
        phenotype_name = phenotype,
        n_variables = read_BETAs(snpnet_BETAs_f) %>% nrow()
    ) %>%
    left_join(split_cnt_df, by='split') %>%
    select(phenotype_name, split, geno, covar, geno_covar, geno_delta, n_variables, setdiff(colnames(split_cnt_df), 'split'))
}
