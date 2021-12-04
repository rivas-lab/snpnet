## We keep older functions for backward compatibility

read_PRS <- function(sscore_f){
    fread(
        cmd=paste(cat_or_zcat(sscore_f), sscore_f),
        select=c("#FID", "IID", "SCORE1_SUM"),
        colClasses=c("#FID"="character", "IID"="character"),
        data.table=F
    ) %>%
    rename("FID"="#FID", "geno_score"="SCORE1_SUM")
}

read_BETAs <- function(snpnet_BETAs_f){
    fread(snpnet_BETAs_f, colClasses=c("ID"="character"), data.table=F)
}

read_predicted_scores <- function(
    phe_df, sscore_f, snpnet_covar_BETAs_f, covariates
){
    as.matrix(
        phe_df %>%
        mutate(ID = paste(FID, IID, sep="_")) %>%
        column_to_rownames("ID") %>%
        select(all_of(covariates))
    ) %*% as.matrix(
        read_BETAs(snpnet_covar_BETAs_f) %>%
        column_to_rownames("ID")
    ) %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    separate(ID, c("FID", "IID"), sep = "_") %>%
    rename("covar_score"="BETA") -> df

    if(file.exists(sscore_f)){
        df %>%
        left_join(read_PRS(sscore_f), by=c("FID", "IID")) -> df
    }else{
        # if the sscore file does not exist
        df %>% mutate(geno_score = 0) -> df
    }

    df %>%
    select(FID, IID, geno_score, covar_score)
}

fit_covar_model <- function(df, phe, covariates, family){
    # fit a regression model with covariates
    # this is useful when one wants to refit such model for external validation set
    stats::as.formula(sprintf("%s ~ 1 + %s", phe, paste(covariates, collapse =" + "))) %>%
    glm(family=family, data=df %>% mutate(phe = phe - ifelse(family=="binomial", 1, 0)))
}


compute_covar_score <- function(phe_df, phe, covariates, family){
    # fit a regression model with covariates and compute covariate-only score
    as.matrix(
        phe_df %>%
        mutate(ID = paste(FID, IID, sep="_")) %>%
        column_to_rownames("ID") %>%
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
    rownames_to_column("ID") %>%
    separate(ID, c("FID", "IID"), sep = "_") %>%
    rename("covar_score"="estimate")
}

perform_eval <- function(response, pred, metric.type) {
    if(metric.type == "r2"){
        summary(lm(response ~ 1 + pred))$r.squared
    }else{
#         pROC::auc(pROC::roc(response, pred))
        pred.obj <- ROCR::prediction(pred, factor(response - 1))
        auc.obj <- ROCR::performance(pred.obj, measure = "auc")
        auc.obj@y.values[[1]]
    }
}

build_eval_df <- function(phe_score_df, split_strs, metric.type) {
    # build a data frame of the preformance metric
    split_strs %>%
    lapply(function(s) {
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
    }) %>%
    bind_rows()
}

compute_phe_score_df <- function(
    phe_df, phenotype, sscore_f, snpnet_covar_BETAs_f,
    covariates, family, refit_split_strs = NULL
){
    # join the given dataframe of phenotype with PRS and covariate-based scores
    phe_df %>%
    read_predicted_scores(sscore_f, snpnet_covar_BETAs_f, covariates) %>%
    drop_na(geno_score, covar_score) %>%
    inner_join(
        phe_df %>%
        rename(!!"phe":= all_of(phenotype)) %>%
        select(FID, IID, phe, split, all_of(covariates)),
        by = c("FID", "IID")
    ) %>%
    drop_na(split, phe) %>% filter(phe != -9) -> phe_score_before_refit_df

    if (is.null(refit_split_strs)) {
        phe_score_before_refit_df %>%
        select(FID, IID, split, phe, geno_score, covar_score) -> phe_score_df
    } else {
        # refit covar models for the specified (typically non-WB) populations
        refit_split_strs %>%
        lapply(function(split_str){
            phe_score_before_refit_df %>%
            filter(split == split_str) %>%
            compute_covar_score("phe", covariates, family)
        }) %>%
        bind_rows() %>%
        select(FID, IID, covar_score) -> refit_df

        bind_rows(
            # the ones without refit
            phe_score_before_refit_df %>%
            filter(!split %in% refit_split_strs) %>%
            select(FID, IID, split, phe, geno_score, covar_score),

            # the ones from refit
            phe_score_before_refit_df %>%
            select(FID, IID, split, phe, geno_score) %>%
            inner_join(refit_df, by = c("FID", "IID"))
        )  -> phe_score_df
    }

    phe_score_df
}

eval_performance <- function(phe_score_df, phenotype, snpnet_BETAs_f, family) {
    if (family == "binomial") {
        # count the number of cases and controls
        phe_score_df %>%
        count(split, phe) %>%
        mutate(phe = if_else(phe == 2, "case_n", "control_n")) %>%
        spread(phe, n) %>%
        filter(control_n > 0, case_n > 0) %>%
        arrange(-case_n) -> split_cnt_df
    } else {
        # count the number of non-NA individuals
        phe_score_df %>%
        count(split) %>%
        arrange(-n) -> split_cnt_df
    }

    phe_score_df %>%
    build_eval_df(
        (split_cnt_df %>% pull(split)),
        ifelse(family == "binomial", "auc", "r2")
    ) %>%
    mutate(
        geno_delta = geno_covar - covar,
        phenotype_name = phenotype,
        n_variables = read_BETAs(snpnet_BETAs_f) %>% nrow()
    ) %>%
    left_join(split_cnt_df, by = "split") %>%
    select(
        phenotype_name, split, geno, covar, geno_covar, geno_delta,
        n_variables, setdiff(colnames(split_cnt_df), "split")
    )
}