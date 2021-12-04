#' @importFrom magrittr %>%
#'
compute_mean <- function(df, percentile_col, phe_col, l_bin, u_bin){
    # Compute the mean and sd of the trait value (phe_col), based on the
    # binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
    stratified_df <- df %>%
    rename(!!'Percentile' := all_of(percentile_col), !!'phe' := all_of(phe_col)) %>%
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
        mean_str = sprintf("%.3f (%.3f,%.3f)", mean, l_err, u_err),
        bin_str = paste0(100 * l_bin, "% - ", 100 * u_bin, "%"),
        stringsAsFactors=F
    ) %>%
    mutate(mean_str = as.character(mean_str))
}


#' @importFrom magrittr %>%
#'
filter_by_percentile_and_count_phe <- function(df, percentile_col, phe_col, l_bin, u_bin){
    # a helper function for compute_OR.
    # This provides the counts of the descrete phenotype value (phe_col)
    # for the specified bin (l_bin, u_bin], based on the percentile of PRS (percentile_col)
    df %>%
    rename(!!'Percentile' := all_of(percentile_col), !!'phe' := all_of(phe_col)) %>%
    filter(l_bin < Percentile, Percentile <= u_bin) %>%
    count(phe) %>%
    # To cope with sparse bins where case or control counts are zero,
    # we add the following dummy counts (zeros)
    bind_rows(data.frame(
        phe=c(1, 2),
        n=as.integer(c(0,0))
    )) %>%
    group_by(phe) %>%
    summarise(n = sum(n)) %>%
    ungroup
}


#' @importFrom magrittr %>%
#'
compute_OR <- function(df, percentile_col, phe_col, l_bin, u_bin, cnt_middle){
    # Compute the OR and sd of the trait value (phe_col), based on the
    # binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
    # The odds ratio is defined against the "middle" of the PRS distribution, and
    # we assume to have the phenotype counts in that bin (cnt_middle)
    cnt_tbl <- df %>%
    filter_by_percentile_and_count_phe(percentile_col, phe_col, l_bin, u_bin) %>%
    bind_rows(cnt_middle %>% mutate(n = 0) %>% mutate(n = as.integer(n))) %>%
    group_by(phe) %>% summarise(n = sum(n), .groups="drop") %>%
    inner_join(cnt_middle, by="phe") %>%
    gather(bin, cnt, -phe) %>% arrange(-phe, bin)

    cnt_res <- cnt_tbl %>% mutate(cnt = as.numeric(cnt)) %>% select(cnt) %>% pull()
    names(cnt_res) <- c("n_TP", "n_FN", "n_FP", "n_TN")

    OR <- (cnt_res[["n_TP"]] * cnt_res[["n_TN"]]) / (cnt_res[["n_FP"]] * cnt_res[["n_FN"]])
    LOR <- log(OR)
    se_LOR <- cnt_tbl %>% select(cnt) %>% pull() %>%
    lapply(function(x){1/x}) %>% reduce(function(x, y){x+y}) %>% sqrt()
    l_OR = exp(LOR - 1.96 * se_LOR)
    u_OR = exp(LOR + 1.96 * se_LOR)

    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        n_TP = cnt_res[["n_TP"]],
        n_FN = cnt_res[["n_FN"]],
        n_FP = cnt_res[["n_FP"]],
        n_TN = cnt_res[["n_TN"]],
        OR   = OR,
        SE_LOR = se_LOR,
        l_OR = l_OR,
        u_OR = u_OR,
        OR_str = sprintf("%.3f (%.3f,%.3f)", OR, l_OR, u_OR),
        bin_str = paste0(100 * l_bin, "% - ", 100 * u_bin, "%"),
        stringsAsFactors=F
    )
}


#' @importFrom magrittr %>%
#'
compute_summary_OR_df <- function(df, percentile_col, phe_col, bins=NULL){
    if(is.null(bins)) bins <- ((0:10)/10)
    cnt_middle <- df %>%
    filter_by_percentile_and_count_phe(percentile_col, phe_col, 0.4, 0.6) %>%
    rename("n_40_60" = "n")

    1:(length(bins)-1) %>%
    lapply(function(i){
        compute_OR(df, percentile_col, phe_col, bins[i], bins[i+1], cnt_middle)
    }) %>%
    bind_rows()
}



#' @importFrom magrittr %>%
#'
compute_summary_mean_df <- function(df, percentile_col, phe_col, bins=NULL){
    if(is.null(bins)) bins <- c(0, .0005, .01, (1:19)/20, .99, .9995, 1)
    1:(length(bins)-1) %>%
    lapply(function(i){
        compute_mean(df, percentile_col, phe_col, bins[i], bins[i+1])
    }) %>%
    bind_rows()
}


compute_summary_df <- function(df, percentile_col, phe_col, bins=NULL, family="gaussian"){
    if(family == "gaussian"){
        compute_summary_mean_df(df, percentile_col, phe_col, bins)
    }else if(family == "binomial"){
        compute_summary_OR_df(df, percentile_col, phe_col, bins)
    }else{
        stop(sprintf("%s family is not supported!", family))
    }
}


#' Get r-squared value from glm fit object
#'
#' @param glm_fit A fit object from glm()
#' @return  The r-squared value
#'
glm_fit_to_r2 <- function(glm_fit){
    with(summary(glm_fit), 1 - deviance/null.deviance)
}


#' Compose a regression formula give a response variable and the predictors
#'
#'  response ~ 1 + predictor[1] + predictor[2] + ... + predictor[n]
#'
#' The function quotes variable names using quote_char
#'
#' @param response The response variable
#' @param predictors The predictor variables
#' @param quote_char A character for quotation of terms
#' @return  A string representing the regressino formula
#' @examples
#' compose_regression_formula_str('HC326', c('covar', 'PRS.score-1'))
#' compose_regression_formula_str('HC326', c('covar', 'PRS.score-1'), quote_char='')
#'
#' @export
compose_regression_formula_str <- function(response, predictors, quote_char='`'){
    return(sprintf('%s ~ 1 + %s', paste0(quote_char, response, quote_char), paste(sapply(predictors, function(term){paste0(quote_char, term, quote_char)}), collapse=' + ')))
}


#' Compute AUC with confidence interval.
#'
#' Given a data frame containing a binary response variable and predictor(s),
#' compute the AUC with confidence interval and the p-value.
#' If multiple predictors are specified, we return the minimum p-value
#'
#' @param data A data frame containing the response and predictor(s)
#' @param response The name of the response variable
#' @param predictors The predictor variable(s)
#' @param level The confidence level of the interval
#' @return a data frame containing AUC (eval), confidence interval (l_eval, u_eval), p-value along with the information of the specified response and predictors
#' @examples
#' \dontrun{
#' eval_AUC_CI(data_df, 'HC326', c('covar_score'))
#' }
#'
#' @importFrom magrittr %>%
#'
#' @export
eval_AUC_CI <- function(data, response, predictors, level=.95){
    # regression formula
    formula_str <- compose_regression_formula_str(response, predictors)

    # get the p-value
    data %>% fit_glm(formula_str, 'binomial') %>% fit_to_df() %>%
    mutate(variable = str_replace_all(variable, '`', '')) %>%
    filter(variable %in% predictors) %>% pull(P) %>%
    # we extract the smallest p-values across multiple predictors for now
    min() -> P_val

    # pROC::roc does not handle the quoted column names
    # this causes issue when predictor and response variable names contain characters like '-'
    # because it is interpreted as a minus operater not a part of the variable name.
    # To handle this, we rename column names when computing AUC-ROC with pROC::roc()
    # Note the glm() function handles such variable names as long as it's quoted.
    col_rename_func <- function(col_name){str_replace_all(col_name, '[-]', '.')}

    # call pROC::roc and pROC::ci.auc to get the AUC with CI
    pROC::ci.auc(
        pROC::roc(
            formula = stats::as.formula(
                compose_regression_formula_str(response, col_rename_func(predictors), quote_char = '')
            ),
            data = (data %>% rename_with(col_rename_func)),
            direction='<', levels=c('control'=0, 'case'=1)
        ), method='delong',  conf.level=level
    ) -> AUC_ci_list

    # format the results as a data frame
    data.frame(
        response=response,
        predictors=paste(predictors, collapse='+'),
        metric='auc',
        `eval`=AUC_ci_list[2],
        l_eval=AUC_ci_list[1],
        u_eval=AUC_ci_list[3],
        P=P_val
    )
}


#' Compute R-squared with confidence interval.
#'
#' Given a data frame containing a continupus response variable and predictor(s),
#' compute the R-squared with confidence interval and the p-value.
#' If multiple predictors are specified, we return the minimum p-value
#'
#' @param data A data frame containing the response and predictor(s)
#' @param response The name of the response variable
#' @param predictors The predictor variable(s)
#' @param level The confidence level of the interval
#' @return a data frame containing r-squared (eval), confidence interval (l_eval, u_eval), p-value along with the information of the specified response and predictors
#' @examples
#' \dontrun{
#' eval_r2_CI(data_df, 'INI50', c('covar_score'))
#' }
#'
#' @export
eval_r2_CI <- function(data, response, predictors, level=.95){
    data %>% fit_glm(
        compose_regression_formula_str(response, predictors),
        'gaussian'
    ) -> glm_fit
    # get the p-value
    glm_fit %>% fit_to_df() %>%
    mutate(variable = str_replace_all(variable, '`', '')) %>%
    filter(variable %in% predictors) %>%
    pull(P) %>%
    # we extract the smallest p-values across multiple predictors for now
    min() -> P_val
    # compute the r-squared value
    glm_fit %>% glm_fit_to_r2() -> rsq
    # call psychometric::CI.Rsq() to compute confidence interval
    # https://rdrr.io/cran/psychometric/man/CI.Rsq.html
    # https://rdrr.io/cran/psychometric/src/R/CI.Rsq.R
    CI.Rsq(rsq, n=nrow(data), k=length(predictors), level=level) %>%
    # format the resulting table
    select(-SErsq) %>% mutate(
        metric='r2',
        response=response,
        predictors=paste(predictors, collapse='+'),
        P= P_val
    ) %>%
    rename('eval'='Rsq', 'l_eval'='LCL', 'u_eval'='UCL') %>%
    select(response, predictors, metric, `eval`, l_eval, u_eval, P)
}


#' Compute R-squared or AUC with confidence interval.
#'
#' Given a data frame containing a response variable and predictor(s),
#' compute R-squared or AUC with confidence interval and the p-value.
#' We currently supports gaussian and binomial models.
#' If multiple predictors are specified, we return the minimum p-value
#'
#' @param data A data frame containing the response and predictor(s)
#' @param response The name of the response variable
#' @param predictors The predictor variable(s)
#' @param level The confidence level of the interval
#' @return a data frame containing r-squared or AUC (eval), confidence interval (l_eval, u_eval), p-value along with the information of the specified response and predictors
#' @examples
#' \dontrun{
#' eval_CI(data_df, 'INI50', c('covar_score'), 'gaussian')
#' eval_CI(data_df, 'HC326', c('covar_score'), 'binomial')
#' }
#'
#' @export
eval_CI <- function(data, response, predictors, family, level=.95){
    stopifnot(family %in% c('gaussian', 'binomial'))
    if(family == 'gaussian'){
        eval_r2_CI(data, response, predictors, level)
    }else{
        eval_AUC_CI(data, response, predictors, level)
    }
}
