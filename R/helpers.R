#' select zstdcat, zcat, or cat
#'
#' Check the suffix of the filename and return one of zstdcat, zcat, or cat
#'
#' @param f A file path.
#' @return One of the following 3 commands: zstdcat, zcat, or cat.
#' @examples
#' cat_or_zcat("file.txt")
#' cat_or_zcat("file.txt.gz")
#' cat_or_zcat("file.txt.zst")
#' \dontrun{
#' cat_or_zcat(file_path)
#' fread(cmd = paste(cat_or_zcat(file_path), file_path))
#' }
#'
#' @export
#'
cat_or_zcat <- function(filename, configs=list(zstdcat.path='zstdcat', zcat.path='zcat')){
    if(stringr::str_ends(basename(filename), '.zst')){
        return(configs[['zstdcat.path']])
    }else if(stringr::str_ends(basename(filename), '.gz')){
        return(configs[['zcat.path']])
    }else{
        return('cat')
    }
}


#' split a list string into a list
#'
#' This is a wrapper for stringr::str_split()
#'
#' @param list_str A string object
#' @return list
#' @examples
#' split_list_str('a,b,c')
#' split_list_str('a;b;c', ';')
#'
#' @importFrom stringr str_split
#' @export
split_list_str <- function(list_str, pattern=','){
    str_split(list_str, pattern)[[1]]
}


#' split a list string into a named list
#'
#' This is a wrapper for stringr::str_split()
#'
#' @param named_list_str A string object
#' @return named list
#' @examples
#' lapply(split_named_list_str('a=1,b=2,c=3'), as.numeric)
#' lapply(split_named_list_str('x=1;y=2;z=3', pattern1=';'), as.numeric)
#' lapply(split_named_list_str('p-1;q-2;r-3', pattern1=';', pattern2 = '-'), as.numeric)
#'
#' @importFrom stringr str_split
#' @export
split_named_list_str <- function(named_list_str, pattern1=',', pattern2='='){
    str_split(named_list_str, pattern1)[[1]] -> l
    setNames(
        lapply(l, function(x){(str_split(x, pattern2))[[1]][[2]]}),
        lapply(l, function(x){(str_split(x, pattern2))[[1]][[1]]})
    )
}


#' convert the fit object to a data frame
#'
#' @param fit a fit object from glm
#'
#' @export
fit_to_df <- function(fit){
    fit_df <- summary(fit)$coeff %>%
    as.data.frame() %>% rownames_to_column("variable")
    colnames(fit_df) <- c("variable", "estimate", "SE", "z_or_t_value", "P")
    fit_df
}


#' This is a helper function for read_config_from_file().
#' Here, we define the default data types so that we can properly parse
#' the parameters specified in the config file.
#'
config_params_data_type <- function() {
    list(
        double <- c(
            "missing.rate",
            "MAF.thresh",
            "glmnet.thresh",
            "lambda.min.ratio",
            "alpha",
            "p.factor"
        ),
        integer <- c(
            "nCores",
            "prevIter",
            "niter",
            "mem",
            "nlambda",
            "nlams.init",
            "nlams.delta",
            "num.snps.batch",
            "increase.size",
            "stopping.lag"
        ),
        logical <- c(
            "use.glmnetPlus",
            "standardize.variant",
            "validation",
            "early.stopping",
            "vzs",
            "save",
            "save.computeProduct",
            "verbose",
            "KKT.verbose",
            "KKT.check.aggressive.experimental",
            "rank"
        )
    )
}


#' Parse a string containing the set of covariates.
#'
#' We assume the covariate_str contains set of covariates separated with ','.
#'
#' @param covariates_str a string containing a list of covariates.
#'
#' @return List of covariates.
#'
parse_covariates <- function(covariates_str) {
    if (is.null(covariates_str) || covariates_str == "None") {
        covariates <- c()
    }else{
        covariates <- split_list_str(covariates_str, ",")
    }
    covariates
}


#' Get ID_ALT of genetic variants
#'
#' read a table (text) file that contains "ID" and "ALT" column
#' and return as a list containing "ID_ALT"
#'
#' @param ID_ALT_file a path to table (flat text) file
#'
#' @return list of ID_ALT values
#'
#' @export
get_ID_ALTs <- function(ID_ALT_file) {
    fread(
        cmd = paste(cat_or_zcat(ID_ALT_file), ID_ALT_file),
        select = c("ID", "ALT")
    ) %>%
    mutate(ID_ALT = paste0(ID, "_", ALT)) %>%
    pull(ID_ALT)
}


#' convert the snpnet regression coeffieicnets (betas) into a data frame.
#'
#' @param beta A list of vevtors containing beta values
#' at different lambda indices
#' @param lambda_idx A lambda index for which we extract the beta values
#' @param covariates A list of covariates
#' @param verbose A boolean variable indicating whether we should dump logs
#'
#' @return a data frame containing betas
#'
#' @export
snpnet_fit_to_df <- function(
    beta, lambda_idx, covariates = NULL, verbose=FALSE
) {
    if(verbose) snpnet::snpnetLogger(
        sprintf("Extracting the BETAs (lambda idx: %d)..", lambda_idx),
        funcname = "snpnetWrapper"
    )
    # extract BETAs from snpnet(glmnet) fit as a data frame
    df <- beta[lambda_idx] %>% data.frame()
    colnames(df) <- "BETA"
    df <- df %>%
    rownames_to_column("ID")

    df$BETA <- format(df$BETA, scientific = T)

    non.zero.BETAs <- union(
        covariates,
        df %>% filter(as.numeric(BETA) != 0) %>% select(ID) %>% pull()
    )
    if (verbose) snpnet::snpnetLogger(sprintf(
        "The BETAs are extracted for %d variables (%d covariates and %d genetic variants).",
        length(non.zero.BETAs), length(intersect(non.zero.BETAs, covariates)),
        length(setdiff(non.zero.BETAs, covariates))
    ), indent = 1, funcname = "snpnetWrapper")
    df %>%
    filter(ID %in% non.zero.BETAs)
}


#' Recode phenotype values
#'
#' Given a data frame of phenotypes, we apply the following operations:
#'  1) replace -9 with NA
#'  2) for binary phenotypes, we subtract 1 so that
#'     the phenotype values takes [1, 0] instead of [2, 1]
#' Note: you need to specify the list of binary and quantitative phenotypes
#' using the phes_binary and phes_quantitative args.
#'
#' @param pheno_df A phenotype data frame
#' @param phes_binary A list of binary phenotypes
#' @param phes_quantitative A list of quantitative phenotypes
#' @return an updated phenotype data drame
#' @examples
#' \dontrun{
#' recode_pheno_values(pheno_df, phes_binary=c('HC326'), phes_quantitative=c('INI50')
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate across
#'
#' @export
recode_pheno_values <- function(pheno_df, phes_binary=NULL, phes_quantitative=NULL){
    pheno_df %>%
    dplyr::mutate(
        across(c(all_of(phes_quantitative)), function(x){na_if(x, -9)}),
        across(c(all_of(phes_binary)), function(x){na_if(x, -9) - 1})
    )
}


#' Assign rownames of phenotype data frame based on FID and IID
#'
#' Given a data frame of the phenotype data, concatenate FID and IID
#' columns with a separater character (sep) and use it as a rowname
#'
#' @param pheno_df A phenotype data frame with columns named 'FID' and 'IID'
#' @param sep A separater to concatenate FID and IID
#' @return an updated phenotype data drame with rownames
#' @examples
#' \dontrun{
#' FID_IID_to_rownames(pheno_df)
#' }
#'
#' @export
FID_IID_to_rownames <- function(pheno_df, sep ="_"){
    pheno_df %>%
    mutate(FID_IID = paste(FID, IID, sep=sep)) %>%
    select(-FID, -IID) %>%
    column_to_rownames("FID_IID")
}


#' Recover FID and IID from rownames of the data frame
#'
#' Given a data frame of the phenotype data, recover FID and IID from
#' rownames
#'
#' @param pheno_df A phenotype data frame with rownames
#' @param sep A separater used to concatenate FID and IID
#' @return an updated phenotype data drame with FID and IID columns
#'
#' @examples
#' \dontrun{
#' FID_IID_from_rownames(pheno_df)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#'
#' @export
FID_IID_from_rownames <- function(pheno_df, sep="_"){
    pheno_df %>%
    rownames_to_column("FID_IID") %>%
    separate(FID_IID, c("FID", "IID"), sep=sep)
}


#' Collapse the training (train) and validation (val) sets into 'train_val'
#'
#' Given a data frame of the phenotype data, look at the 'split' column,
#' and replace 'train' and 'val' with 'train_val'
#'
#' @param pheno_df A phenotype data frame with a column named 'split'
#'
#' @return an updated phenotype data drame
#'
#' @examples
#' \dontrun{
#' update_split_column_for_refit(pheno_df)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#'
#' @export
update_split_column_for_refit <- function(pheno_df){
    pheno_df %>%
    mutate(
        split = if_else(split %in% c("train", "val"), "train_val", split)
    )
}


#' Count the number of individuals in each popluation split
#'
#' @param df A phenotype data frame with a column named 'split'
#' @return a data frame containing the number of individuals
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename count mutate
#' @importFrom tidyr replace_na spread
#'
#' @export
count_n_per_split <- function(df, pheno_col, family, case_code = 1){
    if (family == "binomial") {
        # count the number of cases and controls
        df %>%
        rename(!! "phe__" := pheno_col) %>%
        count(split, phe__) %>%
        mutate(phe__ = if_else(phe__ == case_code, "case_n", "control_n")) %>%
        spread(phe__, n) %>%
        replace_na(list(case_n = 0, control_n = 0))
    } else {
        # count the number of non-NA individuals
        df %>%
        count(split)
    }
}


#' Take a matrix product of data matrix (X_df) and coefficient matrix (beta_df)
#'
#' Given a data matrix (X_df) and a coefficient matrix (beta_df),
#' compute the matrix product.
#' You can specify the set of variables (variables) to focus
#' if you don't want to use all the data in the beta_df.
#'
#' @param X_df A data frame containing the dataset (X)
#' @param beta_df A string representing a regression formula
#' @param variables A subset of variables you'd like to take the product on
#' @param beta_estimate_cols The set of column names of the beta_df that contains the values of the coefficients.
#' @param beta_variable_col The column name of the beta_df that contains the variable name.
#' @return a data frame containing the results of XB
#' @examples
#' \dontrun{
#' compute_matrix_product(data_df, glm_fit_df, c('age', 'sex'), c('estimate'), 'variable')
#' }
#'
#' @export
compute_matrix_product <- function(
    X_df, beta_df, variables=NULL,
    beta_estimate_cols=c('estimate'),
    beta_variable_col='variable'
){
    if(is.null(variables)){
        # if set of variables are not specified,
        # we use all of the provided variables in beta_df
        beta_df %>%
        pull(all_of(beta_variable_col)) -> variables
    }

    # take the matrix product XB
    as.matrix(
        # prepater matrix X
        X_df %>% select(all_of(variables))
    ) %*% as.matrix(
        # prepare matrix B
        beta_df %>%
        rename(!!'variable_' := all_of(beta_variable_col)) %>%
        filter(variable_ %in% variables) %>%
        column_to_rownames('variable_') %>%
        select(all_of(beta_estimate_cols))
    ) %>%
    as.data.frame()
}


#' A wrapper function for glm()
#'
#' Given a data frame, a regression formula in string, and a model family,
#' call glm function.
#'
#' @param data_df A data frame containing the dataset for glm analysis
#' @param formula_str A string representing a regression formula
#' @param family The GLM family
#' @return a glm fit object
#' @examples
#' \dontrun{
#' fit_glm(data, 'response ~ 1 + x + y', 'binomial')
#' }
#'
#' @export
fit_glm <- function(data_df, formula_str, family){
    glm(stats::as.formula(formula_str), family=family, data=data_df)
}


#' Fit a specified regression model for each split independently
#' and aggregate the results into one data frame
#' @export
fit_glm_across_splits <- function(pheno_df, splits, glm_formula_str, family, split_col = "split"){
    lapply(splits, function(s){
        pheno_df %>%
        rename(!!"split__" := split_col) %>%
        filter(split__ == s) %>%
        fit_glm(glm_formula_str, family) %>%
        fit_to_df() %>%
        mutate(split = s) %>%
        select(split, variable, estimate, SE, z_or_t_value, P)
    }) %>%
    bind_rows()
}

#' @export
compute_covariate_score <- function(pheno_df, beta_df, variables, beta_cols){
    pheno_df %>%
    FID_IID_to_rownames() %>%
    compute_matrix_product(
        beta_df,
        variables,
        beta_estimate_cols = beta_cols
    ) %>%
    FID_IID_from_rownames()
}


#' @export
compute_covariate_score_across_splits <- function(
    pheno_df, covar_BETAs_df, covariates, splits, covar_score_name, split_col = "split"
){
    lapply(unique(names(splits)), function(s){

        # use the BETAs on a split specified in named list
        covar_score_split <- splits[[s]]

        # get BETAs
        covar_BETAs_df %>%
        rename(!!"split__" := split_col) %>%
        filter(split__ == covar_score_split) %>%
        rename(!!covar_score_name := 'estimate') -> covar_betas_split_df

        pheno_df %>%
        rename(!!"split__" := split_col) %>%
        filter(split__ == s) %>%
        compute_covariate_score(
            covar_betas_split_df,
            covar_betas_split_df %>%
            pull(variable) %>%
            intersect(covariates),
            c(covar_score_name)
        ) %>%
        mutate(covar_score_computed_on = covar_score_split)
    }) %>%
    bind_rows()
}


#' @export
eval_CI_across_splits <- function(
    pheno_df, predictors, pheno_col, family, splits, split_col = "split", verbose=F
){
    lapply(splits, function(split_str){
        pheno_df %>%
        rename(!!"split__" := split_col) %>%
        filter(split__ == split_str) -> filtered_df
        n_uniq_pheno_vals <- length(unique(
            pull(filtered_df, all_of(pheno_col))
        ))
        if(n_uniq_pheno_vals <= 1){
            if(verbose){
                message(sprintf(
                    " .. skip (the phenotype value is constant in %s)",
                    split_str
                ))
            }
            NULL
        }else{
            lapply(predictors, function(predictor){
                if(verbose){
                    message(sprintf('--%s %s', split_str, predictor))
                }
                tryCatch({
                    filtered_df %>%
                    eval_CI(pheno_col, c(all_of(predictor)), family) %>%
                    mutate(split = split_str)
                }, error=function(e){print(e)})
            }) %>% bind_rows()
        }
    }) %>%
    bind_rows()
}
