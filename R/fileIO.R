#' Load a lambda sequence from a R Data file
#'
#' Read a R data file (saved from a snpnet run with validation set)
#' and load the lambda sequence up to the optimal lambda index
#' This helper function is useful to for refit of the snpnet model
#' using the individuals in combined set of training and validation sets.
#'
#' @param rdata_f a file path to a R Data file
#'
#' @return lambda sequence
#'
#' @export
#'
read_lambda_sequence_from_RData_file <- function(rdata_f) {
    load(rdata_f)
    fit$full.lams[1:which.max(fit$metric.val)]
}


#' Read a configuration from a file
#'
#' Parse the 2-column config file for snpnet wrapper script
#' and return as a named list.
#'
#' We parse the covariate string using parse_covariates() function
#' and cast the data type based on the type definitions stored in
#' config_params_data_type() function.
#'
#' @param config_file a path to a config file
#'
#' @return configugarion stored in a named list
#'
read_config_from_file <- function(config_file) {
    null_strs <- c("null", "Null", "NULL", "none", "None", "NONE")

    config_df <- config_file %>% fread(header = T, sep = "\t") %>%
    setnames(c("key", "val"))

    config <- as.list(setNames(config_df$val, config_df$key))
    config_dt <- config_params_data_type()
    for (k in intersect(config_dt[["double"]], names(config))) {
        config[[k]] <- as.double(config[[k]])
    }
    for (k in intersect(config_dt[["integer"]], names(config))) {
        config[[k]] <- as.integer(config[[k]])
    }
    for (k in intersect(config_dt[["logical"]], names(config))) {
        config[[k]] <- as.logical(config[[k]])
    }
    for (key in c("status.col", "covariates", "split.col", "keep")) {
        if ((! key %in% names(config)) |
            (config[[key]] %in% null_strs) |
            is.na(config[[key]])) {
            config[[key]] <- NULL
        }
    }
    if (!"validation" %in% names(config)) config[["validation"]] <- FALSE
    if (!config[["validation"]]) config[["split.col"]] <- NULL
    config[["covariates"]] <- parse_covariates(config[["covariates"]])

    if (
        # read penalty factor from a RDS file
        ("p.factor.file" %in% names(config)) &&
        (! is.null(config[["p.factor.file"]])) &&
        (! config[["p.factor.file"]] %in% null_strs)
    ) {
        config[["p.factor"]] <- readRDS(
            config[["p.factor.file"]]
        )
    }else if (! "p.factor" %in% names(config)) {
        config[["p.factor"]] <- NULL
    }

    if (
        # for refit, read the lambda sequence from the specified RData file
        ("refit_from_RData" %in% names(config)) &&
        (! is.null(config[["refit_from_RData"]])) &&
        (! config[["refit_from_RData"]] %in% null_strs)
    ) {
        config[["lambda"]] <- read_lambda_sequence_from_RData_file(
            config[["refit_from_RData"]]
        )
        config[["split.col"]] <- NULL
    }else if (! "lambda" %in% names(config)) {
        config[["lambda"]] <- NULL
    }

    if (
        # exclude a list of SNPs
        ("exclude" %in% names(config)) &&
        (! is.null(config[["exclude"]])) &&
        (! config[["exclude"]] %in% null_strs)
    ) {
        config[["excludeSNP"]] <- get_ID_ALTs(config[["exclude"]])
    } else if (
        # extract a list of SNPs
        ("extract" %in% names(config)) &&
        (! is.null(config[["extract"]])) &&
        (! config[["extract"]] %in% null_strs)
    ) {
        config[["excludeSNP"]] <- setdiff(
            get_ID_ALTs(sprintf("%s.pvar.zst", config[["genotype.pfile"]])),
            get_ID_ALTs(config[["extract"]])
        )
    } else if (! "excludeSNP" %in% names(config)) {
        config[["excludeSNP"]] <- NULL
    }

    config
}


#' save the regression coefficients (betas) into file(s).
#'
#' We use up to 2 files corresponding to the betas for genetic variants
#' as well as the ones for covariates (optional).
#'
#' For the genetic variants, we recover the information
#' stored in the pvar[.zst] file and sort the results
#' based on CHROM and POS columns.
#'
#' @param df a dataframe containing beta values
#' @param out.file.head a prefix of the output file name
#' @param pvar a path to pvar file containing the list of variants
#' @param covariates a list of covariates
#' @param verbose A boolean variable indicating whether we should dump logs
#'
save_BETA <- function(df, out.file.head, pvar, covariates = NULL, verbose=FALSE) {
    file.geno   <- paste0(out.file.head, ".tsv")
    file.covars <- paste0(out.file.head, ".covars.tsv")
    file.timestamp <- format(Sys.time(), format="%Y%m%d-%H%M%S")

    if (! is.null(covariates)) {
        if (verbose) snpnet::snpnetLogger(
            sprintf("Saving results to %s", file.covars),
            funcname = "snpnetWrapper"
        )
        if (file.exists(file.covars)) file.rename(
            file.covars,
            paste0(out.file.head, ".", file.timestamp, ".covars.tsv")
        )
        df %>%
        filter(ID %in% covariates) %>%
        fwrite(file.covars, sep = "\t")
    }

    pvar_df <- fread(
        cmd = paste(cat_or_zcat(pvar), pvar, '| sed -e "s/^#//g"', sep = " ")
    ) %>%
    mutate(varID = paste(ID, ALT, sep = "_"))

    if (verbose) snpnet::snpnetLogger(
        sprintf("Saving results to %s", file.geno),
        funcname = "snpnetWrapper"
    )

    if (file.exists(file.geno)) file.rename(
        file.geno,
        paste0(out.file.head, ".", file.timestamp, ".tsv")
    )

    df %>%
    filter(! ID %in% covariates) %>%
    rename("varID" = "ID") %>%
    left_join(pvar_df, by = "varID") %>%
    arrange(CHROM, POS) %>%
    select(c(colnames(pvar_df %>% select(-varID)), "BETA")) %>%
    fwrite(file.geno, sep = "\t")
}


#' Read a PRS file as a data frame
#'
#'
#' Read a sscore file (output from plink2 --score function) as a data frame
#' for the specified columns along with 'FID' and 'IID'
#'
#' @param sscore_file_path A path to the sscore file
#' @param columns a list of columns to read
#' @return data frame containing phenotype data
#' @examples
#' \dontrun{
#' read_sscore_file(sscore_file_path)
#' }
#'
#' @export
read_sscore_file <- function(sscore_file_path, columns=c('SCORE1_SUM')){
    fread(
        cmd=paste(
            cat_or_zcat(sscore_file_path),
            sscore_file_path
        ),
        colClasses = c('#FID'='character', 'IID'='character'),
        select=c('#FID', 'IID', columns)
    ) %>%
    rename_with(
        function(x){str_replace(x, '#', '')}, starts_with("#")
    )
}

#' Read a phenotype file as a data frame
#'
#'
#' Read a phenotype file as a data frame for the selected columns
#' along with 'FID', 'IID', 'population', and 'split' columns
#'
#' We recode -9 as NA in population, split, age, and sex columns
#'
#' @param phenotype_file_path A path to the phenotype file
#' @param columns a list of columns to read
#' @return data frame containing phenotype data
#' @examples
#' \dontrun{
#' read_phenotype_file(phenotype_file_path, 'INI50')
#' read_phenotype_file(phenotype_file_path, c('INI50', 'HC382'))
#' read_phenotype_file(phenotype_file_path, c('age', 'sex', 'INI50', 'HC382'))
#' }
#'
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_with mutate across
#'
#' @export
read_phenotype_file <- function(phenotype_file_path, columns){
    fread(
        cmd = paste(
            cat_or_zcat(phenotype_file_path),
            phenotype_file_path
        ),
        colClasses = c("#FID"="character", "IID"="character", "population"="character", "split"="character"),
        select=c("#FID", "IID", "population", "split", columns)
    ) %>%
    dplyr::rename_with(
        function(x){str_replace(x, "#", "")}, starts_with("#")
    ) %>%
    dplyr::mutate(
        across(all_of(c("population", "split")), function(x){ifelse(x == "-9", NA, x)}),
        across(all_of(intersect(columns, c("age", "sex"))), function(x){na_if(x, -9)})
    )
}
