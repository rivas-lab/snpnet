% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{predict_snpnet}
\alias{predict_snpnet}
\title{Predict from the Fitted Object or File}
\usage{
predict_snpnet(fit = NULL, saved_path = NULL, new_genotype_file, new_phenotype_file,
  phenotype, gcount_path = NULL, meta_dir = NULL, meta_suffix = ".rda", covariate_names = NULL,
  split_col = NULL, split_name = NULL, idx = NULL, family = NULL, snpnet_prefix = "output_iter_",
  snpnet_suffix = ".RData", snpnet_subdir = "results", configs = list(zstdcat.path = "zstdcat",
  zcat.path='zcat'))
}
\arguments{
\item{fit}{Fitted object returned from the snpnet function. If not specified, `saved_path` has to
be provided.}

\item{saved_path}{Path to the file that saves the fit object. The full path is constructed as
${saved_path}/${snpnet_subdir}/${snpnet_prefix}ITER${snpnet_suffix}, where ITER will be the
maximum index found in the snpnet subdirectory. If not specified, `fit` has to be provided.}

\item{new_genotype_file}{Path to the new suite of genotype files. new_genotype_file.{pgen, psam,
pvar.zst}.
must exist.}

\item{new_phenotype_file}{Path to the phenotype. The header must include FID, IID. Used for extracting covaraites and computing metrics.}

\item{phenotype}{Name of the phenotype for which the fit was computed.}

\item{gcount_path}{Path to the saved gcount file on which the meta statistics can be computed. Only if `saved_path` is specified.}

\item{meta_dir}{(Depreciated) Path to the saved meta statistics object. The full path is constructed as ${meta_dir}/${STAT}${meta_suffix}, where such files should exist for STAT = pnas, means and optionally sds. Only if `saved_path` is specified.}

\item{meta_suffix}{(Depreciated) Extension suffix of the meta statistics files. Only if `saved_path` is specified.}

\item{covariate_names}{Character vector of the names of the adjustment covariates.}

\item{split_col}{Name of the split column. If NULL, all samples will be used.}

\item{split_name}{Vector of split labels where prediction is to be made. Should be a combination of "train", "val", "test".}

\item{idx}{Vector of lambda indices on which the prediction is to be made. If not provided, will predict on all lambdas found.}

\item{family}{Type of the phenotype: "gaussian" for continuous phenotype and "binomial" for binary phenotype.}

\item{status_col}{Name of the status column for Cox model.}

\item{snpnet_prefix}{Prefix of the snpnet result files used to construct the full path. Only if `saved_path` is specified.}

\item{snpnet_suffix}{Extension suffix of the snpnet result files used to construct the full path. Only if `saved_path` is specified.}

\item{snpnet_subdir}{Name of the snpnet result subdirectory holding multiple result files for one phenotype. Only if `saved_path` is specified.}

\item{configs}{Additional list of configs including path to either zstdcat or zcat.}
}
\value{
A list containing the prediction and the resopnse for which the prediction is made.
}
\description{
Predict from the Fitted Object or File
}
