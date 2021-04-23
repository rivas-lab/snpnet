# Snpnet - Efficient Lasso Solver for Large-scale SNP Data

License: GPL-2

### Reference:
  - Ruilin Li, Christopher Chang, Yosuke Tanigawa, Balasubramanian Narasimhan, Trevor Hastie, Robert Tibshirani, Manuel A. Rivas. "Fast Numerical Optimization for Genome Sequencing Data in Population Biobanks." doi: https://doi.org/10.1101/2021.02.14.431030.
  - Junyang Qian, Yosuke Tanigawa, Wenfei Du, Matthew Aguirre, Chris Chang, Robert Tibshirani, Manuel A. Rivas, Trevor Hastie. "A Fast and Scalable Framework for Large-Scale and Ultrahigh-Dimensional Sparse Regression with Application to the UK Biobank." PLOS Genetics. 16, e1009141 (2020). https://doi.org/10.1371/journal.pgen.1009141
  - Ruilin Li, Christopher Chang, Johanne Marie Justesen, Yosuke Tanigawa, Junyang Qiang, Trevor Hastie, Manuel A. Rivas, Robert Tibshirani. “Fast Lasso Method for Large-Scale and Ultrahigh-Dimensional Cox Model with Applications to UK Biobank.” Biostatistics, 2020. https://doi.org/10.1093/biostatistics/kxaa038.


### Installation:
To use this package you will need to install the following softwares and packages:
- [zstd(>=1.4.4)](https://github.com/facebook/zstd). Make sure the zstd binaries path is in the environment variable PATH, which should be automatically done if you install using [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd).
- [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/). Make sure the plink2 binary path is in the environment variable PATH.
- Most of the R dependencies can be downloaded from CRAN. Some dependencies are not available on CRAN. These R packages can be installed with:
```r
library(devtools)
install_url("https://github.com/RuilinLi/myglmnet/archive/main.zip")
install_github("chrchang/plink-ng", subdir="/2.0/cindex")
install_github("RuilinLi/plink-ng", subdir="/2.0/pgenlibr")
```
Finally, to install this package, run
```r
devtools::install_github("rivas-lab/snpnet", ref="compact")
```

### Usage
The main exported function of this package are `snpnet` and `snpnet2Base`. Both functions solve large-scale and high-dimensional regularized regression for quantitative, binary (case/control), and survival responses.
#### When should I use `snpnet`?
`snpnet` should be preferred when
- Variable screening is needed. This would be the case if the number of genetic variants is much larger than the number of individuals in the training set, and when the size of the data matrix is larger than (or close to) the size of system memory. For example, when the training data has 200,000 individuals and 1,000,000 variants, the input matrix takes about 200,000 * 1,000,000 * 2 bits = 50 GB. In this case, if your system has less than 50GB of memory then `snpnet` should be used.
#### When should I use `snpnet2Base`?
`snpnet2Base` should be preferred when
- you would like to solve a group Lasso problem, which is not available in `snpnet`
- the input data matrix is sufficiently sparse. This usualy means a large number of variants are rare variants (say MAF < 0.5%). The solver used in this function will exploit the sparsity to accelerate fitting.
- `snpnet2Base` also offers a slightly more flexible, user-defined filter of genetic variants, see below.

#### How to use `snpnet`?
The input arguments of this version of `snpnet` is essentially the same as the earlier version. See the [vignette](https://github.com/rivas-lab/snpnet/blob/compact/vignettes/vignette.pdf).

#### How to use `snpnet2Base`?
Here are the input arguments of `snpnet2Base`
- `genotype.pfile`, the prefix of PLINK2's pgen format of the genetic variants used for prediction. In particular, these three files should exist `genotype.pfile.{pgen,pvar.zst,psam}`. If you would like to utilize sparsity in this matrix, make sure that the reference allele in this file is the major allele. This can be done using Plink2 with the flag `--maj-ref`.
- `phenotype.file` the path of the file that contains the phenotype values and can be read as a table. There should be FID (family ID) and IID (individual ID) columns containing the identifier for each individual, and the phenotype column(s). (optional) some covariate columns and a column specifying the training/validation split can be included in this file.
- `phenotype` the column name of the phenotype in `phenotype.file`
- `VariantFilter`, a function that filters genetic variants based on a table of the following format (this table is generated within `snpnet2Base`)

|#CHROM |   POS|original_ID |REF |ALT | HOM_REF_CT| HET_REF_ALT_CTS| TWO_ALT_GENO_CTS| HAP_REF_CT| HAP_ALT_CTS| MISSING_CT| OBS_CT|ID            | stats_pNAs| NON_REF_CT| miss_over_non_ref| stats_means| index|
|:------|-----:|:-----------|:---|:---|----------:|---------------:|----------------:|----------:|-----------:|----------:|------:|:-------------|----------:|----------:|-----------------:|-----------:|-----:|
|1      | 69081|1:69081:G:C |G   |C   |      95048|             500|                0|          0|           0|       1035|  95548|1:69081:G:C_C |  0.0107162|        500|          2.070000|   0.0052330|     1|
|1      | 69134|1:69134:A:G |A   |G   |      96117|              55|                0|          0|           0|        411|  96172|1:69134:A:G_G |  0.0042554|         55|          7.472727|   0.0005719|     2|
|1      | 69149|1:69149:T:A |T   |A   |      96073|               6|                0|          0|           0|        504|  96079|1:69149:T:A_A |  0.0052183|          6|         84.000000|   0.0000624|     3|
|1      | 69217|1:69217:G:A |G   |A   |      96353|               1|                0|          0|           0|        229|  96354|1:69217:G:A_A |  0.0023710|          1|        229.000000|   0.0000104|     4|
|1      | 69224|1:69224:A:T |A   |T   |      96443|              17|                0|          0|           0|        123|  96460|1:69224:A:T_T |  0.0012735|         17|          7.235294|   0.0001762|     5|
|1      | 69231|1:69231:C:T |C   |T   |      96344|               1|                0|          0|           0|        238|  96345|1:69231:C:T_T |  0.0024642|          1|        238.000000|   0.0000104|     6|

`VariantFilter` (optional) should modify this data table and return the modified table. We encourage the users to look at the [default Variant filter](https://github.com/rivas-lab/snpnet/blob/a7c95cceded3bfb0881f77d92fd6a24ac17f7171/R/sparse.R#L1) and a [more advanced one](https://github.com/rivas-lab/snpnet/blob/a7c95cceded3bfb0881f77d92fd6a24ac17f7171/R/sparse.R#L301) that uses external data. The allowed operations are:
  1. Removing rows from it (based on the column values). This corresponds to exluding a genetic variant from the analysis.
  2. Adding columns to provide additional information (for example the gene symbol. This modified data table will be part of the output of `snpnet2Base`.
  3. Reordering the rows. This is necessary for group Lasso, where the variants in the same group must be in adjacent rows in this data table (see the next bullet point).

- `GroupMap` (required only for group Lasso), which is again a function (or NULL if running Lasso) takes the input from `VariantFilter` and returns an strictly increasing integer vector `out` such that. The variants in group `i` occupy the rows `out[i] + 1, ..., out[i+1]` in the output of `VariantFilter`. See [here](https://github.com/rivas-lab/snpnet/blob/a7c95cceded3bfb0881f77d92fd6a24ac17f7171/R/sparse.R#L343) for an example. 
- `family` (optional), can be one of `'gaussian', 'binomial', 'cox'`.
- `covariates`, the column names of the covariates in `phenotype.file`.
- `sparse` (optional), if set TRUE, this function will try to use a sparse representation of the input genetic matrix for the training data.
- `nlambda` (optional), number of regularization parameters to try. 
- `lambda.min.ratio` (optional), the ratio of the largest regularization parameter (so that the coefficients just become non-zero) and the smallest one.
-  `lambda` (optional), user-specified regularization parameters.
-  `split.col` (optional), column name of the variable that specifies the training/validation split in `phenotype.file`.
-  `status.col` (required for Cox model), column name of the survival status variable in `phenotype.file`.
-  `mem` (only for preprocessing), amount of memory to allocate for Plink2 for preprocessing.
-  `configs` (optional), a list with the following elements
    1. `zstdcat.path` the path where the zstdcat binary is located. Default assumes the binary is added to PATH.
    2. `gcount.full.prefix` a prefix (including path) where the intermediate data the preprocessing step generates should be stored
    3. `nCores` number of threads to use in the preprocessing step.

#### Select the number of threads
Both exported function use OpenMP to parallellize their solvers. To change the number of cores to use (say to 16) when fitting run `OMP_NUM_THREADS=16 R` or `OMP_NUM_THREADS=16 Rscript Your_script_name.R` (this will also change the number of threads of other functions that uses OpenMP in the same session). The default will use all core avaliable in the machine. On our machine we observe when the number of threads exceeds 12-16 the performance improvement plateaus, but this could vary depending on the problem size, the processor, and memory system. 
