# Snpnet - Efficient Lasso Solver for Large-scale SNP Data

License: GPL-2

### Reference:
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
