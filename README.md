# Snpnet - Efficient Lasso Solver for Large-scale SNP Data

License: GPL-2

### Reference:
  - Junyang Qian, Yosuke Tanigawa, Wenfei Du, Matthew Aguirre, Chris Chang, Robert Tibshirani, Manuel A. Rivas, Trevor Hastie. "A Fast and Scalable Framework for Large-Scale and Ultrahigh-Dimensional Sparse Regression with Application to the UK Biobank." PLOS Genetics. 16, e1009141 (2020). https://doi.org/10.1371/journal.pgen.1009141
  - Ruilin Li, Christopher Chang, Johanne M. Justesen, Yosuke Tanigawa, Junyang Qiang, Trevor Hastie, Manuel A. Rivas, Robert Tibshirani. "Fast Lasso Method for Large-Scale and Ultrahigh-Dimensional Cox Model with Applications to UK Biobank." Biostatistics. kxaa038, 2020. https://doi.org/10.1093/biostatistics/kxaa038.

### Installation:
Most of the requirements of snpnet are available from CRAN. It also depends on the `pgenlibr`, `glmnet/glmnetPlus` and `cindex` (for survival analysis) packages. One can install them by running the following commands in R. Notice that the installation of `pgenlibr` requires [zstd(>=1.4.4)](https://github.com/facebook/zstd). It can be built from source or simply available from [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd).

```r
library(devtools)
install_github("junyangq/glmnetPlus")
install_github("chrchang/plink-ng", subdir="/2.0/cindex")
install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
```
We assume the users already have PLINK 2.0. Otherwise it can be installed from https://www.cog-genomics.org/plink/2.0/.
