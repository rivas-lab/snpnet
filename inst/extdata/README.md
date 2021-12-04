# sample data

Here, we have example genotype and phenotype data that we use in our [vignette](/vignettes/vignette.pdf).

For genotype, we have the following 3 files:

- [`sample.pgen`](sample.pgen): this is a binary file storing the genotype matrix in a compressed format.
- [`sample.psam`](sample.psam): this is a text file containing the list of individuals
- [`sample.pvar.zst`](sample.pvar.zst): this is a [`zstd` (Zstandard)](https://facebook.github.io/zstd/)-compressed text file containing the list of genetic variants.

For phenotype (response variable, Y,  in the regression model), we have [`sample.phe`](sample.phe), which is a table stored in a text file.

Lastly, we have [`vars.rds`](vars.rds) as a R language's `.rds` format that stores a list object. One can read it with `readRDS` function in `R`.

Please check our [vignette](/vignettes/vignette.pdf) for more information.
