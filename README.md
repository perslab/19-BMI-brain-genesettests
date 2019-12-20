# Geneset tests

## Usage

Run one of a selection of statistical tests for gene set enrichment on a matrix of gene weights.

## Quick start

1. clone the repo and its submodule `git clone --recurse-submodules https://github.com/perslab/19-BMI-brain-genesettests.git`
2. go to the directory: `cd 19-BMI-brain-genesettests`
3. start R session: `R`
4. install `renv`: `ìnstall.packages("renv")`
5. install R package dependencies for the code: `renv::restore()`. If this fails, install the packages manually using the R command `ìnstall.packages()`. See the full list of packages as well as R version under 'R session info' below.
6. quit R: `quit("no")`  
7. add paths to the data, test and other parameters in `call_run_geneset_tests_celltypes_vs_BMI.sh` 
8. run the analysis: `bash call_run_geneset_tests_celltypes_vs_BMI.sh`

## Input

Takes a gene x annotation table of weights, where the first column contains gene names and column names are celltype or similar.

## Output 

outputs a csv file with a row for each column in the input table and the following columns:

* statistic: if applicable, test statistics
* parameter: if applicable, parameter values (e.g. for `t.test` degrees of freedom. NB: For  `GSEA`, the edge value.)
* alternative: the parameter provided, e.g. `greater` or `two.sided`.
* p.value: **unadjusted** p-values
* p.value_emp: **unadjusted** empirical p-values. If not computed, NAs

## Run cell type versus rare/mendelian variant tests

To reproduce the results from the paper, adjust the file path parameters in `calL_run_geneset_tests_celltypes_vs_BMI.sh` and leave the other parameters as they are, then `bash call_run_geneset_tests_celltypes_vs_BMI.sh`

## Run cell type versus WGCNA module tests

Adjust the file path parameters in `calL_run_geneset_tests_celltypes_vs_modules.sh` and leave the other parameters as they are, then `bash call_run_geneset_tests_celltypes_vs_modules.sh`

## Arguments

`Rscript ./code/run_geneset_tests.R --help`

## Note

* empirical p-values can be impractically slow especially for the `t.test` and `wilcox.test` 
* the `GSEA` function uses the [liger](https://rdrr.io/cran/liger/man/gsea.html) package, which computes p-values by permuting gene labels on the input weights. This is faster than the original GSEA algorithm but may introduce false positives when testing against genesets where some genes are co-expressed.

## R session info

**R version 3.5.3 (2019-03-11)**

**Platform:** x86_64-pc-linux-gnu (64-bit)

**locale:**
_LC_CTYPE=en_US.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_US.UTF-8_, _LC_COLLATE=en_US.UTF-8_, _LC_MONETARY=en_US.UTF-8_, _LC_MESSAGES=en_US.UTF-8_, _LC_PAPER=en_US.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_US.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:**
_parallel_, _stats_, _graphics_, _grDevices_, _datasets_, _utils_, _methods_ and _base_

**other attached packages:**
_pander(v.0.6.3)_, _liger(v.1.0)_, _here(v.0.1)_, _optparse(v.1.6.2)_, _Matrix(v.1.2-15)_, _usethis(v.1.4.0)_, _devtools(v.2.0.1)_, _magrittr(v.1.5)_ and _workflowr(v.1.4.0)_

**loaded via a namespace (and not attached):**
_Rcpp(v.1.0.2)_, _compiler(v.3.5.3)_, _prettyunits(v.1.0.2)_, _base64enc(v.0.1-3)_, _remotes(v.2.0.2)_, _tools(v.3.5.3)_, _digest(v.0.6.20)_, _pkgbuild(v.1.0.2)_, _pkgload(v.1.0.2)_, _evaluate(v.0.14)_, _memoise(v.1.1.0)_, _lattice(v.0.20-38)_, _rlang(v.0.3.3)_, _cli(v.1.1.0)_, _xfun(v.0.8)_, _withr(v.2.1.2)_, _knitr(v.1.24)_, _desc(v.1.2.0)_, _fs(v.1.2.6)_, _rprojroot(v.1.3-2)_, _grid(v.3.5.3)_, _getopt(v.1.20.3)_, _glue(v.1.3.1)_, _R6(v.2.4.0)_, _processx(v.3.2.0)_, _rmarkdown(v.1.14)_, _sessioninfo(v.1.1.1)_, _callr(v.3.0.0)_, _backports(v.1.1.2)_, _ps(v.1.2.1)_, _htmltools(v.0.3.6)_, _assertthat(v.0.2.1)_, _renv(v.0.6.0-108)_ and _crayon(v.1.3.4)_> plibrary("pander")


A [workflowr](https://github.com/jdblischak/workflowr) project.

