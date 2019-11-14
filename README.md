# Geneset tests

## Usage

Run one of a selection of statistical tests for gene set enrichment on a matrix of gene weights.

## Quick start

1. clone the repo and its submodule `git clone --recurse-submodules https://github.com/perslab/19-BMI-brain-genesettests.git`
2. go to the directory: `cd 19-BMI-brain-genesettests`
3. start R session: `R`
4. install `renv`: `ìnstall.packages("renv")`
5. install R package dependencies for the code: `renv::restore()`
6. quit R: `quit("no")`  
7. add paths to the data, test and other parameters in `call_run_geneset_tests.sh` 
8. run the analysis: `bash call_run_geneset_tests.sh`

## Output 

outputs a csv file with a row for each column in the input table and the following columns:

* statistic: if applicable, test statistics
* parameter: if applicable, parameter values (e.g. for t.test degrees of freedom. NB: For  GSEA, the edge value.)
* alternative: the parameter provided, e.g. "greater" or "two.sided".
* p.value: **unadjusted** p-values
* p.value_emp: **unadjusted** empirical p-values. If not computed, NAs

## Example

To reproduce the results from the paper, adjust the file path parameters in `calL_run_geneset_tests.sh` and leave the other parameters as they are, then 

`bash call_run_geneset_tests.sh &> geneset_tests_log_191007.txt`

## Docs

`Rscript ./code/run_geneset_tests.R --help`

## Note

* empirical p-values can be impractically slow especially for the wilcoxon and wilcoxon test 
* the GSEA function uses the [liger](https://rdrr.io/cran/liger/man/gsea.html) package, which computes p-values by permuting gene labels on the input weights. This is faster than the original GSEA algorithm but may introduce false positives when testing against genesets where some genes are co-expressed.

A [workflowr](https://github.com/jdblischak/workflowr) project.

