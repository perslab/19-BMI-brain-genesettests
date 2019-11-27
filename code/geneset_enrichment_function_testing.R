#' @title geneset enrichment function tests
#' @author Jonatan Thompson, Tune Pers lab, rkm916 at ku dot dk

require("data.table")
library("here")
library("optparse")
require("Matrix")
require("magrittr")
require("parallel")
require("here")

dt_ES_full = fread("/projects/jonatan/pub-perslab/timshel-bmicelltypes2019/out/es/campbell2017_lvl2.mu.csv.gz")
vec_geneset <-readRDS("/projects/jonatan/data/genesets/BMI_rareMendelianVariants_combined.RDS")
source(here("code","geneset_enrichment_function.R"))


df_geneScore = dt_ES_full
vec_geneset = vec_geneset
testUse = "wilcoxon"
alternative="two.sided"
empPval=F
doPar = F
nCores = 0
nRep=0
vec_genesBackground = NULL
fishersCutoff = NULL
fishersTopNgenes = NULL
randomSeed = 12345


mat_wilcoxon_twosided_analyt <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = vec_geneset,
  testUse = "wilcoxon",
  alternative="two.sided",
  empPval=F,#T,
  nRep=0,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = NULL,
  doPar = F,
  nCores = 0,
  randomSeed = 12345)


# Fishers 0.75 empirical, greater
mat_fishers0.75cutoff_emp <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "fishers",
  alternative="greater",
  empPval=T,#T,
  nRep=1000L,
  vec_genesBackground = NULL,
  fishersCutoff = 0.75,
  fishersTopNgenes = NULL,
  doPar = T,
  nCores = 10L,
  randomSeed = 12345)


# Fishers 0.75 empirical, greater
mat_fishers0.75cutoff_analyt<- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "fishers",
  alternative="greater",
  empPval=F,#T,
  nRep=0,
  vec_genesBackground = NULL,
  fishersCutoff = 0.75,
  fishersTopNgenes = NULL,
  doPar = F,
  nCores = 10L,
  randomSeed = 12345)

# Fishers topNgenes analytical, greater
mat_fishers100cutoff_analyt <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "fishers",
  alternative="greater",
  empPval=F,#T,
  nRep=0,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = 100,
  doPar = F,
  nCores = 0,
  randomSeed = 12345)

# Fishers topNgenes analytical, less
mat_fishers100cutoff_analyt <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "fishers",
  alternative="less",
  empPval=F,#T,
  nRep=0,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = 100,
  doPar = F,
  nCores = 0,
  randomSeed = 12345)

# wilcoxon analytical greater
mat_wilcoxon_greater_analyt <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "wilcoxon",
  alternative="greater",
  empPval=F,#T,
  nRep=0,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = NULL,
  doPar = F,
  nCores = 0,
  randomSeed = 12345)

# wilcoxon analytical less
mat_wilcoxon_less_analyt <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "wilcoxon",
  alternative="less",
  empPval=F,#T,
  nRep=0,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = NULL,
  doPar = F,
  nCores = 0,
  randomSeed = 12345)

# wilcoxon analytical, two-sided
mat_wilcoxon_twosided_analyt <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "wilcoxon",
  alternative="two.sided",
  empPval=F,#T,
  nRep=0,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = NULL,
  doPar = F,
  nCores = 0,
  randomSeed = 12345)


# wilcoxon empirical two-sided
mat_wilcoxon_twosided_emp <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "wilcoxon",
  alternative="two.sided",
  empPval=T,#T,
  nRep=250,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = NULL,
  doPar = T,
  nCores = 20,
  randomSeed = 12345)

# t.test empirical nRep = 50, two-sided
mat_t.test_twosided_emp <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "t.test",
  alternative="two.sided",
  empPval=T,#T,
  nRep=50,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = NULL,
  doPar = T,
  nCores = 25,
  randomSeed = 12345)

# GSEA empirical
mat_GSEA_twosided_emp <- fnc_geneset_test(
  df_geneScore = dt_ES_full,
  vec_geneset = list_vec_geneset[[1]],
  testUse = "GSEA",
  alternative="two.sided",
  empPval=T,#T,
  nRep=1000,
  vec_genesBackground = NULL,
  fishersCutoff = NULL,
  fishersTopNgenes = NULL,
  doPar = T,
  nCores = 30,
  randomSeed = 12345)

## Review results

vec_mat <- vec_mat[!vec_mat %in% c("dt_ES_full", "df_geneScore")]

mat_meanCor <- sapply(vec_mat, function(mat) {
  sapply(vec_mat, function(othermat) {
    mat_cor <- cor(eval(parse(text=mat)),
                eval(parse(text=othermat)))
    #mat_cor <- mat_cor - diag(x=NA_real_,nrow=nrow(mat_cor))
    meanCor <- mean(mat_cor,na.rm = T)
    return(meanCor)
  })
})


