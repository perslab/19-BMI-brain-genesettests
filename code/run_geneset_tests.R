#' @title run geneset enrichment tests
#' @author Jonatan Thompson, Tune Pers lab, rkm916 at ku dot dk
#' @usage run the script from terminal to call fnc_geneset_test() with the parameters provided
#' @example Rscript run_geneset_tests.R --path_mat_geneScore ${MAT_ES_DIR}/${DATASET}.mu.csv.gz \
#'                                      --prefixData ${DATASET} \
#'                                      --prefixRun ${PREFIX_RUN} \
#'                                      --dirOut ${DIR_OUT} \
#'                                      --path_vec_geneset ${PATH_GENESET} \
#'                                      --testUse ${TESTUSE} \
#'                                      --alternative ${ALTERNATIVE} \
#'                                      --empPval ${EMPPVAL} \
#'                                      --doPar ${DOPAR} \
#'                                      --nCores ${NCORES} \
#'                                      --nRep ${NREP}

######################################################################
##################### INITIAL PACKAGES ###############################
######################################################################

#require("utils")
#require("devtools") # for devtools::session_info()
library("here")
library("optparse")
#require("liger") # only needed if performing GSEA test

######################################################################
####################### SOURCE FUNCTIONS #############################
######################################################################

source(here("code", "geneset_enrichment_function.R"))
source(file=here("perslab-sc-library","utility_functions.R"))

######################################################################
########################### OptParse #################################
######################################################################

option_list <- list(
  make_option("--path_mat_geneScore", type="character",
              help = "Path to matrix with genes in rows and some grouping value, such as cell type, in columns [default %default]"),
  make_option("--prefixData", type="character",
              help = "data prefix for output [default %default]"),
  make_option("--prefixRun", type="character",
              help = "run prefix for output [default %default]"),
  make_option("--dirOut", type="character", default = here("output"),
              help = "output directory [default %default]"),
  make_option("--path_vec_geneset", type="character", default=here("data", "BMI_rareMendelianVariants_combined.RDS"),
              help = "Provide path to gene sets saved as a list of character vectors in an RData or RDS object. [default %default]"),
  make_option("--testUse", type="character", default="fishers",
              help = "type of enrichment test: one of 'fishers', 't.test', 'wilcoxon', 'GSEA' [default %default]"),
  make_option("--alternative", type="character", default=c("two-sided", "greater","less"),
              help = "Which tails of the distribution to consider when computing p-values for rejecting the null hypothesis [default %default]"),
  make_option("--path_vec_genesBackground", type="character", default= NULL,
              help = "Path to a custom background. By default the script uses the row names of mat_geneScore [default %default]"),
  make_option("--fishersCutoff", type="numeric", default=NULL,
              help = "An absolute gene weight lower threshold at which to cut off genes to include in a fisher's test [default %default]"),
  make_option("--fishersTopNgenes", type="integer", default= NULL,
              help = "A gene weight rank lower threshold at which to cut off genes to include in a fisher's test [default %default]"),
  make_option("--empPval", type="logical", default= TRUE,
              help = "Use resampling to determine empirical p-values? [default %default]"),
  make_option("--doPar", type="logical", default= TRUE,
              help = "Use parallel processing (with the parallel package)? [default %default]"),
  make_option("--nCores", type="integer", default= 10L,
              help = "If using parallel processing, how many cores? [default %default]"),
  make_option("--nRep", type="integer", default= 10000,
              help = "If using resampling to determine empirical p-values, number of resampling repetitions [default %default]"),
  make_option("--randomSeed", type="integer", default= 12345L,
              help = "Set a random seed, used when generating empirical p-values [default %default]")
  )

opt <- parse_args(OptionParser(option_list=option_list))


path_mat_geneScore <- opt$path_mat_geneScore
prefixData <- opt$prefixData
prefixRun <- opt$prefixRun
dirOut <- opt$dirOut
path_vec_geneset <- opt$path_vec_geneset
testUse = opt$testUse
alternative = opt$alternative
path_vec_genesBackground <- opt$path_vec_genesBackground
fishersCutoff <- opt$fishersCutoff
fishersTopNgenes <- opt$fishersTopNgenes
empPval <- opt$empPval
doPar = opt$doPar
nCores = opt$nCores
nRep = opt$nRep
randomSeed = opt$randomSeed

######################################################################
######################### REMAINING PACKAGES #########################
######################################################################

require("Matrix")
require("magrittr")
require("parallel")

######################################################################
########################## LOAD DATA #################################
######################################################################

# load gene set(s)#
mat_geneScore <-  read.csv(path_mat_geneScore, header = T, quote = "", row.names = 1, stringsAsFactors = F)
vec_geneset <- readRDS(path_vec_geneset)
vec_genesBackground <- if (!is.null(path_vec_genesBackground)) load_obj(path_vec_genesBackground) else NULL

######################################################################
########################### RUN TEST #################################
######################################################################

message(paste0("running geneset test for ", prefixData))

df_out <- fnc_geneset_test(
  mat_geneScore = mat_geneScore,
  vec_geneset = vec_geneset,
  testUse = testUse,
  alternative=alternative,
  empPval=empPval,
  nRep=nRep,
  vec_genesBackground = vec_genesBackground,
  fishersCutoff = fishersCutoff,
  fishersTopNgenes = fishersTopNgenes,
  doPar = doPar,
  nCores = nCores,
  randomSeed = randomSeed)

df_out <- cbind(rownames(df_out), df_out)
colnames(df_out)[1] <- "cell_type"
rownames(df_out) <- NULL

######################################################################
########################### SAVE OUTPUT ##############################
######################################################################

flagDate <- substr(gsub("-","",as.character(Sys.Date())),3,1000)

saveMeta(savefnc=write.csv, x=df_out, file = gzfile(paste0(dirOut, prefixData, "_", prefixRun, "_", testUse, "_genetestOuts_", flagDate,".csv.gz")), quote = F, row.names = F)

message(paste0("geneset test done for ", prefixData), "!")
