#' @title run geneset enrichment tests
#' @author Jonatan Thompson, Tune Pers lab, rkm916 at ku dot dk
#' @usage run the script from terminal to call fnc_geneset_test() with the parameters provided
#' @example Rscript run_geneset_tests.R --path_df_geneScore ${MAT_ES_DIR}/${DATASET}.mu.csv.gz \
#'                                      --prefixData ${DATASET} \
#'                                      --prefixRun ${PREFIX_RUN} \
#'                                      --dirOut ${DIR_OUT} \
#'                                      --path_list_vec_genesets ${PATH_GENESET} \
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
library("broom")
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
  make_option("--path_df_geneScore", type="character",
              help = "Path to matrix with genes in rows and some grouping value, such as cell type, in columns [default %default]"),
  make_option("--prefixData", type="character",
              help = "data prefix for output [default %default]"),
  make_option("--prefixRun", type="character",
              help = "run prefix for output [default %default]"),
  make_option("--dirOut", type="character", default = here("output"),
              help = "output directory [default %default]"),
  make_option("--path_list_vec_genesets", type="character", default=here("data", "BMI_rareMendelianVariants_combined.RDS"),
              help = "Provide path to gene sets saved as a **named** list of character vectors in an RData or RDS object. [default %default]"),
  make_option("--testUse", type="character", default="fishers",
              help = "type of enrichment test: one of 'fishers', 't.test', 'wilcoxon', 'GSEA' [default %default]"),
  make_option("--alternative", type="character", default=c("two-sided", "greater","less"),
              help = "Which tails of the distribution to consider when computing p-values for rejecting the null hypothesis [default %default]"),
  make_option("--path_vec_genesBackground", type="character", default= NULL,
              help = "Path to a custom background. By default the script uses the row names of df_geneScore [default %default]"),
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


path_df_geneScore <- opt$path_df_geneScore
prefixData <- opt$prefixData
prefixRun <- opt$prefixRun
dirOut <- opt$dirOut
path_list_vec_genesets <- opt$path_list_vec_genesets
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
require("data.table")
require("Matrix")
require("magrittr")
require("parallel")
if (testUse=="GSEA") require("liger")

######################################################################
########################## LOAD DATA #################################
######################################################################

# load gene set(s)#
#df_geneScore <-  read.csv(path_df_geneScore, header = T, quote = "", row.names = 1, stringsAsFactors = F)
df_geneScore <- fread(path_df_geneScore)
list_vec_genesets <- readRDS(path_list_vec_genesets)
vec_genesBackground <- if (!is.null(path_vec_genesBackground)) load_obj(path_vec_genesBackground) else NULL

######################################################################
########################### RUN TEST #################################
######################################################################

message(paste0("running geneset test for ", prefixData))

fun <- function(vec_geneset, geneset_name) {
  df_out <- fnc_geneset_test(
    df_geneScore = df_geneScore,
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

  dt_out <- data.table(
    "geneset_name" = geneset_name,
    "cell_type"=rownames(df_out),
    df_out)
}

list_dt_out <-  if (!doPar) {
  safeParallel(fun = fun, list_iterable = list("vec_geneset" = list_vec_genesets,
                                               geneset_name = names(list_vec_genesets)),
               simplify = F,
               timeout = 600,
               n_cores = 20)
} else {
  mapply(FUN = fun,
         vec_geneset = list_vec_genesets,
         geneset_name = names(list_vec_genesets),
         SIMPLIFY=F)
}

if (length(list_dt_out)>1) {
  dt_out <- Reduce(x=list_dt_out, rbind.data.frame)
} else {
  dt_out <- list_dt_out[[1]]
}
#colnames(df_out)[1] <- "cell_type"
#rownames(df_out) <- NULL

######################################################################
########################### SAVE OUTPUT ##############################
######################################################################

flagDate <- substr(gsub("-","",as.character(Sys.Date())),3,1000)

saveMeta(savefnc=fwrite, x=dt_out, file = paste0(dirOut, prefixData, "_", prefixRun, "_", testUse, "_genetestOuts_", flagDate,".csv"))
system2(command = "gzip", args = c("-f", paste0(dirOut, prefixData, "_", prefixRun, "_", testUse, "_genetestOuts_", flagDate,".csv")))

message(paste0("geneset test done for ", prefixData), "!")
