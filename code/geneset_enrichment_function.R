#' @title geneset test function
#' @author Jonatan Thompson, Tune Pers lab, rkm916 at ku dot dk

# define geneset test function
fnc_geneset_test <- function(mat_geneScore,
                            vec_geneset,
                            testUse,
                            alternative = c("two.sided", "less", "greater"),
                            vec_genesBackground = NULL,
                            fishersCutoff = NULL,
                            fishersTopNgenes = NULL,
                            empPval=F,
                            nRep=1000L,
                            doPar = T,
                            nCores = 10L,
                            randomSeed = 12345L) {
#' @usage Run one of a selection of statistical tests for gene set enrichment on a matrix of gene weights.
#' Generate empirical p-values.
#' @example list_genesetResults <- fnc_geneset_test(mat_geneScore = mat_geneScore,
#'                                        vec_geneset = vec_geneset,
#'                                        testUse = "fishers",
#'                                        alternative = "greater",
#'                                        empPval = TRUE,
#'                                        nRep = 1e4,
#'                                        fishersCutoff = 0.75,
#'                                        doPar = TRUE,
#'                                        nCores = 10,
#'                                        randomSeed = 12345,
#'                                        vec_genesBackground = NULL)
#' @param mat_geneScore Gene score matrix, or type that can be coerced to a matrix, with gene names in the first column and grouping variable column names
#' @param vec_geneset geneset, human (ensembl or symbol). Character vector
#' @param testUse geneset test to use, one of 'fishers', 't.test', 'wilcoxon', 'GSEA'.
#' @param alternative which tails to return in the statistical test. Must be one of "two.sided", "less", "greater". "two.sided" sums the tails above the
#' value of the statistic and below minus the absolute value of the statistic
#' @param empPval compute empirical p-values? Boolean, default if (testUse==GSEA) TRUE else FALSE
#' @param nRep number of samples of NULL distribution for empirical p-values. Integer, default 10000L
#' @param vec_genesBackground custom gene background, human (ensembl or symbol). Character vector. If NULL (default) use all mat_geneScore rownames
#' @param fishersCutoff numeric. if testUse=="fishers", fishersCutoff sets the cutoff for gene scores to include.  If fishersTopNgenes is given, must be NULL
#' @param fishersTopNgenes integer if testUse=="fishers", fishersTopNgenes sets the number of genes to include after sorting by descending gene score. If fishersCutoff is given, must be NULL
#' @param doPar use parallel computation? Uses 'FORK' clusters from the parallel package. Highly recommended for wilcoxon and GSEA in particular. Boolean, default TRUE.
#' @param nCores how many cores to use for parallel computation. Integer, default 10
#' @param randomSeed passed to set.seed. Integer, defaults to 12345L
#' @return data.frame containing
#'  statistic: if applicable, test statistics
#'  parameter: if applicable, parameter values (e.g. for t.test degrees of freedom. NB: For  GSEA, the edge value.)
#'  alternative: the value provided
#'  p.value: **unadjusted** p-values
#'  p.value_emp: **unadjusted** empirical p-values. If not computed, NAs
#'
#' @note the Gene Set Enrichment Analysis (GSEA) uses liger, which generates p-values by permuting gene names to generate null enrichment scores.
#' This method doesn't address gene co-expression. Hence the method will generate false positives. The original GSEA paper avoids this by permuting sample condition labels, which isn't always applicable https://www.pnas.org/content/102/43/15545

  set.seed(seed=randomSeed)

  ######################################################################
  ######################## CHECK INPUTS ################################
  ######################################################################

  stopifnot(!is.na(as.matrix(mat_geneScore)),
            typeof(vec_geneset)=="character",
            testUse %in% c('fishers', 't.test', 'wilcoxon', 'GSEA'),
            length(alternative)==1,
            typeof(fishersCutoff) %in% c("NULL", "integer", "numeric", "double"),
            typeof(fishersTopNgenes) %in% c("NULL", "integer", "numeric", "double"),
            typeof(empPval) == "logical",
            typeof(doPar) =="logical",
            !is.na(as.integer(nCores)),
            !is.na(as.integer(nRep)),
            nCores >= 0,
            !is.na(as.integer(randomSeed)))

  stopifnot(alternative %in% c("greater", "less", "two.sided"))
  if (!is.null(fishersTopNgenes)) stopifnot(fishersTopNgenes <= nrow(mat_geneScore))

  if (testUse=="fishers") stopifnot(!is.null(fishersCutoff) | !is.null(fishersTopNgenes))
  if (testUse=="fishers") if (alternative=="two.sided") stop("Fisher's test requires alternative to be 'less' or 'greater'")

  if (testUse=="GSEA") {
    if (empPval ==F) {
      empPval <- TRUE
      warning("GSEA test always uses empirical p-values")
    }
    if (alternative!="two.sided") {
      warning("The liger implementation of GSEA always generates a two-sided p-value.")
    }
  }

  ######################################################################
  ####################### CONVERT INPUTS ###############################
  ######################################################################

  # get background genes
  if (is.null(vec_genesBackground)) vec_genesBackground <- rownames(mat_geneScore)

  # remove NAs in geneset
  vec_geneset <- vec_geneset[!is.na(vec_geneset)]

  # CONVERT GENE WEIGHT MATRIX TO LIST OF NAMED VECTORS
  list_vec_geneWeight <- lapply("X"=mat_geneScore, FUN = function(vec_geneWeight) {
    names(vec_geneWeight) <- rownames(mat_geneScore)
    vec_geneWeight <- sort(vec_geneWeight, decreasing = T)
    return(vec_geneWeight)
  })

  # Fisher's: cut off gene weight vectors
  if (testUse=="fishers") {
    if (!is.null(fishersCutoff)) {
      list_vec_geneWeight <- lapply(list_vec_geneWeight, FUN = function(vec_geneWeight) vec_geneWeight[vec_geneWeight >= fishersCutoff])
    } else if (!is.null(fishersTopNgenes)) {
      list_vec_geneWeight <- lapply(list_vec_geneWeight, FUN = function(vec_geneWeight) vec_geneWeight[1:fishersTopNgenes])
    }
  }

  ######################################################################
  #################### DEFINE TEST P-VALUE FUNCTION ####################
  ######################################################################

  if (testUse=="fishers") {
    fnc_test <- function(vec_geneWeight,
                             vec_geneset,
                             vec_genesBackground) {
      m = length(vec_geneWeight) # number of white balls in the urn
      n = length(vec_genesBackground)-m  # number of black balls in the urn
      k = length(vec_geneset) # number of balls drawn from the urn
      q = length(intersect(vec_geneset, names(vec_geneWeight))) # number of white balls drawn without replacement from the urn
      pval <- tryCatch({
        phyper(q=q, m=m, n=n, k=k,
                       lower.tail= if (alternative=="less") T else if (alternative=="greater") F)},
        error=function(err1) {
          c("statistic"=NA_real_, "p.value"=NA_real_, "parameter"=NA_integer_)
        }
        ) # returns a p-value
      return(c("statistic"=NA_real_, "parameter"=NA_integer_,"p.value"=pval))
    }
  } else if (testUse == "t.test") {
    fnc_test <- function(vec_geneWeight,
                             vec_geneset,
                             vec_genesBackground) {
      vec_logicalGeneset <- names(vec_geneWeight) %in% vec_geneset
      testout <- tryCatch({
      t.test(x=vec_geneWeight[vec_logicalGeneset],
            y=vec_geneWeight[!vec_logicalGeneset],
            alternative= alternative,
            paired=F,
            conf.level = 0.95) # returns t.test object
      }, error= function(err1) {
        c("statistic"=NA_real_, "p.value"=NA_real_, "parameter"=NA_integer_)
      })
    return(c("statistic"=as.numeric(testout[["statistic"]]), "parameter"=testout[["parameter"]],"p.value"= testout[["p.value"]]))
    }
  } else if (testUse == "wilcoxon") {
    fnc_test <- function(vec_geneWeight,
                             vec_geneset,
                             vec_genesBackground) {
      vec_logicalGeneset <- names(vec_geneWeight) %in% vec_geneset
      testout <- tryCatch({
      wilcox.test(x=vec_geneWeight[vec_logicalGeneset],
                  y=vec_geneWeight[!vec_logicalGeneset],
                  alternative= alternative,
                  conf.level = 0.95) # returns wilcoxon.test obj
    }, error= function(err1) {
      c("statistic"=NA_real_, "p.value"=NA_real_, "parameter"=NA_integer_)
    })
      return(c("statistic"=as.numeric(testout[["statistic"]]),"parameter"=NA_real_, "p.value"=testout[["p.value"]]))
    }
  } else if (testUse =="GSEA") {
    fnc_test <- function(vec_geneWeight,
                             vec_geneset,
                             vec_genesBackground){
      gseaout <- tryCatch({
        # we use bulk gsea on a single geneset because the liger::gsea function
        # just returns a p-value irrespective of whether sscore is positive or neg
        mat_gsea <- bulk.gsea(values = vec_geneWeight,
                              set.list = list(vec_geneset),
                              rank = F,
                              n.rand = nRep-1,
                              mc.cores = if (!doPar) 1 else nCores)
        }, error = function(err1) {
          c("statistic"=NA_real_, "p.value"=NA_real_, "parameter"=NA_integer_)
        })
      return(c("statistic"=mat_gsea[,"sscore"],
               "parameter"=mat_gsea[,"edge"],
               "p.value"=mat_gsea[,"p.val"]))
    }
  }

  ######################################################################
  ################## GENERATE ANALYTICAL P VALUES ######################
  ######################################################################

  sapply(list_vec_geneWeight, function(vec_geneWeight) {
    fnc_test(vec_geneset=vec_geneset,
                 vec_geneWeight=vec_geneWeight,
                 vec_genesBackground=vec_genesBackground)}) %>%
    t %>%
    data.frame ->
    df_testOut

  ######################################################################
  ######################## GENERATE NULL P VAL #########################
  ######################################################################

  if (doPar & nCores>1) {
    fnc <- "parSapply"
    cl <- makeCluster(spec=max(1,nCores),
                      type="FORK",
                      timeout=120)
  } else {
    fnc <- "sapply"
  }

  if (empPval & testUse == "fishers") {
    # generate null p.values and test the obtained p.value against them
    df_testOut[["p.value_emp"]] <- parMapply(FUN=function(vec_geneWeight, pVal){
      # generate shuffled ES vector replicates
      list_vec_geneWeightrand <- lapply(1:(nRep-1), function(i) {
        vec_geneWeightrand <- vec_geneWeight
        names(vec_geneWeightrand) <- sample(x = vec_genesBackground,size = length(vec_geneWeight),replace = F)
        vec_geneWeightrand
      })

      fun = function(vec_geneset,
               vec_geneWeight,
               vec_genesBackground) {
        vec_out <- fnc_test(vec_geneset=vec_geneset,
                 vec_geneWeight=vec_geneWeight,
                 vec_genesBackground=vec_genesBackground)
        return(vec_out["p.value"])
      }

      # set up computation, parallel or vectorized
      list_args <- list("X" = list_vec_geneWeightrand,
                        "FUN" = fun,
                        "vec_geneset" = vec_geneset,
                        "vec_genesBackground" = vec_genesBackground,
                        "simplify"=T)

      if (doPar & nCores>1) {
        # fnc <- "parSapply"
        # cl <- makeCluster(spec=max(1,nCores),
        #                   type="FORK",
        #                   timeout=45)
        list_args[["cl"]] <- cl
      } #else {
      #  fnc <- "sapply"
      #}

      # compute NULL p values
      vec_pValNULL <- tryCatch({
        do.call(what=fnc, args = list_args) #suppressWarnings(do.call(what=fnc, args = list_args))
      }, error = function(err) {
        message(paste0("the following error occurred during parallel computation: ", err, ", switching to non-parallel"))
        list_args <- list_args[!names(list_args)=="cl"]
        fnc = "sapply"
        doPar = F
        nRep = 0
        do.call(what=fnc, args = list_args)
        })

      # if (doPar) stopCluster(cl)

      # return empirical p value
      # NB: pVal is already based on greater/less/two.sided. By checking the empirical p-value of getting such a p-value, we are just using one tail.
      pVal_out <- (sum(vec_pValNULL<pVal)+1) / (nRep)

      return(pVal_out)
    },
    vec_geneWeight = list_vec_geneWeight,
    pVal = df_testOut[["p.value"]],
    SIMPLIFY=T)

  } else if (empPval & testUse %in% c('t.test', 'wilcoxon')) {

    df_testOut[["p.value_emp"]] <- mapply(FUN=function(vec_geneWeight, statistic){
      # generate shuffled ES vector replicates
      list_vec_geneWeightrand <- lapply(1:(nRep-1), function(i) {
        vec_geneWeightrand <- vec_geneWeight
        names(vec_geneWeightrand) <- sample(x = vec_genesBackground,size = length(vec_geneWeight),replace = F)
        vec_geneWeightrand
      })

      fun = function(vec_geneset,
                     vec_geneWeight,
                     vec_genesBackground) {
        vec_out <- fnc_test(vec_geneset=vec_geneset,
                            vec_geneWeight=vec_geneWeight,
                            vec_genesBackground=vec_genesBackground)
        return(vec_out["statistic"])
      }

      # set up computation, parallel or vectorized
      list_args <- list("X" = list_vec_geneWeightrand,
                        "FUN" = fun,
                        "vec_geneset" = vec_geneset,
                        "vec_genesBackground" = vec_genesBackground,
                        "simplify"=T)

      if (doPar & nCores>1) {
        # fnc <- "parSapply"
        # cl <- makeCluster(spec=max(1,nCores),
        #           type="FORK",
        #           timeout=45)
        list_args[["cl"]] <- cl
      } #else {
      #  fnc <- "sapply"
      #}

      vec_statisticNULL <- tryCatch({
        do.call(what=fnc, args = list_args) #suppressWarnings(do.call(what=fnc, args = list_args))
      }, error = function(err) {
        message(paste0("the following error occurred during parallel computation: ", err, ", switching to non-parallel"))
        list_args <- list_args[!names(list_args)=="cl"]
        fnc = "sapply"
        doPar = F
        nRep = 0
        do.call(what=fnc, args = list_args)
      })

      # if (doPar) stopCluster(cl)

      # return empirical p value
      # NB: pVal is already based on greater/less/two.sided. By checking the empirical p-value of getting such a p-value, we are just using one tail.
      if (alternative == "greater") {
        pVal_out <-  (sum(vec_statisticNULL > statistic)+1) / (nRep)
      } else if (alternative == "less") {
        pVal_out <- (sum(vec_statisticNULL < statistic)+1) / (nRep)
      } else if (alternative == "two.sided") {
          # instead of adjusting the signifince threshold
          # we add up the integrals above the 0.05 threshold from both tails
        pVal_out <- (sum(vec_statisticNULL < -abs(statistic)) + sum(vec_statisticNULL > abs(statistic))) / nRep
        }
      return(pVal_out)
    },
    vec_geneWeight = list_vec_geneWeight,
    statistic = df_testOut[,"statistic"],
    SIMPLIFY=T)

  } else if (testUse  == "GSEA") {
    df_testOut[["p.value_emp"]] <- df_testOut[["p.value"]]
  } else if (!empPval) {
    df_testOut[["p.value_emp"]] <-  NA_real_
  }

  if (doPar & nCores>1) try(stopCluster(cl))

  df_testOut$alternative <- alternative

  return(df_testOut)

}
