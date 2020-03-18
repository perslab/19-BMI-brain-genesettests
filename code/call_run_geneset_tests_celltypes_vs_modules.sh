#!/bin/bash
# usage: bash call_run_geneset_tests_celltype_vs_modules.sh &> log_genesettests_celltype_vs_modules_191007.txt


ES_DIR=../data
PREFIX_RUN=run1
DIR_OUT=../output/
PATH_GENESETS=..data/list_vec_WGCNA_modgenes_ENSG.RDS
TESTUSE=wilcoxon
ALTERNATIVE=greater
#two.sided
EMPPVAL=FALSE
DOPAR=FALSE
NCORES=0
NREP=0

declare -a DATASETS=("mousebrain_subCELLECTprior")

for (( i=0; i<${#DATASETS[@]}; i++ )); do
	DATASET=${DATASETS[$i]}
	echo $DATASET
	time Rscript run_geneset_tests.R --path_df_geneScore ${ES_DIR}/${DATASET}.mu.csv.gz \
					 		  --prefixData ${DATASET} \
					 	 	  --prefixRun ${PREFIX_RUN} \
					 		  --dirOut ${DIR_OUT} \
				         		  --path_list_vec_genesets ${PATH_GENESETS} \
		   					  --testUse ${TESTUSE} \
					 		  --alternative ${ALTERNATIVE} \
				  	 		  --empPval ${EMPPVAL} \
					   		  --doPar ${DOPAR} \
					 		  --nCores ${NCORES} \
				 	 		  --nRep ${NREP}
done

