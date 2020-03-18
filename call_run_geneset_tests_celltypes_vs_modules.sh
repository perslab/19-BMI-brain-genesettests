#!/bin/bash
# usage: bash call_run_geneset_tests_celltype_vs_modules.sh &> log_genesettests_celltype_vs_modules_191007.txt


ES_DIR=/projects/jonatan/pub-perslab/19-BMI-brain-genesettests/data
PREFIX_RUN=wgcna_modules_testrun
DIR_OUT=/projects/jonatan/test/
PATH_GENESETS=/projects/jonatan/pub-perslab/19-BMI-brain-genesettests/data/list_vec_WGCNA_modgenes_ENSG.RDS
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
	time Rscript ./code/run_geneset_tests.R --path_df_geneScore ${ES_DIR}/${DATASET}.mu.csv.gz \
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

