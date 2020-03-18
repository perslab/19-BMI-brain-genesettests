#!/bin/bash
# usage: bash call_run_geneset_tests_celltypes_vs_BMIs.sh &> log_genesettests_celltypes_vs_BMI_191007.txt


ES_DIR=/projects/timshel/sc-genetics/timshel-bmicelltypes2019/out/es/
PREFIX_RUN=run4
DIR_OUT=/projects/jonatan/test/
PATH_GENESETS=/projects/jonatan/pub-perslab/19-BMI-brain-genesettests/data/list_BMI_rareMendelianVariants_combined.RDS
TESTUSE=wilcoxon
ALTERNATIVE=greater
#two.sided
EMPPVAL=FALSE
DOPAR=FALSE
NCORES=0
NREP=0

declare -a DATASETS=("campbell2017_lvl2" "chen2017" "romanov2017" "mousebrain" "moffitt2018" "mikkelsen2019" "kimVMH2019_smartseq" "kimVMH2019_10x")

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

