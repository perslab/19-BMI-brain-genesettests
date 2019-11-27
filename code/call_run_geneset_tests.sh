#!/bin/bash
# usage: bash call_run_geneset_tests.sh &> geneset_tests_log_191007.txt

ES_DIR=/projects/timshel/sc-genetics/timshel-bmicelltypes2019/out/es
PREFIX_RUN=run3
DIR_OUT=/projects/jonatan/pub-perslab/19-BMI-brain-genesettests/output/
PATH_GENESET=/projects/jonatan/pub-perslab/19-BMI-brain-genesettests/data/BMI_rareMendelianVariants_combined.RDS
TESTUSE=wilcoxon
ALTERNATIVE=two.sided
EMPPVAL=FALSE
DOPAR=FALSE
NCORES=0
NREP=0

#declare -a DATASETS=("kimVMH2019_smartseq" "kimVMH2019_10x" "campbell2017_lvl1" "campbell2017_lvl2" "chen2017" "moffitt2018" "mousebrain" "romanov2017" "tabula_muris")
declare -a DATASETS=("chen2017")
for (( i=0; i<${#DATASETS[@]}; i++ )); do
	DATASET=${DATASETS[$i]}
	echo $DATASET
	time Rscript run_geneset_tests.R --path_df_geneScore ${ES_DIR}/${DATASET}.mu.csv.gz \
					 --prefixData ${DATASET} \
					 --prefixRun ${PREFIX_RUN} \
					 --dirOut ${DIR_OUT} \
				         --path_vec_geneset ${PATH_GENESET} \
					 --testUse ${TESTUSE} \
					 --alternative ${ALTERNATIVE} \
				  	 --empPval ${EMPPVAL} \
					 --doPar ${DOPAR} \
					 --nCores ${NCORES} \
				 	 --nRep ${NREP}
done

