#!/bin/bash
############################################################################
# Script : CN-Learn_pipeline.sh                                   #
# Author : Renjie Tan                                                      #
# Date   : 5/22/2020                                                       #
#                                                                          #
# This script is part of the XHMM pipeline. It uses the read depth info    #
# extracted in the script xhmm_extract.sh to predict CNVs.                 #
#                                                                          #
# (c) 2020 - Renjie Tan                                                    #
# Licenced under the GNU General Public License 3.0.                       #
############################################################################

############################################
# STEP 0: Declare variables and file names #
############################################
DATA_LOGS_DIR='/home/rt2776/SPARK/temp/'
GATK_SW_DIR='/home/rt2776/softwares/gatk-3.5/'
TARGET_PROBES='/home/rt2776/source/capture_kits/xgen_plus_spikein.b38.bed'
REF_GENOME='/share/data/RGC_b38/genome.hg38rg.fa'

##############################################################################################
# STEP1 : Limited to 27270 samples for downstream analysis
##############################################################################################
data_path='/home/rt2776/CN_Learn/data5_spark27k_rare/'
cd ${data_path}
cat canoes_calls.csv |grep -w -f /home/rt2776/CN_Learn/source/spark_sample_27270.txt >canoes_calls_27k.csv
cat clamms_calls.txt |grep -w -f /home/rt2776/CN_Learn/source/spark_sample_27270.txt>clamms_calls_27k.txt
cat xhmm_calls.txt |grep -w -f /home/rt2776/CN_Learn/source/spark_sample_27270.txt>xhmm_calls_27k.txt

cat ${data_path}canoes_calls_27k.csv |cut -f 2 -d " "|sort|uniq|wc
cat ${data_path}clamms_calls_27k.txt |cut -f 5|sort|uniq|wc
cat ${data_path}xhmm_calls_27k.txt |cut -f 1|sort|uniq|wc

python /home/rt2776/CN_Learn/cnv_analysis/scripts/cn_learn_sample_checking.py

##############################################################################################
# STEP 1: Get rare CNVs for each callset
#         Annotate CNV frequency in the current cohort and filter >1%
##############################################################################################
data_path='/home/rt2776/CN_Learn/data5_spark27k_rare/'
qsub -t 1-7272 /home/rt2776/cnv_analysis/scripts/cluster_af_give_cohort.sh \
     ${data_path}xhmm_calls_27k.txt \
     ${data_path}xhmm_calls_27k.txt \
     ${data_path}af_cohort/xhmm_af/xhmm_calls_af.cnv
cat ${data_path}af_cohort/xhmm_af/* >${data_path}af_cohort/xhmm_calls_af.cnv
python /home/rt2776/cnv_analysis/scripts/cluster_results_checking.py

qsub -t 1-9047 /home/rt2776/cnv_analysis/scripts/cluster_af_give_cohort.sh \
     ${data_path}clamms_calls_27k.txt \
     ${data_path}clamms_calls_27k.txt \
     ${data_path}af_cohort/clamms_af/clamms_calls_af.cnv 
cat ${data_path}af_cohort/clamms_af/* >${data_path}af_cohort/clamms_calls_af.cnv
python /home/rt2776/cnv_analysis/scripts/cluster_results_checking.py

python ~/cnv_analysis/scripts/5_annotation.py cnvfrequency_given_cohort \
    --input ${data_path}canoes_calls_27k.csv \
    --given_cohort ${data_path}canoes_calls_27k.csv \
    --output ${data_path}af_cohort/canoes_calls_af.cnv

## Filtering af >1% in cohort
NUM_SAMPLES=`cat ${data_path}clamms_calls_27k.txt |awk '{print $5}'|sort|uniq|wc -l`
THRESH_NUM=`echo "0.01 * $NUM_SAMPLES"|bc|awk '{x = $1; if (x != int(x)) {x = int(x)+1} print (x-1)}'`
cat ${data_path}af_cohort/clamms_calls_af.cnv |awk -F '\t' \
    '{if($20 == "#Carriers(inThisCohort)" || $20 <= '$THRESH_NUM') print $0}' >${data_path}clamms_calls_rare.cnv

NUM_SAMPLES=`cat ${data_path}xhmm_calls_27k.txt |awk '{print $1}'|sort|uniq|wc -l`
THRESH_NUM=`echo "0.01 * $NUM_SAMPLES"|bc|awk '{x = $1; if (x != int(x)) {x = int(x)+1} print (x-1)}'`
cat ${data_path}af_cohort/xhmm_calls_af.cnv |awk -F '\t' \
    '{if($17 == "#Carriers(inThisCohort)" || $17 <= '$THRESH_NUM') print $0}' >${data_path}xhmm_calls_rare.cnv

NUM_SAMPLES=`cat ${data_path}af_cohort/canoes_calls_af.cnv |cut -f2|sort|uniq|wc -l`
THRESH_NUM=`echo "0.01 * $NUM_SAMPLES"|bc|awk '{x = $1; if (x != int(x)) {x = int(x)+1} print (x-1)}'`
cat ${data_path}af_cohort/canoes_calls_af.cnv |awk -F '\t' \
    '{if($15 == "#Carriers(inThisCohort)" || $15 <= '$THRESH_NUM') print $0}' >${data_path}canoes_calls_rare.cnv


##############################################################################################
# STEP : CN-Learn pipeline
##############################################################################################
script_path='/home/rt2776/CN_Learn/scripts/'
## Step 5 | Measure the overlap among callers
## Run calculate_CNV_overlap.sh to measure the CNV overlap among all the callers used.
bash ${script_path}calculate_CNV_overlap.sh
#Job started on compute-0-83.local at Tue May  5 23:49:55 EDT 2020
# Wed May  6 00:30:57 EDT 2020

## Step 6A | Extract basepair level coverage info
## Run extract_bp_coverage.sh to extract the basepair level coverage for each sample. Since this information can be extracted independently for each sample, make the necessary changes to this script to parallelize the process.
bash ${script_path}extract_bp_coverage.sh


## Step 6B | Resolve breakpoints
## Run merge_overlapping_CNVs_readdepth.sh to resolve breakpoint conflicts of concordant CNVs.
bash ${script_path}merge_overlapping_CNVs_readdepth_step1_prechecks_RT.sh
qsub -t 1-27270 ${script_path}merge_overlapping_CNVs_readdepth_step2_core_cluster_RT.sh
#check the results and process the missing samples after cluster running
bash ${script_path}merge_overlapping_CNVs_readdepth_step2_core_singleRun_RT.sh 1 27270
bash ${script_path}merge_overlapping_CNVs_readdepth_step2_core_sample_checking.sh 1 27270

bash ${script_path}merge_overlapping_CNVs_readdepth_step34_generate_results_RT.sh

##Step 7 | Extract GC content and mappability in breakpoint-resolved CNV regions
bash ${script_path}extract_gc_map_vals.sh

##Step 8 | Label CNVs based on gold-standard validations
#goto
cn_learn_validated_cnv_selection.sh

# modify the file location and run
python ${script_path}rt_lable_cnvs_step1_transfer_format.py
python ${script_path}rt_lable_cnvs_step2_annotate_with_gsd.py

##Step 9 | Classify CNVs
# divid predicted CNVs into three groups. Samples in the same family will send to the same group.
data_path='/home/rt2776/CN_Learn/data5_spark27k_rare/classify_cnvs_cross_validation/'
cat ${data_path}final_preds_with_lable.txt |grep -w -f ${data_path}spark_28k_group1.txt >${data_path}final_preds_with_lable_group1.txt
cat ${data_path}final_preds_with_lable.txt |grep -w -f ${data_path}spark_28k_group2.txt >${data_path}final_preds_with_lable_group2.txt
cat ${data_path}final_preds_with_lable.txt |grep -w -f ${data_path}spark_28k_group3.txt >${data_path}final_preds_with_lable_group3.txt

#group1
cat ${data_path}header ${data_path}final_preds_with_lable_group1.txt ${data_path}final_preds_with_lable_group2.txt| \
    awk '{if($16=="LABEL_VAL"||$16==0||$16==1) print $0}'>${data_path}training12_test3/training_data.txt
cat ${data_path}header ${data_path}final_preds_with_lable_group3.txt |cut -f1-15>${data_path}training12_test3/test_data.txt
sed -i '/^[Y|X]/d' ${data_path}training12_test3/training_data.txt
sed -i '/^[Y|X]/d' ${data_path}training12_test3/test_data.txt

#group2
cat ${data_path}header ${data_path}final_preds_with_lable_group1.txt ${data_path}final_preds_with_lable_group3.txt| \
    awk '{if($16=="LABEL_VAL"||$16==0||$16==1) print $0}'>${data_path}training13_test2/training_data.txt
cat ${data_path}header ${data_path}final_preds_with_lable_group2.txt |cut -f1-15>${data_path}training13_test2/test_data.txt
sed -i '/^[Y|X]/d' ${data_path}training13_test2/training_data.txt
sed -i '/^[Y|X]/d' ${data_path}training13_test2/test_data.txt

#group3
cat ${data_path}header ${data_path}final_preds_with_lable_group2.txt ${data_path}final_preds_with_lable_group3.txt| \
    awk '{if($16=="LABEL_VAL"||$16==0||$16==1) print $0}'>${data_path}training23_test1/training_data.txt
cat ${data_path}header ${data_path}final_preds_with_lable_group1.txt |cut -f1-15>${data_path}training23_test1/test_data.txt
sed -i '/^[Y|X]/d' ${data_path}training23_test1/training_data.txt
sed -i '/^[Y|X]/d' ${data_path}training23_test1/test_data.txt

##########################################
# STEP 2: Run CN-Learn for each strategy #
##########################################
source /home/rt2776/CN_Learn/config.params
LOWER_SIZE_LIMIT=1
UPPER_SIZE_LIMIT=50000000

GROUP_DATA_DIR=${data_path}'training12_test3/'
python -u ${SCRIPTS_DIR}cn_learn.py  ${GROUP_DATA_DIR}  training_data.txt  test_data.txt \
                                ${CLASSIFIER}  ${LOWER_SIZE_LIMIT}  ${UPPER_SIZE_LIMIT}  ${NUM_TREES} \
                                                            ${CALLER_COUNT}  ${CALLER_LIST}

GROUP_DATA_DIR=${data_path}'training13_test2/'
python -u ${SCRIPTS_DIR}cn_learn.py  ${GROUP_DATA_DIR}  training_data.txt  test_data.txt \
                                ${CLASSIFIER}  ${LOWER_SIZE_LIMIT}  ${UPPER_SIZE_LIMIT}  ${NUM_TREES} \
                                                            ${CALLER_COUNT}  ${CALLER_LIST}

GROUP_DATA_DIR=${data_path}'training23_test1/'
python -u ${SCRIPTS_DIR}cn_learn.py  ${GROUP_DATA_DIR}  training_data.txt  test_data.txt \
                                ${CLASSIFIER}  ${LOWER_SIZE_LIMIT}  ${UPPER_SIZE_LIMIT}  ${NUM_TREES} \
                                                            ${CALLER_COUNT}  ${CALLER_LIST}
## Combine Results
cat ${data_path}'training12_test3/'CNV_list_with_predictions.csv \
    ${data_path}'training13_test2/'CNV_list_with_predictions.csv \
    ${data_path}'training23_test1/'CNV_list_with_predictions.csv \
    >${data_path}CNV_list_with_predictions.csv 
delete two redundant top rows

##########################################
# STEP  For the CNVs in offspring, get quality scores in parents
##########################################
script_path='/home/rt2776/CN_Learn/scripts/'
data_path='/home/rt2776/CN_Learn/data5_spark27k_rare/'
python ${script_path}post_cn_learn.py prepare_trio_scores \
    --input ${data_path}'classify_cnvs_cross_validation/training12_test3/'CNV_list_with_predictions.txt \
    --consolidated_calls ${data_path}'consolidated_calls.bed' \
    --output ${data_path}'classify_cnvs_cross_validation/training12_test3/'test_data_trios_cluster/test_data_trios.txt \
    --sge 1

## cluster version
script_path='/home/rt2776/CN_Learn/scripts/'
data_path='/home/rt2776/CN_Learn/data5_spark27k_rare/'
qsub -t 1-1607 ${script_path}post_cn_learn_cluster.sh \
    ${data_path}'classify_cnvs_cross_validation/training12_test3/'CNV_list_with_predictions.txt \
    ${data_path}'consolidated_calls.bed' \
    ${data_path}'classify_cnvs_cross_validation/training12_test3/'test_data_trios_cluster0609/test_data_trios.txt

qsub -t 1-1130 ${script_path}post_cn_learn_cluster.sh \
    ${data_path}'classify_cnvs_cross_validation/training13_test2/'CNV_list_with_predictions.txt \
    ${data_path}'consolidated_calls.bed' \
    ${data_path}'classify_cnvs_cross_validation/training13_test2/'test_data_trios_cluster0609/test_data_trios.txt

qsub -t 1-1204 ${script_path}post_cn_learn_cluster.sh \
    ${data_path}'classify_cnvs_cross_validation/training23_test1/'CNV_list_with_predictions.txt \
    ${data_path}'consolidated_calls.bed' \
    ${data_path}'classify_cnvs_cross_validation/training23_test1/'test_data_trios_cluster0609/test_data_trios.txt

## modify and run the following script to check result file. If missing some files, then redo the upstream step for the 
## specific jobs
python /home/rt2776/cnv_analysis/scripts/cluster_results_checking.py

cat ${data_path}'classify_cnvs_cross_validation/training12_test3/'test_data_trios_cluster0609/test_data_trios* > \
    ${data_path}'classify_cnvs_cross_validation/training12_test3/'test_data_trio.txt

cat ${data_path}'classify_cnvs_cross_validation/training13_test2/'test_data_trios_cluster0609/test_data_trios* > \
    ${data_path}'classify_cnvs_cross_validation/training13_test2/'test_data_trio.txt

cat ${data_path}'classify_cnvs_cross_validation/training23_test1/'test_data_trios_cluster0609/test_data_trios* > \
    ${data_path}'classify_cnvs_cross_validation/training23_test1/'test_data_trio.txt

## Run CN-Learn Machine learning approach to get the prediction scores for the trios
source /home/rt2776/CN_Learn/config.params
LOWER_SIZE_LIMIT=1
UPPER_SIZE_LIMIT=50000000

GROUP_DATA_DIR=${data_path}'classify_cnvs_cross_validation/training12_test3/'
python -u ${SCRIPTS_DIR}cn_learn_trios.py  ${GROUP_DATA_DIR}  training_data.txt  test_data_trio.txt \
                                ${CLASSIFIER}  ${LOWER_SIZE_LIMIT}  ${UPPER_SIZE_LIMIT}  ${NUM_TREES} \
                                                            ${CALLER_COUNT}  ${CALLER_LIST}

GROUP_DATA_DIR=${data_path}'classify_cnvs_cross_validation/training13_test2/'
python -u ${SCRIPTS_DIR}cn_learn_trios.py  ${GROUP_DATA_DIR}  training_data.txt  test_data_trio.txt \
                                ${CLASSIFIER}  ${LOWER_SIZE_LIMIT}  ${UPPER_SIZE_LIMIT}  ${NUM_TREES} \
                                                            ${CALLER_COUNT}  ${CALLER_LIST}

GROUP_DATA_DIR=${data_path}'classify_cnvs_cross_validation/training23_test1/'
python -u ${SCRIPTS_DIR}cn_learn_trios.py  ${GROUP_DATA_DIR}  training_data.txt  test_data_trio.txt \
                                ${CLASSIFIER}  ${LOWER_SIZE_LIMIT}  ${UPPER_SIZE_LIMIT}  ${NUM_TREES} \
                                                            ${CALLER_COUNT}  ${CALLER_LIST}
## Combine Results
cat ${data_path}'classify_cnvs_cross_validation/training12_test3/'trio_CNV_list_with_predictions.csv \
    ${data_path}'classify_cnvs_cross_validation/training13_test2/'trio_CNV_list_with_predictions.csv \
    ${data_path}'classify_cnvs_cross_validation/training23_test1/'trio_CNV_list_with_predictions.csv \
    >${data_path}trio_CNV_list_with_predictions.csv

cat ${data_path}'classify_cnvs_cross_validation/training12_test3/'trio_CNV_list_with_predictions.txt \
    ${data_path}'classify_cnvs_cross_validation/training13_test2/'trio_CNV_list_with_predictions.txt \
    ${data_path}'classify_cnvs_cross_validation/training23_test1/'trio_CNV_list_with_predictions.txt \
    >${data_path}trio_CNV_list_with_predictions.txt

## don't forget to delete two redundant top rows
