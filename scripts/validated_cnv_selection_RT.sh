#!/bin/bash
############################################################################
# Script : inherited_cnvs_for_CNLearn.sh                                   #
# Author : Renjie Tan                                                      #
# Date   : 4/27/2020                                                       #
#                                                                          #
# Here just using CANOES results as GSD. XHMM results haven't employed.    #
# Using inherited CNVs as positive and de novo CNVs as negative in the     #
# training data                                                            #
#
#
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
# STEP 1: XHMM
##############################################################################################
# unique and resort !!!
: sort u
data_path='/home/rt2776/CN_Learn/data/'
python /home/rt2776/cnv_analysis/scripts/5_annotation.py batch \
    --input ${data_path}xhmm_calls.txt \
    --output ${data_path}gsd_data/xhmm_calls_batch.txt \
    --project spark

# Annotate pedigree info and divide samples into parents and offsprings
python /home/rt2776/cnv_analysis/scripts/5_annotation.py pedigree \
    --input ${data_path}gsd_data/xhmm_calls_batch.txt \
    --output ${data_path}gsd_data/xhmm_calls_ped.cnv \
    --ped /home/rt2776/source/pedigree/SPARK30K_20190729_extention_ancestry.ped

# Divide samples into Parent and Offspring
cat  ${data_path}gsd_data/xhmm_calls_ped.cnv \
    |awk '{if ($22 == "Role" || $22 == "Parent") print $0}' \
    > ${data_path}gsd_data/xhmm_calls_parent.cnv 
    
cat  ${data_path}gsd_data/xhmm_calls_ped.cnv \
    |awk '{if ($22 == "Role" || $22 == "Offspring") print $0}' \
    > ${data_path}gsd_data/xhmm_calls_offspring.cnv

# some of samples don't have pedigree info
cat xhmm_calls_ped.cnv|cut -f 22|sort|uniq

# Genotyping for CNVs in offsprings
cat ${data_path}gsd_data/xhmm_calls_offspring.cnv | \
    awk '{if(($23=="ParentalID"||$23 != "0") && ($24=="MaternallID"||$24!="0")) {print $0}}' \
    > ${data_path}gsd_data/xhmm_calls_offspring_trio.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py genotype \
    --input ${data_path}gsd_data/xhmm_calls_offspring_trio.cnv \
    --output ${data_path}gsd_data/xhmm_calls_offspring_trio_gt.cnv \
    --vcf /home/rt2776/SPARK/2_xhmm/3_xhmm_cnv_results/byChromosomesGender \
    --sge 2652

qsub -t 1-2651 /home/rt2776/cnv_analysis/scripts/cluster_genotyping_scores.sh \
    ${data_path}gsd_data/xhmm_calls_offspring_trio.cnv \
    ${data_path}gsd_data/xhmm_genotyping/xhmm_calls_offspring_trio_gt.cnv \
    /home/rt2776/SPARK/2_xhmm/3_xhmm_cnv_results/byChromosomesGender

cat ${data_path}gsd_data/xhmm_genotyping/* > ${data_path}gsd_data/xhmm_calls_offspring_trio_gt.cnv
# above here 04/29/2020

# Filter CNVs with number of target region < 3
data_path='/home/rt2776/CN_Learn/gsd_data/xhmm/'
input_file=${data_path}xhmm_calls_offspring_trio_gt.cnv
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="NUM_TARG") {print i}}}'

cat ${data_path}xhmm_calls_offspring_trio_gt.cnv |awk '{if($8=="NUM_TARG" || $8 >= 3) {print $0}}' \
    > ${data_path}xhmm_calls_offspring_3targets.cnv

# Select Mendelian erros (de novo CNVs) as FALSE label. Offspring.SQ>=10 & parents.NQ>=60
# 2020.6.13 updated.
input_file=${data_path}xhmm_calls_offspring_3targets.cnv
head -1 ${input_file} | awk '{for (i=1;i<=NF;i++) {if ($i=="Offspring_SQ") {print i}}}'
head -1 ${input_file} | awk '{for (i=1;i<=NF;i++) {if ($i=="Father_NQ") {print i}}}'
head -1 ${input_file} | awk '{for (i=1;i<=NF;i++) {if ($i=="Mother_NQ") {print i}}}'

Offspring_SQ=10
Parent_NQ=60
cat ${input_file} | grep -v "None" \
    | awk -F '\t' '{if(($26=="Offspring_SQ"||$26>='${Offspring_SQ}') && ($29=="Father_NQ"||$29>='${Parent_NQ}') \
    &&($31=="Mother_NQ"||$31>='${Parent_NQ}')) print $0 }' \
    > ${data_path}xhmm_calls_offspring_denovo.cnv

cat ${data_path}xhmm_calls_offspring_denovo.cnv|cut -f 1,2,3,26,28,29,30,31|less

# Select inherited CNVs as TRUE label. Offspring.SQ>=60 & OneParent.SQ>=60 & OtherParent.NQ>=60
head -1 ${input_file} | awk '{for (i=1;i<=NF;i++) {if ($i=="CNV") {print i}}}'
head -1 ${input_file} | awk '{for (i=1;i<=NF;i++) {if ($i=="Offspring_SQ") {print i}}}'
head -1 ${input_file} | awk '{for (i=1;i<=NF;i++) {if ($i=="Father_SQ") {print i}}}'
head -1 ${input_file} | awk '{for (i=1;i<=NF;i++) {if ($i=="Mother_SQ") {print i}}}'

cat ${input_file} |grep -v "None"| awk -F '\t' \
    '{if( ($26=="Offspring_SQ"||$26>=60) && \
    (($28=="Father_SQ"||$28>=60) && ($31=="Mother_NQ"||$31>=60) || \
    ($30=="Mother_SQ"||$30>=60) && ($29=="Father_NQ"||$29>=60)) \
    ) print $0}'>${data_path}xhmm_calls_offspring_inherited.cnv

cat ${data_path}xhmm_calls_offspring_inherited.cnv|cut -f 1,2,3,26,28,29,30,31|less

# annotate SD
python /home/rt2776/cnv_analysis/scripts/5_annotation.py sd \
    --input ${data_path}xhmm_calls_offspring_denovo.cnv \
    --output ${data_path}xhmm_calls_offspring_denovo_sd.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py sd \
    --input ${data_path}xhmm_calls_offspring_inherited.cnv \
    --output ${data_path}xhmm_calls_offspring_inherited_sd.cnv

# annotate Mappability
python /home/rt2776/cnv_analysis/scripts/5_annotation.py mappability \
    --input ${data_path}xhmm_calls_offspring_denovo_sd.cnv \
    --output ${data_path}xhmm_calls_offspring_denovo_sd_mapp.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py mappability \
    --input ${data_path}xhmm_calls_offspring_inherited_sd.cnv \
    --output ${data_path}xhmm_calls_offspring_inherited_sd_mapp.cnv

# annotate GC
python /home/rt2776/cnv_analysis/scripts/5_annotation.py gc \
    --input ${data_path}xhmm_calls_offspring_denovo_sd_mapp.cnv \
    --output ${data_path}xhmm_calls_offspring_denovo_sd_mapp_gc.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py gc \
    --input ${data_path}xhmm_calls_offspring_inherited_sd_mapp.cnv \
    --output ${data_path}xhmm_calls_offspring_inherited_sd_mapp_gc.cnv

# Filtering denovo
input_file=${data_path}xhmm_calls_offspring_denovo_sd_mapp_gc.cnv
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="SD_region") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="Mappability") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="GC") {print i}}}'

cat ${input_file}|awk -F '\t' '{if(($32=="SD_region"||$32=="-") && ($33=="Mappability"||$33>=0.75) && ($34=="GC"||($34>=0.3&&$34<=0.7))) print $0}' \
    > ${data_path}xhmm_calls_offspring_denovo_final.cnv

cat ${input_file}|awk -F '\t' '{if(($32=="SD_region"||$32!="-") || ($33=="Mappability"||$33<0.75) || ($34=="GC"||$34<0.3||$34>0.7)) print $0}' \
    > ${data_path}xhmm_calls_offspring_denovo_final_Removed.cnv

cat ${data_path}xhmm_calls_offspring_denovo_final.cnv|cut -f 3,2,1 > ${data_path}validation_xhmm_denovo.txt
cat ${data_path}validation_xhmm_denovo.txt|awk -F "\t|-|:" '{print $3"\t"$4"\t"$5"\t"$2"\t"$5-$4"\t"$1"\tFalse\tXHMM"}' \
    >${data_path}validated_xhmm_false_cnvs.txt

# Filtering inherited
input_file=${data_path}xhmm_calls_offspring_inherited_sd_mapp_gc.cnv
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="SD_region") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="Mappability") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="GC") {print i}}}'

cat ${input_file}|awk -F '\t' '{if(($32=="SD_region"||$32=="-") && ($33=="Mappability"||$33>=0.75) && ($34=="GC"||($34>=30&&$34<=70))) print $0}' \
    > ${data_path}xhmm_calls_offspring_inherited_final.cnv

cat ${input_file}|awk -F '\t' '{if(($32=="SD_region"||$32!="-") || ($33=="Mappability"||$33<0.75) || ($34=="GC"||$34<30||$34>70)) print $0}' \
    > ${data_path}xhmm_calls_offspring_inherited_final_Removed.cnv

cat ${data_path}xhmm_calls_offspring_inherited_final.cnv|cut -f 3,2,1 > ${data_path}validation_xhmm_inherited.txt
cat ${data_path}validation_xhmm_inherited.txt|awk -F "\t|-|:" '{print $3"\t"$4"\t"$5"\t"$2"\t"$5-$4"\t"$1"\tTrue\tXHMM"}' \
    >${data_path}validated_xhmm_true_cnvs.txt


##############################################################################################
# STEP 2.1: CANOES CNVs 
##############################################################################################
data_path='/home/rt2776/CN_Learn/gsd_data/canoes/'
cd ${data_path}
## Add reference
cat ${data_path}canoes_calls_offspring_trio_gt_target3.cnv | awk '{print $0"\thg38"}' \
    >${data_path}canoes_calls_offspring_trio_gt_target3_hg38.cnv

## Got inherited CNVs as TP and de novo CNVs as TN
#replace all " in file by vim :%s/"//g
cat ${data_path}canoes_calls_offspring_trio_gt_target3_hg38.cnv | \
    awk -F '\t' '{if(($11=="Q_SOME"||$11>=70) && ($19!="0") && ($20!="0")) print $0}' \
    > ${data_path}canoes_calls_offspring_trio_gt_target3_inherited.tmp

cat ${data_path}canoes_calls_offspring_trio_gt_target3_inherited.tmp | awk -F '\t' \
    '{if((($3 == "CNV" || $3 == "DEL") && ($24 == "Paternal_SQDel" || $24 >= 70) && ($28 == "Maternal_NQDel" || $28 >= 70) ||
         ($3 == "CNV" || $3 == "DEL") && ($23 == "Paternal_NQDel" || $23 >= 70) && ($29 == "Maternal_SQDel" || $29 >= 70) ||
         ($3 == "CNV" || $3 == "DUP") && ($25 == "Paternal_NQDup" || $25 >= 70) && ($31 == "Maternal_SQDup" || $31 >= 70) ||
         ($3 == "CNV" || $3 == "DUP") && ($26 == "Paternal_SQDup" || $26 >= 70) && ($30 == "Maternal_NQDup" || $30 >= 70)) &&
         $23!="NA" && $28!="NA") print $0}'\
    >${data_path}canoes_calls_offspring_trio_gt_target3_inherited.cnv 

## Select Mendelian erros (de novo CNVs) as FALSE label. Offspring.SQ>=10 & parents.NQ>=60
## 2020.6.13 updated.
SQ_threshold=10
cat ${data_path}canoes_calls_offspring_trio_gt_target3_hg38.cnv | \
    awk -F '\t' '{if(($11=="Q_SOME"||$11>='${SQ_threshold}') && ($19!="0") && ($20!="0")) print $0}' \
    > ${data_path}canoes_calls_offspring_trio_gt_target3_denovo.tmp

NQ_threshold=70
cat ${data_path}canoes_calls_offspring_trio_gt_target3_denovo.tmp | awk -F '\t' \
    '{if((($3 == "CNV" || $3 == "DEL") && ($23 == "Paternal_NQDel" || $23 >= '${NQ_threshold}') && ($28 == "Maternal_NQDel" || $28 >= '${NQ_threshold}') ||
        ($3 == "CNV" || $3 == "DUP") && ($25 == "Paternal_NQDup" || $25 >= '${NQ_threshold}') && ($30 == "Maternal_NQDup" || $30 >= '${NQ_threshold}')) &&
         $23!="NA" && $28!="NA") print $0}' >${data_path}canoes_calls_offspring_trio_gt_target3_denovo.cnv

## annotate SD
python /home/rt2776/cnv_analysis/scripts/5_annotation.py sd \
    --input ${data_path}canoes_calls_offspring_trio_gt_target3_inherited.cnv \
    --output ${data_path}canoes_calls_offspring_trio_w_gt_t3_inherited_sd.cnv

cat ${data_path}canoes_calls_offspring_trio_w_gt_t3_inherited_sd.cnv|awk -F '\t' '{if($33=="SD_region" || $33=="-") print $0}' \
        > ${data_path}canoes_calls_offspring_trio_inherited_sd_filtered.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py sd \
    --input ${data_path}canoes_calls_offspring_trio_gt_target3_denovo.cnv \
    --output ${data_path}canoes_calls_offspring_trio_gt_target3_denovo_sd.cnv
cut -f33 ${data_path}canoes_calls_offspring_trio_gt_target3_denovo_sd.cnv|grep -w '-'|wc

## annotate mappability
python /home/rt2776/cnv_analysis/scripts/5_annotation.py mappability \
    --input ${data_path}canoes_calls_offspring_trio_inherited_sd_filtered.cnv \
    --output ${data_path}canoes_calls_offspring_trio_inherited_sd_mappability.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py mappability \
    --input ${data_path}canoes_calls_offspring_trio_gt_target3_denovo_sd.cnv \
    --output ${data_path}canoes_calls_offspring_trio_denovo_sd_mappability.cnv

## annotate GC
python /home/rt2776/cnv_analysis/scripts/5_annotation.py gc \
    --input ${data_path}canoes_calls_offspring_trio_inherited_sd_mappability.cnv \
    --output ${data_path}canoes_calls_offspring_trio_inherited_sd_mappability_gc.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py gc \
    --input ${data_path}canoes_calls_offspring_trio_denovo_sd_mappability.cnv \
    --output ${data_path}canoes_calls_offspring_trio_denovo_sd_mappability_gc.cnv

## Filter inherited CNVs
input_file=${data_path}canoes_calls_offspring_trio_inherited_sd_mappability_gc.cnv
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="Mappability") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="GC") {print i}}}'

cat ${input_file}|awk -F '\t' '{if(($34=="Mappability"||$34>=0.75) && ($35=="GC"||($35>30&&$35<=70))) print $0}' \
    > ${data_path}canoes_calls_offspring_trio_inherited_final.cnv 

cat ${data_path}canoes_calls_offspring_trio_inherited_final.cnv|cut -f 4,3,2 > ${data_path}validation_canoes_inherited.txt
cat ${data_path}validation_canoes_inherited.txt|awk -F "\t|-|:" '{print $3"\t"$4"\t"$5"\t"$2"\t"$5-$4+1"\t"$1"\tTrue\tCANOES"}' \
    >${data_path}validated_canoes_true_cnvs.txt

## Filter de novo CNVs
input_file=${data_path}canoes_calls_offspring_trio_denovo_sd_mappability_gc.cnv
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="SD_region") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="Mappability") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="GC") {print i}}}'

cat ${input_file}|awk -F '\t' '{if(($33=="SD_region"||wc$33=="-") && ($34=="Mappability"||$34>=0.75) && ($35=="GC"||($35>=0.3&&$35<=0.7))) print $0}' \
    > ${data_path}canoes_calls_offspring_trio_denovo_final.cnv 

cat ${input_file}|awk -F '\t' '{if(($33=="SD_region"||$33!="-") || ($34=="Mappability"||$34<0.75) || ($35=="GC"||$35<0.3||$35>0.7)) print $0}' \
    > ${data_path}canoes_calls_offspring_trio_denovo_final_filtered.cnv 

cat ${data_path}canoes_calls_offspring_trio_denovo_final.cnv|cut -f 4,3,2 > ${data_path}validation_canoes_denovo.txt
cat ${data_path}validation_canoes_denovo.txt|awk -F "\t|-|:" '{print $3"\t"$4"\t"$5"\t"$2"\t"$5-$4+1"\t"$1"\tFalse\tCANOES"}' \
    >${data_path}validated_canoes_false_cnvs.txt


##############################################################################################
# STEP 2.3: CLAMMS CNVs 
##############################################################################################
data_path='/home/rt2776/CN_Learn/gsd_data/clamms/'
cd ${data_path}
python /home/rt2776/cnv_analysis/scripts/5_annotation.py pedigree \
    --input ${data_path}clamms_calls.txt \
    --output ${data_path}clamms_calls_ped.txt \
    --ped /home/rt2776/source/pedigree/SPARK30K_20190729_extention_ancestry.ped 

## Filter low quality calls. "Any call with Q_EXACT < 0 is of questionable quality."
cat ${data_path}clamms_calls_ped.txt |awk '{if($10=="Q_EXACT"||$10>0) print $0}' >${data_path}clamms_calls_ped_EQabove0.txt

## Filter CNVs with number of target region < 3
cat ${data_path}clamms_calls_ped_EQabove0.txt |awk '{if($8=="num_target"||$8>3) print $0}'>${data_path}clamms_calls_ped_EQabove0_tar3.txt

## Divide samples into Parent and Offspring
cat ${data_path}clamms_calls_ped_EQabove0_tar3.txt |awk '{if($23=="Role"||$23=="Offspring") print $0}'>${data_path}clamms_calls_offspring.cnv
cat ${data_path}clamms_calls_ped_EQabove0_tar3.txt |awk '{if($23=="Role"||$23=="Parent") print $0}'>${data_path}clamms_calls_parent.cnv

## Stratify de novo and inherited CNVs
## output:/home/rt2776/CN_Learn/gsd_data/clamms/clamms_offspring_inheritance_w_SQ.cnv
python /home/rt2776/cnv_analysis/scripts/clamms_denovo_inherited_cnv.py

## annotate SD, Mappability, GC
python /home/rt2776/cnv_analysis/scripts/5_annotation.py sd \
    --input ${data_path}clamms_offspring_inheritance_w_SQ.cnv \
    --output ${data_path}clamms_offspring_inheritance_sd.cnv

python /home/rt2776/cnv_analysis/scripts/5_annotation.py mappability \
    --input ${data_path}clamms_offspring_inheritance_sd_dict.cnv \
    --output ${data_path}clamms_offspring_inheritance_sd_mapp.cnv

## GC method1
python /home/rt2776/cnv_analysis/scripts/5_annotation.py gc \
    --input ${data_path}clamms_offspring_inheritance_sd_mapp.cnv \
    --output ${data_path}clamms_offspring_inheritance_sd_mapp_gc.cnv

## GC calulation method2(much faster)
cnv_gc_region='/home/rt2776/CN_Learn/gsd_data/clamms/clamms_offspring_inheritance_gc.bed'
REF_GENOME='/share/data/RGC_b38/genome.hg38rg.fa'
data_path='/home/rt2776/CN_Learn/gsd_data/clamms/'
bedtools nuc -fi ${REF_GENOME} -bed ${cnv_gc_region} >${data_path}clamms_bedtools_gc.txt
paste ${data_path}clamms_offspring_inheritance_sd_mapp_gc.cnv ....

## Filter by SD, Mappability and GC content
input_file=${data_path}clamms_offspring_inheritance_sd_mapp_gc.cnv
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="SD_region") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="Mappability") {print i}}}'
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="GC") {print i}}}'

cat ${input_file} |awk -F '\t' '{if($7=="SD_region" || $7=="-") print $0}' |\
    awk '{if(($8=="Mappability"||$8>=0.75) && ($9=="GC"||($9>0.30 && $9<=0.70))) print$0}' \
        > ${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv

cut -f9 ${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv|sort|tail
cut -f9 ${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv|sort|head

cut -f8 ${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv|sort|tail
cut -f8 ${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv|sort|head

cut -f7 ${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv|sort|tail
cut -f7 ${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv|sort|head

## Get de novo and inherited CNVs
input_file=${data_path}clamms_offspring_inheritance_sd_mapp_gc_filtered.cnv
head -1 $input_file | awk '{for (i=1;i<=NF;i++) {if ($i=="inheritance") {print i}}}'

cat $input_file |awk '{if($5=="inheritance"||$5=="Denovo") print $0}'>${data_path}clamms_offspring_denovo.cnv
cat $input_file |awk '{if($5=="inheritance"||$5=="Inherited") print $0}'>${data_path}clamms_offspring_inherited.cnv

## check the distribution of SQ and select suitable part of CNVs as True or False
## 2020.06.14
github cnv_toolkit/plot/hist_clamms_sq_distirbution.R

q_threshold=80
cat ${data_path}clamms_offspring_denovo.cnv|awk '{if($6=="Q_SOME"||$6<'${q_threshold}') print $0}' > ${data_path}clamms_offspring_denovo_sq80.cnv

cat ${data_path}clamms_offspring_denovo_sq80.cnv|awk -F "\t|-|:" '{print $3"\t"$4"\t"$5"\t"$2"\t"$5-$4"\t"$1"\tFalse\tCLAMMS"}' \
    >${data_path}validated_clamms_false_cnvs.txt

cat ${data_path}clamms_offspring_inherited.cnv | awk -F "\t|-|:" '{print $3"\t"$4"\t"$5"\t"$2"\t"$5-$4"\t"$1"\tTrue\tCLAMMS"}' \
    >${data_path}validated_clamms_true_cnvs.txt

##############################################################################################
# STEP 3: Combine results
##############################################################################################
data_path='/home/rt2776/CN_Learn/gsd_data/'
cat ${data_path}xhmm/validated_xhmm_true_cnvs.txt ${data_path}canoes/validated_canoes_true_cnvs.txt \
    ${data_path}clamms/validated_clamms_true_cnvs.txt > ${data_path}validated_xhmm_canoes_clamms_true_cnvs.txt

cat ${data_path}xhmm/validated_xhmm_false_cnvs.txt ${data_path}canoes/validated_canoes_false_cnvs.txt \
    ${data_path}clamms/validated_clamms_false_cnvs.txt > ${data_path}validated_xhmm_canoes_clamms_false_cnvs.txt
