########################################################################
# Script : post_cn_learn_cluster.sh                                      #
# Author : Renjie Tan                                                  #
# Date   : 6/8/2020                                                    #
########################################################################
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=4G
#$ -l mem_free=4G
#$ -pe smp 1
#$ -o '/share/terra/rt2776/SPARK/CN_Learn/oe/'
#$ -e '/share/terra/rt2776/SPARK/CN_Learn/oe/'

echo "Job started on `hostname` at `date`"
echo "Input:"$1
echo "Consolidated_calls:"$2
echo "Output:"$3

python /home/rt2776/CN_Learn/scripts/post_cn_learn.py prepare_trio_scores \
    --input $1 \
    --consolidated_calls $2 \
    --output $3 \
    --sge $SGE_TASK_ID

echo "Job ended on `hostname` at `date`"
