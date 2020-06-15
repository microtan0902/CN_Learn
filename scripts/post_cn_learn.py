##############################################################################
# Function: Post CN-Learn data processes 
# Input:   
# Output:  
# Author:   Renjie Tan
# Date:     06/02/2020
##############################################################################
from __future__ import division
import os
import pdb
import datetime
import function as df
import numpy as np
import argparse

def reciprocalOverlap(chrom_data1, start_data1, stop_data1, type_data1, chrom_data2, start_data2, stop_data2, type_data2, thresthold):
    if chrom_data1 == chrom_data2:
        start_data1 = int(start_data1)
        stop_data1 = int(stop_data1)
        start_data2 = int(start_data2)
        stop_data2 = int(stop_data2)
        if type_data1 == type_data2: 
            if start_data1 in range(start_data2, stop_data2+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data1 
                if thresthold == 1:  #1bp overlap
                    if overlap > 0:
                        return overlap
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return overlap
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data2
                if thresthold == 1:
                    if overlap > 0:
                        return overlap
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            elif start_data1 in range(start_data2, stop_data2+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data1 
                if thresthold == 1:
                    if overlap > 0:
                        return overlap
                else:    
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return overlap
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data2
                if thresthold == 1:
                    if overlap > 0:
                        return overlap
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return overlap
            else:
                return 0
        else:
            return 0
    else:
        return 0

def fetch_consolidated_calls(sampleID, cnv_chr, cnv_start, cnv_stop,cnv_type, consolidated_cnv_list):
    overlap_w_canoes,overlap_w_clamms,overlap_w_xhmm = 0,0,0
    for cons_cnv_reader in consolidated_cnv_list:
        cons_cnv_chr = cons_cnv_reader[0]
        cons_cnv_start = int(cons_cnv_reader[1])
        cons_cnv_stop = int(cons_cnv_reader[2])
        cons_cnv_type = cons_cnv_reader[3]
        cons_cnv_caller = cons_cnv_reader[5]

        overlap_region = reciprocalOverlap(cnv_chr, cnv_start, cnv_stop, cnv_type, cons_cnv_chr, cons_cnv_start, 
                cons_cnv_stop, cons_cnv_type, 1)
        if overlap_region != 0:
            # Note: the ratio is divied by the length of CNV prediced by CN-Learn, not each tools' prediction
            # so using "overlap_ratio = round(float(overlap_region/(cons_cnv_stop-cons_cnv_start)),2)" is problematic
            overlap_ratio = round(float(overlap_region/(cnv_stop-cnv_start)),2)
            if overlap_ratio == 1.0:
                overlap_ratio = 1

            if cons_cnv_caller == "XHMM":
                overlap_w_xhmm = overlap_ratio 
            elif cons_cnv_caller == "CANOES":
                overlap_w_canoes = overlap_ratio 
            elif cons_cnv_caller == "CLAMMS":
                overlap_w_clamms = overlap_ratio 
            else:
                pdb.set_trace()

    return [overlap_w_canoes,overlap_w_clamms,overlap_w_xhmm]

def get_targets_num(cnv_chr, cnv_start, cnv_stop, target_probe_list):
    num = -1 #list start from 0
    cnv_target_num_list = []
    target_start_num, target_stop_num = None, None
    for target_reader in target_probe_list:
        num += 1
        target_chr = target_reader[0]
        target_start = int(target_reader[1])
        target_stop = int(target_reader[2])
        if cnv_chr == target_chr:
            if target_start in range(cnv_start, cnv_stop+1) and target_stop in range(cnv_start, cnv_stop+1):
                cnv_target_num_list.append(num)
        else:
            if cnv_target_num_list != []:
                break
    return cnv_target_num_list

def get_targets_num_w_flanking(cnv_target_num_list, adjust_target_num):
    left_flank_target_num_list,right_flank_target_num_list=[],[]
    try:
        left_flank_target_num_list = list(range(cnv_target_num_list[0]-adjust_target_num, cnv_target_num_list[0]))
        right_flank_target_num_list = list(range(cnv_target_num_list[-1]+1, cnv_target_num_list[-1]+adjust_target_num+1))
        #TODO: what should return if the CNV start/stop target is the first/last target in a chromosome?
    except:
        pdb.set_trace()
    return [left_flank_target_num_list,right_flank_target_num_list]

def get_targets_rd(bp_coverage_file, target_num_list, target_probe_list):
    rd_list = []
    for num in target_num_list:
        target_reader = target_probe_list[num]
        target_chr = target_reader[0]
        target_start = int(target_reader[1])
        target_stop = int(target_reader[2])
        #target_rd = df.fetchReadDepth(cram_file, target_chr, target_start, target_stop, MAQ=0)
        target_rd = df.fetchReadDepthFromTabularFile(bp_coverage_file, target_chr, target_start, target_stop)
        rd_list.append(target_rd)
    rd = np.mean(rd_list)
    return rd

def fetch_rd_ratio(bp_coverage_file, cnv_chr, cnv_start, cnv_stop, num_targets, target_probe_list):
    half_num_target = max(int(num_targets/2),1) #TODO: need to check the results with CN-learn
    cnv_target_num_list = get_targets_num(cnv_chr, cnv_start, cnv_stop, target_probe_list)
    left_flank_target_num_list,right_flank_target_num_list = get_targets_num_w_flanking(cnv_target_num_list, half_num_target)
    cnv_target_rd = get_targets_rd(bp_coverage_file, cnv_target_num_list, target_probe_list)
    left_flank_rd = get_targets_rd(bp_coverage_file, left_flank_target_num_list, target_probe_list)
    right_flank_rd = get_targets_rd(bp_coverage_file, right_flank_target_num_list, target_probe_list)
    rd_ratio = round(cnv_target_rd/((left_flank_rd+right_flank_rd)/2),2)
    #print(cnv_target_num_list, left_flank_target_num_list, right_flank_target_num_list)
    #print(cnv_target_rd, left_flank_rd, right_flank_rd, rd_ratio)
    return rd_ratio

def prepare_trio_CNLearn_score(args):
    cn_predict_file = args.input[0]
    consolidated_calls_file = args.consolidated_calls[0]
    
    if args.sge == None:
        sge_task_id = None
    else:
        sge_task_id = int(args.sge[0])

    output_file = args.output[0]
    directory = os.path.dirname(output_file)
    if not os.path.exists(directory):
        os.makedirs(directory)
    file_name, file_extension = os.path.splitext(output_file)
    if sge_task_id != None:
        output_file = file_name+'_'+str(sge_task_id).zfill(4)+file_extension
    if os.path.isfile(output_file):
        print("%s exists. If you want to calculate it again, please delete it at first."%output_file)
        return None

    cn_predict_list = df.fileToList(cn_predict_file)
    consolidated_dict = df.fileToDict(consolidated_calls_file, start_row=0, key_col=4)
    target_probe_list = df.fileToList(target_probe_file)
    result_list = []
    if sge_task_id == 1:
        result_list.append(header)
    ped_indiv_dict, ped_family_dict = df.parse_pedigreeFile(pedigree_file)

    num=0
    total_num = len(cn_predict_list[1:])
    for reader in cn_predict_list[1:]:
        num+=1
        if sge_task_id == None or num in range((sge_task_id-1)*100+1, min(sge_task_id*100+1,len(cn_predict_list)+1)):
            cnv_chr = reader[0]
            cnv_start = int(reader[1])
            cnv_stop = int(reader[2])
            cnv_type = reader[3]
            sampleID = reader[4]
            
            GC = round(float(reader[7]),2)
            cnv_size = int(reader[8])
            MAP = round(float(reader[9]),2)
            num_targets = int(reader[10])
            sizeID = reader[11]
            cn_learn_score = reader[15]

            if sampleID in ped_indiv_dict:
                FamID = ped_indiv_dict[sampleID][0]
                PaternalID = ped_indiv_dict[sampleID][2]
                MaternalID = ped_indiv_dict[sampleID][3]
                Role = ped_indiv_dict[sampleID][6]
                Phenotype = ped_indiv_dict[sampleID][5]
                Sex = ped_indiv_dict[sampleID][4]
                Ancestry = ped_indiv_dict[sampleID][7]
                Family_Member = len(ped_family_dict[FamID])
                OffspringID = ped_family_dict[FamID][0][1]
            else:
                FamID, Role, Phenotype, Sex, Family_Member, PaternalID, MaternalID, Ancestry = '-','-','-','-','-','-','-','-'

            if Role == 'Offspring' and PaternalID in consolidated_dict and MaternalID in consolidated_dict and num_targets >=3:
                # Note: the boundaries have some problems when the num_targets==1.
                bp_coverage_file_offspring = bp_coverage_path + sampleID + '.bpcov.bed.gz'
                bp_coverage_file_paternal = bp_coverage_path + PaternalID + '.bpcov.bed.gz'
                bp_coverage_file_maternal = bp_coverage_path + MaternalID + '.bpcov.bed.gz'
                
                ## RD ratio 
                ofs_rd_prob = fetch_rd_ratio(bp_coverage_file_offspring,cnv_chr,cnv_start,cnv_stop,num_targets,target_probe_list)
                pat_rd_prob = fetch_rd_ratio(bp_coverage_file_paternal,cnv_chr,cnv_start,cnv_stop,num_targets,target_probe_list) 
                mat_rd_prob = fetch_rd_ratio(bp_coverage_file_maternal,cnv_chr,cnv_start,cnv_stop,num_targets,target_probe_list)
                
                ## Overlap condition
                ### -Offspring
                ofs_overlap_w_canoes,ofs_overlap_w_clamms,ofs_overlap_w_xhmm = fetch_consolidated_calls(sampleID, cnv_chr, cnv_start, cnv_stop,
                        cnv_type, consolidated_dict[sampleID])
                ofs_num_overlaps = 0
                for ofs_overlap_reader in [ofs_overlap_w_canoes,ofs_overlap_w_clamms,ofs_overlap_w_xhmm]:
                    if ofs_overlap_reader != 0:
                        ofs_num_overlaps += 1
                ### -Father
                pat_overlap_w_canoes,pat_overlap_w_clamms,pat_overlap_w_xhmm = fetch_consolidated_calls(PaternalID, cnv_chr, cnv_start, cnv_stop,
                        cnv_type, consolidated_dict[PaternalID])
                pat_num_overlaps = 0
                for pat_overlap_reader in [pat_overlap_w_canoes,pat_overlap_w_clamms,pat_overlap_w_xhmm]:
                    if pat_overlap_reader != 0:
                        pat_num_overlaps += 1
                ### -Mother
                mat_overlap_w_canoes,mat_overlap_w_clamms,mat_overlap_w_xhmm = fetch_consolidated_calls(MaternalID, cnv_chr, cnv_start, cnv_stop,
                        cnv_type, consolidated_dict[MaternalID])
                mat_num_overlaps = 0
                for mat_overlap_reader in [mat_overlap_w_canoes,mat_overlap_w_clamms,mat_overlap_w_xhmm]:
                    if mat_overlap_reader != 0:
                        mat_num_overlaps += 1
                
                ## output 
                ofs_result_row = [cnv_chr, cnv_start, cnv_stop, cnv_type, sampleID, ofs_overlap_w_canoes, ofs_overlap_w_clamms,
                    ofs_overlap_w_xhmm, ofs_num_overlaps, ofs_rd_prob, GC, cnv_size, MAP, num_targets, size_label_dict[sizeID], cn_learn_score, 'Offspring', FamID]
                pat_result_row = [cnv_chr, cnv_start, cnv_stop, cnv_type, PaternalID, pat_overlap_w_canoes, pat_overlap_w_clamms,
                    pat_overlap_w_xhmm, pat_num_overlaps, pat_rd_prob, GC, cnv_size, MAP, num_targets, size_label_dict[sizeID], '-', 'Father', FamID]
                mat_result_row = [cnv_chr, cnv_start, cnv_stop, cnv_type, MaternalID, mat_overlap_w_canoes, mat_overlap_w_clamms,
                    mat_overlap_w_xhmm, mat_num_overlaps, mat_rd_prob, GC, cnv_size, MAP, num_targets, size_label_dict[sizeID], '-', 'Mather', FamID]
                result_list.append(ofs_result_row)
                result_list.append(pat_result_row)
                result_list.append(mat_result_row)
                print('[%d|%d]%s----------------------------------------------------------'%(num,total_num,sampleID))
                print(ofs_result_row)
                print(pat_result_row)
                print(mat_result_row)

    df.output_to_file(result_list, file_name = output_file)


if __name__ == '__main__':
    size_label_dict = {"1":"A)<1KB","2":"B)1KB-5KB","3":"C)5KB-10KB","4":"D)10KB-25KB","5":"E)25KB-50KB",
            "6":"F)50KB-75KB","7":"G)75KB-100KB","8":"H)100KB-250KB","9":"I)250KB-500KB",
            "10":"J)500KB-1MB","11":"K)1MB-5MB","12":"L)>5MB"}
    header = ['CHR','PRED_START','PRED_END','TYPE', 'SAMPLE','CANOES','CLAMMS','XHMM','NUM_OVERLAPS','RD_PROP','GC',
                'PRED_SIZE','MAP','NUM_TARGETS','SIZE_LABEL','CN_LEARN_SCORE','ROLE','FAMID']
    ##input
    pedigree_file = '/home/rt2776/source/pedigree/SPARK30K_20190729_extention_ancestry.ped'
    target_probe_file = '/home/rt2776/CN_Learn/source/xgen_plus_spikein.b38.bed'
    cram_file = '/home/rt2776/CN_Learn/source/spark_cram_list_20190729_merged.txt'
    bp_coverage_path = '/home/rt2776/SPARK/CN_Learn/bp_coverage_dir/'
    #data_path = '/home/rt2776/CN_Learn/data5_spark27k_rare/'
    #consolidated_calls_file = data_path+'consolidated_calls.bed'
    #cn_predict_file = data_path+'classify_cnvs_cross_validation/training12_test3/CNV_list_with_predictions.txt'
    
    ##output
    #result_file = data_path+"classify_cnvs_cross_validation/training12_test3/test_data_parents.txt"

    parser = argparse.ArgumentParser(prog="Post CN-Learn process", description="Post CN-Learn process.")
    subparsers = parser.add_subparsers(help='Command to be run.')
    post_parser = subparsers.add_parser('prepare_trio_scores', help="prepare for getting trio CN-Learn scores.")
    post_parser.add_argument('--input',action='store',required=True, nargs=1)
    post_parser.add_argument('--consolidated_calls',required=True, nargs=1)
    post_parser.add_argument('--output',action='store',required=True, nargs=1)
    post_parser.add_argument('--sge',action='store',required=False, default=None, nargs=1)
    post_parser.set_defaults(func=prepare_trio_CNLearn_score)
    
    args = parser.parse_args()
    args.func(args)
