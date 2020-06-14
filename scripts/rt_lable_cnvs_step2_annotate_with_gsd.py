##############################################################################
# Function: Label "True" or "False" by GSD dataset  
#           Take the union of inherited CNVs identified by XHMM, CANOES and CLAMMS as TRUE
#           Take the union of denovo CNVs identified by XHMM, CANOES and CLAMMS as FALSE
#           
# Input:    Prepare_training_testing_file.txt
# Output:   
# Author:   Renjie Tan
# Date:     05/20/2020
##############################################################################
import os
import pdb
import datetime

def fileToList(file_name):
    result_list = []
    fp = open(file_name)
    for row in fp:
        row = row.strip()
        row = row.replace("\"","")
        # row = re.split('\t| ',row)
        row = row.split()
        result_list.append(row)
    fp.close()
    return result_list

def output_to_file(results,file_name):
    path = os.path.dirname(file_name)
    if not os.path.exists(path) and path != '':
        os.mkdir(path)
    fp = open(file_name,'w')
    if type(results) == list:
        for result_line in results:
            result_line = [str(x) for x in result_line]
            fp.write("\t".join(result_line) + os.linesep)
    elif type(results) == dict:
        for key in results.keys():
            for result_line in results[str(key)]:
                result_line = [str(x) for x in result_line]
                fp.write("\t".join(result_line) + os.linesep)
    else:
        print('[ERROR]: Unsupport results type!')
        return
    print('[INFO]: File outputs to %s'%file_name)
    time_stamp = datetime.datetime.now()
    print(time_stamp.strftime('[%Y.%m.%d-%H:%M:%S]'))
    fp.close()

def fileToDict(file_name, start_row=0, key_col=0):
    result_dict = {}
    fp = open(file_name)
    num = 0
    for row in fp:
        num += 1
        if num < start_row:
            continue
        row = row.strip()
        row = row.split()
        row = [s.replace('"','') for s in row]

        if row[key_col] in result_dict:
            result_dict[row[key_col]].append(row)
        else:
            result_dict[row[key_col]] = [row]
    return result_dict

def reciprocalOverlap(chrom_data1, start_data1, stop_data1, type_data1, chrom_data2, start_data2, stop_data2, type_data2, thresthold):
    if chrom_data1 == chrom_data2:
        start_data1 = int(start_data1)
        stop_data1 = int(stop_data1)
        start_data2 = int(start_data2)
        stop_data2 = int(stop_data2)
        if type_data1 == type_data2:
            if start_data1 in range(start_data2, stop_data2+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data1 + 1
                if thresthold == 1:  #1bp overlap
                    if overlap > 0:
                        return True
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data2 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return True
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            elif start_data1 in range(start_data2, stop_data2+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data1 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return True
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data2 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return True
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            else:
                return False
        else:
            return False
    else:
        return False

def fetchValidation(validated_list, cnv_chr, cnv_start, cnv_stop, cnv_type):
    result = None 
    for reader in validated_list:
        try:
            region_chr = reader[0].replace("chr","")
            region_start = int(reader[1])
            region_stop = int(reader[2])
            region_type = reader[3]
            region_caller = reader[7]
            cnv_start = int(cnv_start)
            cnv_stop = int(cnv_stop)
            if reciprocalOverlap(cnv_chr, cnv_start, cnv_stop, cnv_type, region_chr, region_start, region_stop, region_type, 0.5) == True:
                result = region_caller+"_"+str(region_chr)+":"+str(region_start)+"-"+str(region_stop)+"_"+region_type
                break
        except:
            pass
    return result

gsd_path='/home/rt2776/CN_Learn/gsd_data/'
data_path='/home/rt2776/CN_Learn/data5_spark27k_rare/'

validated_true_file = gsd_path+'validated_canoes_xhmm_clamms_true_cnvs.txt'
validated_false_file = gsd_path+'validated_canoes_xhmm_clamms_false_cnvs.txt'
cn_predict_file = data_path+'prepare_training_testing_file.txt'
result_file = data_path+"final_preds_with_lable.txt"

cn_predict_list = fileToList(cn_predict_file)
validated_true_dict = fileToDict(validated_true_file, key_col=5)
validated_false_dict = fileToDict(validated_false_file, key_col=5)
result_list = []

num=0
for reader in cn_predict_list:
    num+=1
    label_result = None
    if num == 1:
        reader.extend(["LABEL_VAL","LABEL_DETAIL"])
    else:
        sampleID=reader[4]
        pre_chr = reader[0]
        pre_start = reader[1] 
        pre_stop = reader[2]
        pre_type = reader[3]
        label_value,label_detail = '-','-'

        if sampleID in validated_true_dict:
            label_result = fetchValidation(validated_true_dict[sampleID], pre_chr, pre_start, pre_stop, pre_type)
            if label_result != None:
                label_value = 1
                label_detail = label_result
                print(num,reader,label_value, label_detail)
        if sampleID in validated_false_dict:
            label_result = fetchValidation(validated_false_dict[sampleID], pre_chr, pre_start, pre_stop, pre_type) 
            if label_result != None:
                if label_value == 1:
                    label_value = 'Conflict'
                    label_detail = 'T:' + label_detail + '__F:' + label_result
                else:
                    label_value = 0
                    label_detail = label_result
                print(num,reader,label_value, label_detail)

        reader.extend([label_value, label_detail])
    
    result_list.append(reader)

output_to_file(result_list, file_name = result_file)
