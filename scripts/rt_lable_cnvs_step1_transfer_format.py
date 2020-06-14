##############################################################################
# Function: Transfer CN-Learn consensus CNV calls with GC, mappability and 
#           target info to training and testing files format 
# Input:    final_preds_GC_Map_Targ.txt
# Output:   Prepare_training_testing_file.txt
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
    if not os.path.exists(path) and path!="":
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

data_path = '/home/rt2776/CN_Learn/data5_spark27k_rare/'
cn_predict_file = data_path+'final_preds_GC_Map_Targ.txt'
cn_predict_list = fileToList(cn_predict_file)
num=0
result_list = []
header = ['CHR','PRED_START','PRED_END','TYPE', 'SAMPLE','CANOES','CLAMMS','XHMM','NUM_OVERLAPS','RD_PROP','GC',
        'PRED_SIZE','MAP','NUM_TARGETS','SIZE_LABEL']
result_file = data_path+"prepare_training_testing_file.txt"
result_list.append(header)

for reader in cn_predict_list:
    num+=1
    GC = round(float(reader[10]),2)
    MAP = round(float(reader[12]),2)
    cnv_size = int(reader[11])
    if cnv_size < 1000:
        size_lable = "A)<1KB"
    elif cnv_size in range(1000, 5000):
        size_lable = "B)1KB-5KB"
    elif cnv_size in range(5000, 10*1000):
        size_lable = "C)5KB-10KB"
    elif cnv_size in range(10*1000, 25*1000):
        size_lable = "D)10KB-25KB"
    elif cnv_size in range(25*1000, 50*1000):
        size_lable = "E)25KB-50KB"
    elif cnv_size in range(50*1000, 75*1000):
        size_lable = "F)50KB-75KB"
    elif cnv_size in range(75*1000, 100*1000):
        size_lable = "G)75KB-100KB"
    elif cnv_size in range(100*1000, 250*1000):
        size_lable = "H)100KB-250KB"
    elif cnv_size in range(250*1000, 500*1000):
        size_lable = "I)250KB-500KB"
    elif cnv_size in range(500*1000, 1000*1000):
        size_lable = "J)500KB-1MB"
    elif cnv_size in range(1000*1000, 5000*1000+1):
        size_lable = "K)1MB-5MB"
    elif cnv_size > 5000*1000:
        size_lable = "L)>5MB"
    else:
        pdb.set_trace()
         
    print([reader[0:5],GC,MAP, cnv_size, size_lable])

    result_list.append([reader[0],reader[1],reader[2],reader[3],reader[4],reader[5],
        reader[6],reader[7],reader[8],reader[9],GC,reader[11],MAP,reader[13],size_lable])

output_to_file(result_list, file_name = result_file)
