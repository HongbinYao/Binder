# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 17:16:38 2022

@author: yhb
"""
import os
from copy import deepcopy
from write_gtf import write_gtf

def filter_gtf(params_output, params_genome):
    pre_file = open(params_output,'r')
    pre_line = [line.rstrip().split('\t') for line in pre_file]; pre_trans = {}; pre_file.close()
    
    tran = pre_line[0][8].split(' ')[3].replace('"','').replace(';', '')
    
    for line in pre_line:
        if line[6] == '.':
            continue
        curr_tran = line[8].split(' ')[3].replace('"','').replace(';', '')
        if curr_tran != tran:
            pre_trans[chrom].update(temp)
            
        if line[2] == "transcript":
            tran = line[8].split(' ')[3].replace('"','').replace(';', '')
            chrom = line[0] + line[6]; 
            temp = {}; seq = []
        else:
            if line[2] == "exon":
                seq.append(int(line[3]))
                seq.append(int(line[4]))
        
        if chrom not in pre_trans:
            pre_trans[chrom] = {}

        if tran not in temp:
            temp[tran] = seq
        else:
            temp[tran] = seq

    pre_trans[chrom].update(temp)#添加最后一个transcript
            
    os.system('chmod +x gffread/gffread')
    os.system('./gffread {params_output} -g {params_genome} -w transcripts.fa')
    
    os.system('chmod +x kallisto/kallisto')
    os.system('./kallisto/estimate.sh')
    
    tsv_file = open('kallisto/out_put/abundance.tsv','r')
    tsv_line = [line.rstrip().split('\t') for line in tsv_file]; remove_trans = {}; tsv_file.close()
#*******提取要去除的tran_id*********
    for line in tsv_line[1::]:
        if float(line[-1]) <= 0.4:
            remove_trans[line[0]] = float(line[-1])
            
#*******筛掉要去除的tran_id*********
    rectify_trans = {}
    for chrom, trans in pre_trans.items():
        rectify_trans[chrom] = {}
        for tran, exons in trans.items():
            if tran not in remove_trans:
                rectify_trans[chrom].update({tran:exons})

#*******写入rectify_trans**********
    out_put_file = open(params_output,'w')
    write_gtf(rectify_trans, out_put_file)
    
    
def pre_est():
    os.system("chmod +x gffread")
    os.system("gffread Tiglon.gtf -g genome.fa -w transcripts.fa")
    os.system("chmod +x kallisto.sh")
    os.system("./kallisto.sh")
    
    path = r'C:\Users\yhb\Desktop\gtf文件\gtf\\'
    pre_file  = open(path + 'Binder.gtf','r')
    tsv_file = open(path + 'abundance.tsv','r')
    out_put_file = open(path + 'new.gtf','w')
    
    tsv_line = tsv_file.readlines()
    
    tsv_row = [line.rstrip().split('\t') for line in tsv_line]; tsv_row.pop(0)
    
    tpm = {}
    for line in tsv_row:
        tpm[line[0]] = float(line[-1])
    
    return tpm
    
def filter_genome_guided(chrom, tpm, merge_write):
    copy = deepcopy(merge_write)
    value = []
    for k, v in tpm.items():
        value.append(v)
    value.sort()
    median = value[int(len(value)/2)]
    for region in copy[::-1]:
        for transcript in region[::-1]:
            if tpm[transcript[-2]] <= median:
                region.remove(transcript)
        if not region:
            copy.remove(region)
    
    return copy
    # for region in copy:
    #     if len(region) >= 2:
    #         value = []
    #         for transcript in region:
    #             value.append(tpm[transcript[-2]])
    #         high = max(value)
            #只保留一个转录本
            
                
    # for k, v in tpm.items():
    #     if k.split('.')[0] == "chr" + str(chrom):
    #         value.append(v)
    #         trans.append(k.split('.')[1])
            
# if __name__ == '__main__':
#     for k, v in tpm.items():
#         if k.split('.')[0] == "chr" + str(10):
#             value.append(v)
#             trans.append(k.split('.')[1])     
# output_file = open(path + 'remove_id.txt','w')
# re = {}
# value = []

# for line in items:
#     re[line[0]] = float(line[3])
#     value.append(float(line[3]))

# value.sort()
# select = []

# for k, v in re.items():
#     if v >= 531.746:
#         select.append(k)
        
# for region in merge_write:
#     if len(region)>=2:
#         t = []
#         for i in region:
#             t.append(i[])

        
        
# tsv_file.close()
# output_file.close()


# print("filter transtripts which its' id not in remove_id.txt!\n")

# path = r'C:\Users\yhb\Desktop\gtf文件\\'
# de_nove_file = open(path + r'gtf\new.gtf','r')
# # remove_file = open(path + r'tsv\remove_id.txt','r')
# output_file = open(path + r'tsv\new.gtf','w')

# file1_seq = []

# de_nove_line = de_nove_file.readlines()
# # remove_line = remove_file.readlines()


# de_nove_row = [line.rstrip().split('\t') for line in de_nove_line]
# # remove_row = [line.rstrip() for line in remove_line]

    
# def calculate(chain,select):
#     print("Start to calculate the chain: " + chain + "\n")
#     for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M'}:
#         seq1 = []
#         for line in de_nove_row:
#             if line[0] == 'chr' + str(i) and line[6] == chain:
#                 line[8] = line[8].rstrip()#bridger.gtf每行末尾含有转义字符'\t',需删除
#                 seq1.append(line)
#         ##提取de nove 各个chr的转录本分别保存为一个list(区分正负链)
         
#         re1 = {}
        
#         for line in seq1:
#             seq = []
#             if line[2] == 'transcript':
#                 line[8] = line[8].rstrip()#含有'transcript'行末尾含有转义字符'\t',需删除
#                 id = line[8].split(' ')[3].replace('"','').replace(';', '')
#                 seq.append(line)
#             else:
#                 line[8] = line[8].rstrip()
#                 seq.append(line)
            
#             if id not in re1:
#                 re1[id] = seq
#             else:
#                 re1[id] += seq
#         ###提取上述de nove list中的每个转录本，并保存为一个字典，key为转录本的起点和终点
        
        
#         for k1, v1 in re1.items():
#             if k1 in select:
#                 for line in v1:
#                     for ele in line:
#                         output_file.write(ele + '\t')
#                     output_file.write('\n')
                    
#         # for k1, v1 in re1.items():
#         #     if k1 in remove_row:
#         #         for line in v1:
#         #             for ele in line:
#         #                 output_file.write(ele + '\t')
#         #             output_file.write('\n')

                
# for chain in {'+', '-'}:
#     calculate(chain,select)
           
# de_nove_file.close()
# # remove_file.close()
# output_file.close()

# print("The work is finished!")