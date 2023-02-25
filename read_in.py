# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 15:42:54 2022

@author: yhb
"""

from write_in import write_in

path = r'C:\Users\yhb\Desktop\gtf文件\gtf\\'
pre_file = open(path + 'new.gtf','r')
processed_file = open(path + 'Binder.gtf','w')

file1_seq = []
file2_seq = []

pre_line = pre_file.readlines()


pre_row = [line.rstrip().split('\t') for line in pre_line]


tsv_file = open(r'C:\Users\yhb\Desktop\gtf文件\tsv\\' + 'abundance.tsv','r')
lines = []
items = []
text = tsv_file.readlines()
items = [line.rstrip().split('\t') for line in text]
items.pop(0)
re = {}
seq = []
for line in items:
    seq.append(float(line[3]))
    
seq.sort()

median = seq[int(len(seq)/2)]

for line in items:
    
    re[line[0]] = float(line[3])

for chrom in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M"}:
    print("Start to calculate the chrom: " + 'str' + str(chrom) + "\n")
    gene_index = 0
    for chain in {'+', '-'}:
        seq1 = []
        seq2 = []
        for line in pre_row:
            if line[0] == "chr" + str(chrom) and line[6] == chain:
                line[8] = line[8].rstrip()#bridger.gtf每行末尾含有转义字符'\t',需删除
                seq1.append(line)
        ###提取de nove 各个chr的转录本分别保存为一个list(区分正负链)
        
        re1 = {}
        seq = []
        
        for line in seq1:
            
            if line[2] == "transcript":
                line[8] = line[8].rstrip()#含有'transcript'行末尾含有转义字符'\t',需删除
                key = line[8].split(' ')[3].replace('"','').replace(';', '')
                seq = []
            else:
                if line[2] == "exon":
                    line[8] = line[8].rstrip()
                    seq.append(int(line[3]))
                    seq.append(int(line[4]))
            
            if key not in re1:
                re1[key] = seq
            else:
                re1[key] = seq
                
        re1 = {key:val for key, val in re1.items() if len(val)/2 > 1}#删除exon个数为1的转录本 
    
        pre_exon = dict(sorted(re1.items(),key=lambda x:x[1]))
        
        key = list(pre_exon.keys())
        value = list(pre_exon.values())
        if len(value) == 1:
                value[0].append(key[0])
                value[0].append("de_nove")
        else:
            for i in range(len(value)):
                value[i].append(key[i])
                value[i].append("de_nove")
        pre_transcripts = value
        
        pre_temp = []
        temp = []
        if pre_transcripts:
            temp.append(pre_transcripts.pop(0))
            start = temp[0][1]
            end = temp[0][-4]
            # merge.pop(0)
            
            
            for transcript in pre_transcripts:
                start_temp = transcript[1]
                end_temp = transcript[-4]
                if max(start, start_temp) < min(end, end_temp):
                    temp.append(transcript)
                    start = min(start, start_temp)
                    end = max(end, end_temp)
                else:
                    pre_temp.append(temp)
                    temp = []
                    temp.append(transcript)
                    start = temp[0][1]
                    end = temp[0][-4]
                    
            if pre_transcripts.index(transcript) == len(pre_transcripts) - 1:
                pre_temp.append(temp)
            
            for region in pre_temp:
                for i in region[::-1]:
                    if re[i[-2]] < median:
                        # print(i)
                        region.remove(i)
                        
            gene_index = write_in(chrom, chain, gene_index, processed_file, pre_temp)
    
processed_file.close()
pre_file.close()
tsv_file.close()
