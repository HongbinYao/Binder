# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 19:36:41 2022

@author: yhb
"""
from cluster import *
from exon_similar import *
from exon_rectify import *
from extract_de_nove import *
from extract_genome_guided import *
from write_gtf import *

path = r'C:\Users\yhb\Desktop\gtf文件\gtf\\'
de_nove_file = open(path + r'bridger-clean.gtf','r')
genome_guided_file = open(path + r'transcomb-clean.gtf','r')
merge_file = open(path + r'new.gtf','w')


# de_nove_line = de_nove_file.readlines()
# genome_guided_line = genome_guided_file.readlines()


de_nove_row = [line.rstrip().split('\t') for line in de_nove_file]
genome_guided_row = [line.rstrip().split('\t') for line in genome_guided_file]

chrom = 1
chain = '-'
print("Start to calculate the chrom: " + 'str' + str(chrom) + "\n")
gene_index = 0
seq1 = []
seq2 = []
for line in de_nove_row:
    if line[0] == "chr" + str(chrom) and line[6] == chain:
        line[8] = line[8].rstrip()#bridger.gtf每行末尾含有转义字符'\t',需删除
        seq1.append(line)
###提取de nove 各个chr的转录本分别保存为一个list(区分正负链)
for line in genome_guided_row:
    if line[0] == "chr" + str(chrom) and line[6] == chain:
        seq2.append(line)
###提取genome guided 各个chr的转录本分别保存为一个list(区分正负链)         

de_nove_trans = {}
genome_guided_trans = {}
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
    
    if key not in de_nove_trans:
        de_nove_trans[key] = seq
    else:
        de_nove_trans[key] = seq
        
de_nove_trans = {key:val for key, val in de_nove_trans.items() if len(val)/2 > 1}#删除exon个数为1的转录本 
###提取上述de nove list中的每个转录本，并保存为一个字典，key为转录本的起点和终点
for line in seq2:
    
    if line[2] == "transcript":
        line[8] = line[8].rstrip()#含有'transcript'行末尾含有转义字符'\t',需删除
        key = line[8].split(' ')[3].replace('"','').replace(';', '')
        seq = []
    else:
        if line[2] == "exon":
            line[8] = line[8].rstrip()
            seq.append(int(line[3]))
            seq.append(int(line[4]))
    
    if key not in genome_guided_trans:
        genome_guided_trans[key] = seq
    else:
        genome_guided_trans[key] = seq
  ###提取上述genome guided list中的每个转录本，并保存为一个字典，key为转录本的起点和终点
for key, val in genome_guided_trans.items():
    index = 0
    bound = len(val)
    while index <=  bound - 1 - 2:
        right_end = val[index + 1]; left_end = val[index + 2]
        if left_end - right_end <= 5:
            val.pop(index + 1); val.pop(index + 1)
            bound -= 2
        else:
            index += 2
##去除exon之间短的间隙, example 11134 12345 12346 12666
genome_guided_trans = {key:val for key, val in genome_guided_trans.items() if len(val)/2 > 1}#删除exon个数为1的转录本
            
merge_temp = cluster(de_nove_trans, genome_guided_trans)

###对每个基因座进行sort
for region in merge_temp:
    merge_temp[merge_temp.index(region)] = sorted(region,key=lambda x:len(x),reverse = True)

#提取rectify exon = 2 transcripts
# merge_write = exon_similar(merge_temp)
# gene_index = write_gtf(chrom, chain, gene_index, merge_file, merge_write)

#提取 all region transcripts
# merge_write = exon_similar(merge_temp)
# merge_write = exon_rectify(merge_temp, merge_write)
# gene_index = write_rectify(chrom, chain, gene_index, merge_file, merge_write)

#提取de_nove and genome_guided transcripts     
# merge_write = merge_rectify(merge_temp)
# gene_index = write_rectify(chrom, chain, gene_index, merge_file, merge_write)

#提取与genome_guided 有交集的 de_nove transcripts        
merge_write = extract_de_nove(merge_temp)
gene_index = write_rectify(chrom, chain, gene_index, merge_file, merge_write)

#提取de_nove end_point
# merge_ends = extract_end_point(merge_temp)
# write_merge_ends(chrom, chain, merge_file, merge_ends)

#提取与de_nove 有交集的 genome_guided transcripts
# merge_write = extract_genome_guided(merge_temp)
# gene_index = write_gtf(chrom, chain, gene_index, merge_file, merge_write)

#提取基因座只含有genome_guided transcripts
# merge_write = extract_only_genome_guided(merge_temp)
# gene_index = write_gtf(chrom, chain, gene_index, merge_file, merge_write)

#提取基因座只含有de_nove transcripts
# merge_write = extract_only_de_nove(merge_temp)
# gene_index = write_gtf(chrom, chain, gene_index, merge_file, merge_write)

#提取所有genome_guided transcripts
# merge_write = extract_all_genome_guided(merge_temp)
# gene_index = write_gtf(chrom, chain, gene_index, merge_file, merge_write)
        
#提取所有修正后的de_nove_rectify transcripts
# merge_write = extract_de_nove_rectify(merge_temp)
# gene_index = write_rectify(chrom, chain, gene_index, merge_file, merge_write)
        
        
de_nove_file.close()
genome_guided_file.close()
merge_file.close()

print("The work is finished!")