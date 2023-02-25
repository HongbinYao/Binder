# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 15:29:22 2022

@author: yhb
"""
from copy import deepcopy

def merge_cluster(de_nove_trans, genome_guided_trans):
    de_nove_temp = {}; genome_guided_temp = {}; merge_trans = {}
    
    for chrom, trans in de_nove_trans.items():
        de_nove_temp[chrom] = []
        for tran_id, tran in trans.items():
            copy = deepcopy(tran); copy.append(tran_id); copy.append('de_nove')
            de_nove_temp[chrom].append(copy)
    
    for chrom, trans in genome_guided_trans.items():
        genome_guided_temp[chrom] = []
        for tran_id, tran in trans.items():
            copy = deepcopy(tran); copy.append(tran_id); copy.append('genome_guided')
            genome_guided_temp[chrom].append(copy)
    
    for chrom in genome_guided_temp:
        merge_trans[chrom] = deepcopy(genome_guided_temp[chrom])
        if chrom in de_nove_temp:
            merge_trans[chrom] += deepcopy(de_nove_temp[chrom])
            
    for chrom, trans in merge_trans.items():
        merge_trans[chrom] = sorted(trans, key=lambda x:x[1])#以第一个exon的右边界 排序
        merge = deepcopy(merge_trans[chrom]); total = []
        
        if merge:
            region = []#; region.append(merge.pop(0))#region 用于存储每个基因座上的转录本
            start = merge[0][1]; end = merge[0][-4]
        else:
            continue
        
        for tran in merge:
            start_temp = tran[1]
            end_temp = tran[-4]
            if max(start, start_temp) < min(end, end_temp):
                region.append(tran)
                start = min(start, start_temp)
                end = max(end, end_temp)
            else:
                total.append(region)
                region = []; region.append(tran)
                start = tran[1]
                end = tran[-4]
        
        # if merge[-1] not in total[-1]:
        #     total.append(region)
        # else:
        total.append(region)
        
        for region in total:
            total[total.index(region)] = \
                sorted(region,key=lambda x:len(x),reverse = True)
                
        merge_trans[chrom] = total#把total赋值给各个chrom

    return merge_trans

def cluster(de_nove_trans, genome_guided_trans):
    de_nove_temp = []
    # temp = []
    key = deepcopy(list(de_nove_trans.keys()))
    value = deepcopy(list(de_nove_trans.values()))
    if len(value) == 1:
            value[0].append(key[0])
            value[0].append("de_nove")
    else:
        for i in range(len(value)):
            value[i].append(key[i])
            value[i].append("de_nove")
    de_nove_transcripts = value
    # temp.append(de_nove_transcripts[0])
    # start = temp[0][1]
    # end = temp[0][-4]
    # de_nove_transcripts.pop(0)
    
    
    # for de_nove_transcript in de_nove_transcripts:
    #     start_temp = de_nove_transcript[1]
    #     end_temp = de_nove_transcript[-4]
    #     if max(start, start_temp) < min(end, end_temp):
    #         temp.append(de_nove_transcript)
    #         start = min(start, start_temp)
    #         end = max(end, end_temp)
    #     else:
    #         de_nove_temp.append(temp)
    #         temp = []
    #         temp.append(de_nove_transcript)
    #         start = temp[0][1]
    #         end = temp[0][-4]
            
    # if not de_nove_transcripts:
    #     de_nove_temp.append(temp)
            
    genome_guided_temp = []
    # temp = []
    key = deepcopy(list(genome_guided_trans.keys()))
    value = deepcopy(list(genome_guided_trans.values()))
    if len(value) == 1:
            value[0].append(key[0])
            value[0].append("genome_guided")
    else:
        for i in range(len(value)):
            value[i].append(key[i])
            value[i].append("genome_guided")
    genome_guided_transcripts = value
    # temp.append(genome_guided_transcripts[0])
    # genome_guided_transcripts.pop(0)  
    
    merge = de_nove_transcripts + genome_guided_transcripts
    
    merge = sorted(merge,key=lambda x:x[1])#以第一个exon的右边界 排序
    # index = 0
    
    # start = de_nove_temp[index][0][1]
    # end = de_nove_temp[index][-1][-2]
    
    # for genome_guided_transcript in genome_guided_transcripts:
    #     start_temp = genome_guided_transcript[1]
    #     end_temp = genome_guided_transcript[-2]
    #     if max(start, start_temp) < min(end, end_temp):
    #         temp.append(genome_guided_transcript)
    #         start = min(start, start_temp)
    #         end = max(end, end_temp)
    #     else:        
    #         genome_guided_temp.append(temp)
    #         temp = []
    #         temp.append(genome_guided_transcript)
    #         index += 1
    #         if index >= len(de_nove_temp):
    #             index = 0
    #             break
    #             start = de_nove_temp[index][0][1]
    #             end = de_nove_temp[index][-1][-2]
    
    merge_temp = []
    temp = []
    temp.append(merge.pop(0))
    start = temp[0][1]
    end = temp[0][-4]
    # merge.pop(0)
    
    if merge:
        for transcript in merge:
            start_temp = transcript[1]
            end_temp = transcript[-4]
            if max(start, start_temp) < min(end, end_temp):
                temp.append(transcript)
                start = min(start, start_temp)
                end = max(end, end_temp)
            else:
                merge_temp.append(temp)
                temp = []
                temp.append(transcript)
                start = temp[0][1]
                end = temp[0][-4]
        
        if merge[-1] not in merge_temp[-1]:
            merge_temp.append(temp)
    else:
        merge_temp.append(temp)
    # if merge.index(transcript) == len(merge) - 1:
    #     merge_temp.append(temp)
    
    return merge_temp