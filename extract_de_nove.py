# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:24:58 2022

@author: yhb
"""
import copy
from adjacency_list import *
from path_search import *


def extract_de_nove(merge_temp):
    merge_write = []
    deepcopy = copy.deepcopy(merge_temp)
    for region in deepcopy:
        de_nove_temp = []
        genome_guided_temp = []
        for transcript in region:
            if len(transcript) > 2*2 + 2 and transcript[-1] == "genome_guided":
                genome_guided_temp.append(transcript)
            elif len(transcript) > 2*2 + 2 and transcript[-1] == "de_nove":
                de_nove_temp.append(transcript)
        
        if genome_guided_temp and de_nove_temp:
            splice_graph = copy.deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
            de_nove_rectify, not_append_end = copy.deepcopy(rectify(de_nove_temp, splice_graph))#对de_nove exon 进行修正
            genome_guided_end_point = copy.deepcopy(splice_graph.end_point)
            
            splice_graph = copy.deepcopy(append_exon(de_nove_rectify, not_append_end, splice_graph))#添加修正后的de_nove exon
            # end_point = {transcript[0]:[transcript[-1]] for transcript in de_nove_rectify if de_nove_rectify.index(transcript) not in not_append_end}
            end_point = copy.deepcopy(splice_graph.end_point)
            
            only_ends ={}
            for k, v in end_point.items():
                if k not in genome_guided_end_point:
                    only_ends[k] = v
                else:
                    for exon in v:
                        if exon not in genome_guided_end_point[k]:
                            if k not in only_ends:
                                temp = []; temp.append(exon)
                                only_ends[k] = temp
                            else:
                                temp = []; temp.append(exon)
                                only_ends[k] += temp
            if only_ends:
                path = specific_path_search(splice_graph, only_ends)
                temp = []
                for transcript in path:
                    path_temp = []
                    for exon in transcript:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
                merge_write.append(temp)
         
    return merge_write

def extract_only_de_nove(merge_temp):
    merge_write = []
    deepcopy = copy.deepcopy(merge_temp)
    for region in deepcopy:
        de_nove_temp = []
        genome_guided_temp = []
        for transcript in region:
            if len(transcript) > 2*2 + 2 and transcript[-1] == "genome_guided":
                genome_guided_temp.append(transcript)
            elif len(transcript) > 2*2 + 2 and transcript[-1] == "de_nove":
                de_nove_temp.append(transcript)
        
        if not genome_guided_temp and de_nove_temp:
            merge_write.append(de_nove_temp)
            
    return merge_write

def extract_de_nove_rectify(merge_temp):
    merge_write = []
    deepcopy = copy.deepcopy(merge_temp)
    for region in deepcopy:
        de_nove_temp = []
        genome_guided_temp = []
        for transcript in region:
            if len(transcript) > 2*2 + 2 and transcript[-1] == "genome_guided":
                genome_guided_temp.append(transcript)
            elif len(transcript) > 2*2 + 2 and transcript[-1] == "de_nove":
                de_nove_temp.append(transcript)
        
        if genome_guided_temp and de_nove_temp:
            splice_graph = copy.deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
            de_nove_rectify, not_append_end = copy.deepcopy(rectify(de_nove_temp, splice_graph))#对de_nove exon 进行修正
            temp = []
            for transcript in de_nove_rectify:
                if transcript:
                    path_temp = []
                    for exon in transcript:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
            merge_write.append(temp)
            
    return merge_write

def extract_end_point(merge_temp):
    merge_ends = {}
    only_ends = {}
    deepcopy = copy.deepcopy(merge_temp)
    for region in deepcopy:
        de_nove_temp = []
        genome_guided_temp = []
        for transcript in region:
            if len(transcript) > 2*2 + 2 and transcript[-1] == "genome_guided":
                genome_guided_temp.append(transcript)
            elif len(transcript) > 2*2 + 2 and transcript[-1] == "de_nove":
                de_nove_temp.append(transcript)
        
        if genome_guided_temp and de_nove_temp:
            splice_graph = copy.deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
            de_nove_rectify, not_append_end = copy.deepcopy(rectify(de_nove_temp, splice_graph))#对de_nove exon 进行修正
            genome_guided_end_point = copy.deepcopy(splice_graph.end_point)
            
            splice_graph = copy.deepcopy(append_de_nove_exon(de_nove_rectify, not_append_end, splice_graph))#添加修正后的de_nove exon
            # print(de_nove_rectify)
            de_nove_end_point = {transcript[0]:[transcript[-1]] for transcript in de_nove_rectify if de_nove_rectify.index(transcript) not in not_append_end}
            
            for k, v in de_nove_end_point.items():
                if k not in genome_guided_end_point:
                    only_ends[k] = v
                else:
                    for exon in v:
                        if exon not in genome_guided_end_point[k]:
                            if k not in only_ends:
                                temp = []; temp.append(exon)
                                only_ends[k] = temp
                            else:
                                temp = []; temp.append(exon)
                                only_ends[k] += temp
            if only_ends:
                path = specific_path_search(splice_graph, only_ends)
                
                temp = []
                for transcript in path:
                    path_temp = []
                    for exon in transcript:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
                merge_ends.update(only_ends)
    # merge_ends = dict(sorted(merge_ends.items(),key=lambda x:x[1]))
    return merge_ends

if __name__ == '__main__':
    import copy
    from write_gtf import *
    from adjacency_list import *
    region = copy.deepcopy(merge_temp[18])
    
    de_nove_temp = []
    genome_guided_temp = []
    for transcript in region:
        if len(transcript) > 2*2 + 2 and transcript[-1] == "genome_guided":
            genome_guided_temp.append(transcript)
        elif len(transcript) > 2*2 + 2 and transcript[-1] == "de_nove":
            de_nove_temp.append(transcript)
            
    splice_graph = copy.deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
    de_nove_rectify, not_append_end = copy.deepcopy(rectify(de_nove_temp, splice_graph))#对de_nove exon 进行修正
    genome_guided_end_point = copy.deepcopy(splice_graph.end_point)
    splice_graph = copy.deepcopy(append_de_nove_exon(de_nove_rectify, not_append_end, splice_graph))#添加修正后的de_nove exon
    de_nove_end_point = {transcript[0]:[transcript[-1]] for transcript in de_nove_rectify if de_nove_rectify.index(transcript) not in not_append_end}
    
    only_ends ={}
    for k, v in de_nove_end_point.items():
        if k not in genome_guided_end_point:
            only_ends[k] = v
        else:
            for exon in v:
                if exon not in genome_guided_end_point[k]:
                    if k not in only_ends:
                        temp = []; temp.append(exon)
                        only_ends[k] = temp
                    else:
                        temp = []; temp.append(exon)
                        only_ends[k] += temp
                        
    path = specific_path_search(splice_graph, only_ends)
    
    temp = []
    for transcript in path:
        path_temp = []
        for exon in transcript:
            path_temp.append(exon[0]); path_temp.append(exon[1])
        temp.append(path_temp)
    
    test_file = open(r'C:\Users\yhb\Desktop\gtf文件\gtf\new.gtf','w')
    test_write=[]
    test_write.append(temp)
    write_rectify(chrom, chain, 0, test_file, test_write)
    test_file.close()