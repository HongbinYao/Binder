# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 21:32:51 2022

@author: yaohongbin
"""
from copy import deepcopy
from adjacency_list import construct_graph, rectify, append_exon
from path_search import path_search, specific_path_search
from assembly import assembly, assembly_sr

def merge_rectify(merge_trans, merge_write):
    # merge_write = {}
    for chrom, trans in merge_trans.items():
        # merge_write[chrom] = []; 
        index = 0
        for region in trans:
            de_nove_temp = []
            genome_guided_temp = []
            for tran in region:
                if len(tran) > 2*2 + 2 and tran[-1] == "genome_guided":
                    genome_guided_temp.append(tran)
                elif len(tran) > 2*2 + 2 and tran[-1] == "de_nove":
                    de_nove_temp.append(tran)
            
            if genome_guided_temp and de_nove_temp:
                splice_graph = deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
                de_nove_rectify, not_append_end = deepcopy(rectify(de_nove_temp, splice_graph))#对de_nove exon 进行修正
                splice_graph = deepcopy(append_exon(de_nove_rectify, not_append_end, splice_graph))#添加修正后的de_nove exon
                path = path_search(splice_graph)
                
                temp = []
                for tran in path:
                    path_temp = []
                    for exon in tran:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
                    
                merge_write[chrom][index] += temp; index += 1
            elif genome_guided_temp and not de_nove_temp:
                merge_write[chrom][index] += deepcopy(genome_guided_temp); index += 1 
                for tran in merge_write[chrom][index-1]:
                    tran.pop(-1); tran.pop(-1)
            elif not genome_guided_temp and de_nove_temp:#后续进一步分析
                merge_write[chrom][index] += deepcopy(de_nove_temp); index += 1
                for tran in merge_write[chrom][index-1]:
                    tran.pop(-1); tran.pop(-1)
            else:#not genome_guided_temp and not de_nove_temp
                index += 1

    return merge_write

def intersect_rectify(merge_trans):
    merge_write = {}
    for chrom, trans in merge_trans.items():
        merge_write[chrom] = []; 
        for region in trans:
            de_nove_temp = []
            genome_guided_temp = []
            for tran in region:
                if len(tran) > 2*2 + 2 and tran[-1] == "genome_guided":
                    genome_guided_temp.append(tran)
                elif len(tran) > 2*2 + 2 and tran[-1] == "de_nove":
                    de_nove_temp.append(tran)
            
            if genome_guided_temp and de_nove_temp:
                splice_graph = deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
                de_nove_rectify, not_append_end = deepcopy(rectify(de_nove_temp, splice_graph))#对de_nove exon 进行修正
                splice_graph = deepcopy(append_exon(de_nove_rectify, not_append_end, splice_graph))#添加修正后的de_nove exon
                path = specific_path_search(splice_graph)
                
                temp = []
                for tran in path:
                    path_temp = []
                    for exon in tran:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
                    
                merge_write[chrom].append(temp)
                
    return merge_write

def de_nove_rectify(merge_trans, de_nove_write, bam_file):
    # merge_write = {}
    for chrom, trans in merge_trans.items():
        # merge_write[chrom] = []
        index = 0
        for region in trans:
            de_nove_temp = []
            genome_guided_temp = []
            for tran in region:
                if len(tran) > 2*2 + 2 and tran[-1] == "genome_guided":
                    genome_guided_temp.append(tran)
                elif len(tran) > 2*2 + 2 and tran[-1] == "de_nove":
                    de_nove_temp.append(tran)
            
            if genome_guided_temp and de_nove_temp:
                de_nove_rectify, pair_end = assembly(chrom[0:-1], bam_file, de_nove_temp, genome_guided_temp)
                # splice_graph = deepcopy(construct_graph(de_nove_rectify))
                temp = []
                # if pair_end:
                #     path = specific_path_search(splice_graph, pair_end)
                #     for tran in path:
                #         path_temp = []
                #         for exon in tran:
                #             path_temp.append(exon[0]); path_temp.append(exon[1])
                #         temp.append(path_temp)
                    
                #此时只考虑修正后的de_nove exon 不进行path_search
                for tran in de_nove_rectify:
                    path_temp = []
                    for exon in tran:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
                    
                de_nove_write[chrom][index] += temp; index += 1
            elif genome_guided_temp and not de_nove_temp:
                index += 1
            elif not genome_guided_temp and de_nove_temp:
                de_nove_rectify, pair_end = assembly_sr(chrom[0:-1], bam_file, de_nove_temp)
                # splice_graph = deepcopy(construct_graph(de_nove_rectify))
                temp = []
                # path = specific_path_search(splice_graph, pair_end)
                # for tran in path:
                #     path_temp = []
                #     for exon in tran:
                #         path_temp.append(exon[0]); path_temp.append(exon[1])
                #     temp.append(path_temp)
                for tran in de_nove_rectify:
                    path_temp = []
                    for exon in tran:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
                de_nove_write[chrom][index] += temp; index += 1
            # else:#not genome_guided_temp and not de_nove_temp
            #     index += 1
            
    return de_nove_write

if __name__ == '__main__':
    import copy
    from adjacency_list import *
    from path_search import *
    from write_gtf import *
    region = copy.deepcopy(merge_temp[4])
    
    de_nove_temp = []
    genome_guided_temp = []
    for transcript in region:
        if len(transcript) > 2*2 + 2 and transcript[-1] == "genome_guided":
            genome_guided_temp.append(transcript)
        elif len(transcript) > 2*2 + 2 and transcript[-1] == "de_nove":
            de_nove_temp.append(transcript)
            
    splice_graph = copy.deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
    de_nove_rectify, not_append_end = copy.deepcopy(rectify(de_nove_temp, splice_graph))#对de_nove exon 进行修正
    splice_graph = copy.deepcopy(append_de_nove_exon(de_nove_rectify, not_append_end, splice_graph))#添加修正后的de_nove exon
    # path = path_search(splice_graph)
    path = simply_path_search(splice_graph)
    
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
    
#*************写入纠正后的de_nove exon*****************    
    temp = []
    for transcript in de_nove_rectify:
        path_temp = []
        for exon in transcript:
            path_temp.append(exon[0]); path_temp.append(exon[1])
        temp.append(path_temp)
    
    test_file = open(r'C:\Users\yhb\Desktop\gtf文件\gtf\new.gtf','w')
    test_write=[]
    test_write.append(temp)
    write_rectify(chrom, chain, 0, test_file, test_write)
    test_file.close()