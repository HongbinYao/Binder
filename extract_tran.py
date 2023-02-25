# -*- coding: utf-8 -*-
"""
Created on Tue May 31 21:54:39 2022

@author: yhb
"""
from copy import deepcopy
from adjacency_list import construct_graph
from path_search import path_search

def extract_intersect(merge_trans):
    merge_write = {}
    for chrom, trans in merge_trans.items():
        merge_write[chrom] = []; 
        for region in trans:
            de_nove_temp = []
            genome_guided_temp = []
            for tran in region:
                if tran[-1] == "genome_guided":
                    genome_guided_temp.append(tran)
                elif tran[-1] == "de_nove":
                    de_nove_temp.append(tran)
            
            if genome_guided_temp and de_nove_temp:
                merge_write[chrom].append(region)
    
    return merge_write


def construct_genome_guided(merge_trans):
    merge_write = {}
    for chrom, trans in merge_trans.items():
        merge_write[chrom] = []; 
        for region in trans:
            genome_guided_temp = []
            for tran in region:
                de_nove_temp = []
                genome_guided_temp = []
                for tran in region:
                    if tran[-1] == "genome_guided":
                        genome_guided_temp.append(tran)
                    elif tran[-1] == "de_nove":
                        de_nove_temp.append(tran)
            
            if genome_guided_temp:
                splice_graph = deepcopy(construct_graph(genome_guided_temp))#将genome_guided_exon添加到splice graph中
                path = path_search(splice_graph)
                
                temp = []
                for tran in path:
                    path_temp = []
                    for exon in tran:
                        path_temp.append(exon[0]); path_temp.append(exon[1])
                    temp.append(path_temp)
                    
                merge_write[chrom].append(temp)
            
    return merge_write


def extract_genome_guided(merge_trans):
    merge_write = {}
    for chrom, trans in merge_trans.items():
        merge_write[chrom] = []; 
        for region in trans:
            genome_guided_temp = []
            for tran in region:
                de_nove_temp = []
                genome_guided_temp = []
                for tran in region:
                    if tran[-1] == "genome_guided":
                        genome_guided_temp.append(tran)
                    elif tran[-1] == "de_nove":
                        de_nove_temp.append(tran)
            
            if genome_guided_temp:
                temp = []
                for tran in genome_guided_temp:
                    path_temp = []
                    for exon in tran[:-2:]:
                        path_temp.append(exon)
                    temp.append(path_temp)
                    
                merge_write[chrom].append(temp)
            
    return merge_write