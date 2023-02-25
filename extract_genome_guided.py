# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 21:45:41 2022

@author: yhb
"""
import copy

def extract_genome_guided(merge_temp):
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
            merge_write.append(genome_guided_temp)
            
    return merge_write

def extract_only_genome_guided(merge_temp):
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
        
        if genome_guided_temp and not de_nove_temp:
            merge_write.append(genome_guided_temp)
            
    return merge_write

def extract_all_genome_guided(merge_temp):
    merge_write = []
    deepcopy = copy.deepcopy(merge_temp)
    for region in deepcopy:
        genome_guided_temp = []
        for transcript in region:
            if transcript[-1] == "genome_guided":
                genome_guided_temp.append(transcript)
        
        merge_write.append(genome_guided_temp)
            
    return merge_write