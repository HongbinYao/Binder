# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 21:22:25 2022

@author: yhb
"""
from copy import deepcopy
from read_bam import read_bam

def exon_similar(merge_temp, bam_file):
    merge_write = []
    copy = deepcopy(merge_temp) 
    for region in copy:
        temp = []
        genome = []
        for tran in region[::-1]:
            if len(tran) > 2*2 + 2 and tran[-1] == "genome_guided":
                for exon in tran[:-2]:
                    genome.append(exon)
        genome = set(genome)
        
        while len(region) > 1 and len(region[-1]) == 2*2 + 2:
            temp.append(region.pop(-1))
            choose = False
            left = temp[-1][1];right = temp[-1][2]
            
            for tran in region[::-1]:
                if len(tran) == 2*2 + 2:
                    
                    left_temp = tran[1];right_temp = tran[2]
                    if left == left_temp and right == right_temp:
                        choose = True
                        region.remove(tran)
            
            if choose == False:
                if left in genome and right in genome:
                    choose = True
                elif read_bam(bam_file, chrom, start, end) == True:
                    choose = True
                else:
                    temp.pop(-1)
        else:
            while len(region) == 1 and len(region[-1]) == 2*2 + 2 and region[-1][-1] == "genome_guided":
                temp.append(region.pop(-1))
            else:
                merge_write.append(temp)
            
    return merge_write

def merge_similar(merge_trans, bam_file):
    merge_write = {}
    for chrom, trans in deepcopy(merge_trans).items():
        merge_write[chrom] = []
        for region in trans:
            temp = []
            genome_guided_temp = []; genome = []
            for tran in region[::-1]:
                if len(tran) > 2*2 + 2 and tran[-1] == 'genome_guided':
                    genome_guided_temp.append(tran[0:-2])
                    for exon in tran[:-2]:
                        genome.append(exon)
            genome = set(genome)
            
            while len(region) >= 2 and len(region[-1]) == 2*2 + 2:
                temp.append(region.pop(-1))
                choose = False
                left = temp[-1][1]; right = temp[-1][2]
                
                for tran in region[::-1]:
                    if len(tran) == 2*2 + 2:
                        left_temp = tran[1]; right_temp = tran[2]
                        if left == left_temp and right == right_temp:
                            choose = True
                            region.remove(tran)
                
                if choose == False:
                    if read_bam(bam_file, chrom[0:-1], left, right) == True:
                        choose = True
                    else:
                        temp.pop(-1)
            else:
                if len(region) >= 1 and len(region[-1]) == 2*2 + 2:
                    left = region[-1][1]; right = region[-1][2]
                    if read_bam(bam_file, chrom[0:-1], left, right) == True:
                        if left not in genome or right not in genome:
                            temp.append(region.pop(-1))
                        else:
                            choose = False
                            for tran in genome_guided_temp:
                                if tran[0] == left and tran[-1] == right:
                                   choose = True
                            if choose == True:
                                temp.append(region.pop(-1))
                    if temp:
                        for tran in temp:
                            tran.pop(-1); tran.pop(-1)
                        merge_write[chrom].append(temp)
                else:
                    if temp:
                        for tran in temp:
                            tran.pop(-1); tran.pop(-1)
                    merge_write[chrom].append(temp)
            
    return merge_write

def de_nove_similar(merge_trans, bam_file):
    de_nove_write = {}
    for chrom, trans in deepcopy(merge_trans).items():
        de_nove_write[chrom] = []
        for region in trans:
            temp = []
            genome = []
            for tran in region[::-1]:
                if len(tran) > 2*2 + 2 and tran[-1] == "genome_guided":
                    for exon in tran[:-2]:
                        genome.append(exon)
                elif len(tran) == 2*2 + 2 and tran[-1] == "genome_guided":
                    region.remove(tran)
            genome = set(genome)
            
            while len(region) >= 2 and len(region[-1]) == 2*2 + 2:
                temp.append(region.pop(-1))
                choose = False
                left = temp[-1][1]; right = temp[-1][2]
                
                for tran in region[::-1]:
                    if len(tran) == 2*2 + 2:
                        left_temp = tran[1]; right_temp = tran[2]
                        if left == left_temp and right == right_temp:
                            choose = True
                            region.remove(tran)
                
                if choose == False:
                    if left in genome and right in genome: #or temp[-1][-1] == "genome_guided"
                        choose = True
                    elif read_bam(bam_file, chrom[0:-1], left, right) == True:
                        choose = True
                    else:
                        temp.pop(-1)
            else:
                if len(region) == 1 and len(region[-1]) == 2*2 + 2:#and region[-1][-1] == 'genome_guided'
                    left = region[-1][1]; right = region[-1][2]
                    if read_bam(bam_file, chrom[0:-1], left, right) == True:
                        temp.append(region.pop(-1))
                        for tran in temp:
                            tran.pop(-1); tran.pop(-1)
                        de_nove_write[chrom].append(temp)
                else:
                    if temp:
                        for tran in temp:
                            tran.pop(-1); tran.pop(-1)
                    de_nove_write[chrom].append(temp)
            
    return de_nove_write


