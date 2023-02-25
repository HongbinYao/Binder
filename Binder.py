# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 19:51:23 2022

@author: yhb
"""
# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import os
from calculate import read_gtf, pre_process
from cluster import merge_cluster
from exon_similar import merge_similar, de_nove_similar
from exon_rectify import merge_rectify, intersect_rectify, de_nove_rectify
from write_gtf import write_merge
from abundance_estimation import filter_gtf
from utils.os_utils import *
from extract_tran import extract_intersect, extract_genome_guided, \
                         construct_genome_guided
                         

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--de_nove', '-d',
                        help='path of gtf file about de_nove assembler',
                        )
    parser.add_argument('--genome_guided', '-g',
                        help='path of gtf file about genome_guided assembler',
                        )
    parser.add_argument('--bam', '-B',
                        help='BAM file',
                        )
    parser.add_argument('--genome', '-G',
                        help='Genome fasta to generate gtf fasta file',
                        )
    parser.add_argument('--read1', '-r1',
                        help='fasta file of Paired-end left-read',
                        )
    parser.add_argument('--read2', '-r2',
                        help='fasta file of Paired-end right-read',
                        )
    parser.add_argument('--output', '-o',
                        help='Output path',
                        )
    parser.add_argument('--h',
                        help='help infomation',
                        type=int)
    
    params = parser.parse_args()
    
    return params


def Binder():
    params = parse_args()
    
    print('Start to read de_nove gtf file and genome_guided gtf file', '\n')
    de_nove_trans, genome_guided_trans = read_gtf(params.de_nove, params.genome_guided)
    
    print('Start to pre_process de_nove file and genome_guided file', '\n')
    de_nove_trans, genome_guided_trans = pre_process(de_nove_trans, genome_guided_trans)
    
    print('Start to cluster', '\n')
    merge_trans = merge_cluster(de_nove_trans, genome_guided_trans)
    
    print('Start to rectify exon_similar', '\n')
    merge_write = merge_similar(merge_trans)
    
    print('Start to rectify exon_merge', '\n')
    merge_write = merge_rectify(merge_trans, merge_write)
        
    print('write the rectify gtf file', '\n')
    out_put_file = open(params.output,'w')
    write_merge(merge_write, out_put_file)
    
    print('estimate the rectify gtf file', '\n')
    filter(params.output, params.genome)
    
    return None
# Binder()
params = parse_args()
print('Start to read de_nove gtf file and genome_guided gtf file', '\n')
de_nove_trans, genome_guided_trans = read_gtf(params.de_nove, params.genome_guided)
print('Start to pre_process de_nove file and genome_guided file', '\n')
de_nove_trans, genome_guided_trans = pre_process(de_nove_trans, genome_guided_trans)
print('Start to cluster', '\n')
merge_trans = merge_cluster(de_nove_trans, genome_guided_trans)
print('Start to rectify exon_similar', '\n')
merge_write = merge_similar(merge_trans, params.bam)
print('Start to rectify exon_merge', '\n')
merge_write = merge_rectify(merge_trans, merge_write)
print('write the rectify gtf file', '\n')
out_put_file = open(params.output,'w')
write_merge(merge_write, out_put_file)

if __name__ == "__main__":
    path = r'C:\Users\yaohongbin\Desktop\gtf\\'
    de_nove_file = path + 'bridger.gtf'
    genome_guided_file = path + 'Tiglon.gtf'
    out_put_file = open(path + r'test.txt','w')
    bam_file = path + 'STAR.bam'

    de_nove_trans, genome_guided_trans = read_gtf(de_nove_file, genome_guided_file)
    de_nove_trans, genome_guided_trans = pre_process(de_nove_trans, genome_guided_trans)
    merge_trans = merge_cluster(de_nove_trans, genome_guided_trans)
    merge_write = merge_similar(merge_trans, bam_file)
    merge_write = merge_rectify(merge_trans, merge_write)
    write_merge(merge_write, out_put_file)
    
#*********intersect*********
    intersect_write = intersect_rectify(merge_trans)
    write_merge(intersect_write, out_put_file)
    
#*********de_nove***********
    de_nove_write = de_nove_similar(merge_trans, bam_file)
    de_nove_write = de_nove_rectify(merge_trans, de_nove_write, bam_file)
    write_merge(de_nove_write, out_put_file)

#*********提取与de_nove有交集的 genome_guided的转录本并构图***********
    construct_write = extract_genome_guided(merge_trans)
    write_merge(construct_write, out_put_file)


    merge_write = {}
    for chrom, trans in merge_trans.items():
        merge_write[chrom] = []
        for region in trans:
            de_nove_temp = []
            genome_guided_temp = []
            for tran in region:
                if len(tran) > 2*2 + 2 and tran[-1] == "genome_guided":
                    genome_guided_temp.append(tran)
                elif len(tran) > 2*2 + 2 and tran[-1] == "de_nove":
                    de_nove_temp.append(tran)
                    
            prime = []; s_prime = []; e_prime = []
            for tran in genome_guided_temp:
                prime.append((tran[1],tran[-4]))
                s_prime.append(tran[1]); e_prime.append(tran[-4])
            prime = list(set(prime)); prime.sort()
            s_prime = list(set(s_prime)); s_prime.sort()
            e_prime = list(set(e_prime)); e_prime.sort()
            
            splice_junction = []
            for tran in genome_guided_temp:
                for i in tran[1:-3]:
                    splice_junction.append(i)
            splice_junction = list(set(splice_junction)); splice_junction.sort()
            
            de_nove_exon = []
            for tran in de_nove_temp:
                temp = []
                for index in range(0,len(tran)-2, 2):
                    exon = ()
                    exon += (tran[index],); exon += (tran[index + 1],)
                    temp.append(exon)
                de_nove_exon.append(temp)
            temp = []
            if genome_guided_temp and de_nove_temp:
                for tran in de_nove_exon:
                    append = False
                    for exon in tran[1:-1]:
                        if exon[0] in e_prime or exon[-1] in s_prime:
                            index = 0; num = 0; total = len(tran) - 1
                            while index < len(tran) - 1:
                                if tran[index][-1] in splice_junction and tran[index+1][0] in splice_junction:
                                    num += 1
                                index += 1
                            if num/total >= 0.5:
                                append = True
                            
                    if append == True:
                        temp.append(de_nove_temp[de_nove_exon.index(tran)][0:-2])
                merge_write[chrom].append(temp)
            
    write_merge(merge_write, out_put_file)
