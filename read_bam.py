# -*- coding: utf-8 -*-
"""
Created on Mon May 16 21:44:59 2022

@author: yhb
"""

import numpy as np
from pysam import AlignmentFile

def read_bam(bam_file, chrom, start, end):
    bam = AlignmentFile(bam_file, 'r')
    
    splice_junction = {}
    for r in bam.fetch(chrom, start, end):
        data = str(r).split('\t')
        # print(data[3], data[5], data[7], data[8])
        # 121184274 76M 121184048 -302
        cigar = data[5];bound=[]
        read = []
        temp = ''
    
        for char in cigar:
            if str.isdigit(char):
                temp += char
            else:
                read.append(int(temp)); read.append(char); temp = ''
        
        junction = int(data[3])
        encounter = False; bound.append(junction)
        
        for index in range(0,len(read),2):
            if read[index+1] == 'M':
                junction += read[index] - 1
                if index+1 == len(read) - 1:
                    bound.append(junction)
                elif read[index+3] == 'N' or read[index+3] == 'S':
                    bound.append(junction)
            elif read[index+1] == 'N':
                junction += read[index] + 1
                bound.append(junction)
            elif read[index+1] == 'S' or read[index+1] == 'H':#若S出现在第一位不需要添加 因为起点自动加上了S的个数
                continue
            elif read[index+1] == 'I' or read[index+1] == 'P':#插入 即read多出I个碱基 过程中应忽略
                if encounter == False:
                    junction += 1; encounter = True
        
        # print(data[3], data[5], data[7], data[8], bound)
        
        if len(bound) == 4:
            junction = (bound[1],bound[2])
            if junction not in splice_junction:
                splice_junction[junction] = 1
            else:
                splice_junction[junction] += 1
                
    junction = (start,end); result = False#; print(splice_junction)
    if junction in splice_junction:
        if splice_junction[junction] >= 2:
            result = True#; print(True)
        else:#splice_junction[junction] <= 1
            temp = bam.count_coverage(chrom, start, end)
            coverage = [0 for i in range(0,end-start)]
            for index in range(0,len(coverage)):
                for row in temp:
                    coverage[index] += row[index]
                    
    # cigar = '39M67106N12M25S'
    # read = []
    # temp = ''
    
    # for char in cigar:
    #     if str.isdigit(char):
    #         temp += char
    #     else:
    #         read.append(int(temp)); read.append(char); temp = ''
            
    # temp = bam.count_coverage('chr1', 109508246, 109508248)
    # coverage = np.zeros((4,2))
    
    # for row in range(len(temp)):
    #     for col in range(len(temp[row])):
    #         coverage[row][col]=temp[row][col]
    
    bam.close()
    
    return result

def read_coverage(bam_file, chrom, start, end):
    bam = AlignmentFile(bam_file, 'r')
    temp = bam.count_coverage(chrom, start, end)
    
    coverage = [0 for i in range(0,end-start)]
    for index in range(0,len(coverage)):
        for row in temp:
            coverage[index] += row[index]
    aver_coverage = sum(coverage)/len(coverage)
    
    result = False
    if aver_coverage < 1:
        result = True
    
    bam.close()
    
    return result

def out_put_coverage(bam_file, chrom, start, end):
    bam = AlignmentFile(bam_file, 'r')
    temp = bam.count_coverage(chrom, start, end)
    
    coverage = [0 for i in range(0,end-start)]
    for index in range(0,len(coverage)):
        for row in temp:
            coverage[index] += row[index]
    aver_coverage = sum(coverage)/len(coverage)
    
    result = False
    if aver_coverage < 1:
        result = True
    
    bam.close()
    
    return coverage
        
def out_put_bam(bam_file, chrom, start, end):
    bam = AlignmentFile(bam_file, 'r')
    
    splice_junction = {}
    start = int(start); end = int(end)
    for r in bam.fetch(chrom, start, end):
        data = str(r).split('\t')
        cigar = data[5];bound=[]
        read = []
        temp = ''
    
        for char in cigar:
            if str.isdigit(char):
                temp += char
            else:
                read.append(int(temp)); read.append(char); temp = ''
        
        junction = int(data[3])
        encounter = False; bound.append(junction)
        
        for index in range(0,len(read),2):
            if read[index+1] == 'M':
                junction += read[index] - 1
                if index+1 == len(read) - 1:
                    bound.append(junction)
                elif read[index+3] == 'N' or read[index+3] == 'S':
                    bound.append(junction)
            elif read[index+1] == 'N':
                junction += read[index] + 1
                bound.append(junction)
            elif read[index+1] == 'S' or read[index+1] == 'H':#若S出现在第一位不需要添加 因为起点自动加上了S的个数
                continue
            elif read[index+1] == 'I' or read[index+1] == 'P':#插入 即read多出I个碱基 过程中应忽略
                if encounter == False:
                    junction += 1; encounter = True
        
        #print(data[3], data[5], data[7], data[8], bound)
        
        if len(bound) == 4:
            junction = (bound[1],bound[2])
            if junction not in splice_junction:
                splice_junction[junction] = 1
            else:
                splice_junction[junction] += 1
                
    junction = (start,end); result = False#; print(splice_junction)
    if junction in splice_junction:
        if splice_junction[junction] >= 2:
            result = True#; print(True)
    
    bam.close()
    
    return splice_junction

# result = read_bam('in.bam', 'chr1', 121181348, 121183458); print(result)

# if __name__ == "__main__":
#     bam = AlignmentFile('in.bam', 'r')

#     splice_junction = {}
#     for r in bam.fetch('chr1', 120913450, 120953220):
#         data = str(r).split('\t')
#         # print(data[3], data[5], data[7], data[8])
#         # 121184274 76M 121184048 -302
#         cigar = data[5];bound=[]
#         read = []
#         temp = ''

#         for char in cigar:
#             if str.isdigit(char):
#                 temp += char
#             else:
#                 read.append(int(temp)); read.append(char); temp = ''
        
#         junction = int(data[3])
#         encounter = False; bound.append(junction)
        
#         for index in range(0,len(read),2):
#             if read[index+1] == 'M':
#                 junction += read[index] - 1
#                 if index+1 == len(read) - 1:
#                     bound.append(junction)
#                 elif read[index+3] == 'N' or read[index+3] == 'S':
#                     bound.append(junction)
#             elif read[index+1] == 'N':
#                 junction += read[index] + 1
#                 bound.append(junction)
#             elif read[index+1] == 'S' or read[index+1] == 'H':#若S出现在第一位不需要添加 因为起点自动加上了S的个数
#                 continue
#             elif read[index+1] == 'I' or read[index+1] == 'P':#插入 即read多出I个碱基 过程中应忽略
#                 if encounter == False:
#                     junction += 1; encounter = True
        
#         print(data[3], data[5], data[7], data[8], bound)
        
#         if len(bound) == 4:
#             junction = (bound[1],bound[2])
#             if junction not in splice_junction:
#                 splice_junction[junction] = 1
#             else:
#                 splice_junction[junction] += 1
                
#     junction = (start,end); print(splice_junction)
#     if junction in splice_junction:
#         if splice_junction[junction] >= 2:
#             print(True)

#     cigar = '39M67106N12M25S'
#     read = []
#     temp = ''

#     for char in cigar:
#         if str.isdigit(char):
#             temp += char
#         else:
#             read.append(int(temp)); read.append(char); temp = ''
            
#     temp = bam.count_coverage('chr1', 109508246, 109508248)
#     coverage = np.zeros((4,2))

#     for row in range(len(temp)):
#         for col in range(len(temp[row])):
#             coverage[row][col]=temp[row][col]

#     bam.close()