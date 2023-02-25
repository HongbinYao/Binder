# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:35:32 2022

@author: yhb
"""
from copy import deepcopy
from adjacency_list import adjacency_list, Vertex
from read_bam import read_bam, read_coverage, out_put_bam


def decide(index, tran, bam_file, chrom, guided_exon, splice_junction):
    curr_exon = tran[index]; result = False
    for junction in splice_junction:
        if curr_exon[0] < junction[0] and junction[-1] < curr_exon[-1]:#存在junction在curr_exon内
            if read_bam(bam_file, chrom, junction[0], junction[-1]) == True:#判断bam文件中是否存在该剪接位点
            #进一步判断剪接位点内部覆盖度大小以及该exon是否在genome guided exon中
                if curr_exon not in guided_exon or read_coverage(bam_file, chrom, junction[0], junction[-1]) == True:
                    result = True; break
    return result, junction
                    
def junction_rectify(index, tran, curr_exon, temp):
    if temp[-1] < tran[index+1][-1]:#若候选剪接位点右侧在下一个exon右侧边界的左侧
        tran[index+1] = (temp[-1],tran[index+1][-1])
    elif temp[-1] >= tran[index+1][-1]:
        tran.pop(index+1); tran[index+1] = (junction[-1],tran[index+1][-1]) 
        
def deal(index, chrom, curr_exon, tran, bam_file):
    candidate = []; confirm = []; left = []; right = []
    bam_splice_junction = out_put_bam(bam_file, chrom, tran[index][-1], tran[index+1][0])
    aver_frequency = sum(list(bam_splice_junction.values()))/len(bam_splice_junction)
    if aver_frequency > 2:
        for junction, frequency in bam_splice_junction.items():
            if frequency > aver_frequency:
                candidate.append(junction)
        for junction in candidate:
            if junction[0] - tran[index][-1] == 0:#剪接位点左侧与候选剪接位点一致
                confirm.append(junction); confirm.sort()
            elif tran[index][0] < junction[0] < tran[index][-1]:#候选剪接位点在左侧
                left.append(junction); left.sort()
            elif tran[index][-1] < junction[0] < tran[index+1][0]:#候选剪接位点在右侧
                right.append(junction); right.sort()
        if confirm:
            temp = confirm[0]
            if temp[-1] < tran[index+1][-1]:#若候选剪接位点右侧在下一个exon右侧边界的左侧
                tran[index+1] = (temp[-1],tran[index+1][-1])
            elif temp[-1] >= tran[index+1][-1]:
                if index+1 < len(tran) - 1:
                    tran.pop(index+1); tran[index+1] = (junction[-1],tran[index+1][-1])
                else:
                    tran[index+1] = (junction[-1],junction[-1]+20)
        elif left:
            temp = left[0]
            tran[index] = (curr_exon[0],temp[0]); tran[index+1] = (temp[-1],tran[index+1][-1])
        elif right:
            temp = right[0]
            tran[index] = (curr_exon[0],temp[0]); tran[index+1] = (temp[-1],tran[index+1][-1])
        else:#bam文件中不存在可以支撑或修改该剪接位点的信息
            tran[index] = (curr_exon[0],tran[index+1][-1]); tran.pop(index+1)
    else:
        tran[index] = (curr_exon[0],tran[index+1][-1]); tran.pop(index+1)
    return tran

def assembly(chrom, bam_file, de_nove_temp, genome_guided_temp):
    #读取genoem guided transcript 表示成exon的形式
    genome_guided_exon = []
    for tran in genome_guided_temp:
        temp = []
        for index in range(0,len(tran)-2, 2):
            exon = ()
            exon += (tran[index],); exon += (tran[index + 1],)
            temp.append(exon)
        genome_guided_exon.append(temp)
    
    #提取剪接位点
    splice_junction = []
    for tran in genome_guided_exon:
        index = 0
        while index < len(tran) - 1:
            splice_junction.append((tran[index][-1],tran[index+1][0]))
            index += 1
    splice_junction = list(set(splice_junction)); splice_junction.sort()
    
    #提取3_prime 5_prime
    prime = []
    for tran in genome_guided_temp:
        prime.append(tran[1]); prime.append(tran[-4])
    prime = list(set(prime)); prime.sort()
    
    #提取内部剪接位点
    internal = []
    for tran in genome_guided_exon:
        for exon in tran[1:-1]:
            internal.append(exon[0]); internal.append(exon[1])
    internal = list(set(internal)); internal.sort()
    
    #提取exon
    guided_exon = []
    for tran in genome_guided_exon:
        for exon in tran:
            guided_exon.append(exon)
    guided_exon = list(set(guided_exon)); guided_exon.sort()         
    
    #读取de nove transcript 表示成exon的形式
    de_nove_exon = []
    for tran in de_nove_temp:
        temp = []
        for index in range(0,len(tran)-2, 2):
            exon = ()
            exon += (tran[index],); exon += (tran[index + 1],)
            temp.append(exon)
        de_nove_exon.append(temp)
        
    de_nove_rectify = []; pair_end = {}
    for tran in deepcopy(de_nove_exon):
        index = 0; append = False; novel_prime = False; mismatch = 0; match = 0
        #首先判断是否要将该de nove transcript的3' 5'加入到剪接图中
        if tran[0][-1] not in prime or tran[-1][0] not in prime:
            if tran[0][-1] in internal or tran[-1][0] in internal:
                novel_prime = True
        while index < len(tran) - 1:
            curr_junction = (tran[index][-1],tran[index+1][0]); curr_exon = tran[index]
            result, junction = decide(index, tran, bam_file, chrom, guided_exon, splice_junction)
            if result == True:#判断exon内部是否存在剪接位点
                curr_exon = tran[index]; mismatch += 1
                tran[index] = (curr_exon[0],junction[0]); tran.insert(index+1,(junction[-1],curr_exon[-1]))
                index += 1
            elif curr_junction in splice_junction and read_bam(bam_file, chrom, tran[index][-1], tran[index+1][0]) == True:
                index += 1; match += 1
            elif read_bam(bam_file, chrom, tran[index][-1], tran[index+1][0]) == True:
                if 1 <= index < len(tran) - 2:
                    bam_splice_junction = out_put_bam(bam_file, chrom, tran[index-1][-1], tran[index+2][0])
                    curr = bam_splice_junction[(tran[index][-1],tran[index+1][0])]
                    if (tran[index-1][-1],tran[index][0]) in bam_splice_junction:
                        left = bam_splice_junction[(tran[index-1][-1],tran[index][0])]
                    else:
                        left = 5
                    if (tran[index+1][-1],tran[index+2][0]) in bam_splice_junction:
                        right = bam_splice_junction[(tran[index+1][-1],tran[index+2][0])]
                    else:
                        right = left
                    if 2*curr < 2*left*right/(left+right) and curr < 5:#说明当前剪接位点的覆盖度远小于两侧的覆盖度 应删除
                        tran = deal(index, chrom, curr_exon, tran, bam_file)
                        # tran[index] = (curr_exon[0],tran[index+1][-1]); tran.pop(index+1)
                        index += 1; mismatch += 1
                    else:#此时应当保留
                        index += 1; match += 1
                elif index == 0:
                    bam_splice_junction = out_put_bam(bam_file, chrom, tran[index][-1], tran[index+2][0])
                    curr = bam_splice_junction[(tran[index][-1],tran[index+1][0])]
                    if (tran[index+1][-1],tran[index+2][0]) in bam_splice_junction:
                        right = bam_splice_junction[(tran[index+1][-1],tran[index+2][0])]
                    else:
                        right = 5
                    if 2*curr < right and curr < 5:
                        tran = deal(index, chrom, curr_exon, tran, bam_file)
                        index += 1; mismatch += 1
                    else:
                        index += 1; match += 1
                elif index == len(tran) - 2:
                    bam_splice_junction = out_put_bam(bam_file, chrom, tran[index-1][-1], tran[index+1][0])
                    curr = bam_splice_junction[(tran[index][-1],tran[index+1][0])]
                    if (tran[index-1][-1],tran[index][0]) in bam_splice_junction:
                        left = bam_splice_junction[(tran[index-1][-1],tran[index][0])]
                    else:
                        left = 5
                    if 2*curr < left and curr < 5:
                        tran = deal(index, chrom, curr_exon, tran, bam_file)
                        index += 1; mismatch += 1
                    else:
                        index += 1; match += 1
            else:#此时exon内部没有剪接位点且该剪接位点也并没有额外信息去支持
                mismatch += 1; candidate = []; confirm = []; left = []; right = []
                bam_splice_junction = out_put_bam(bam_file, chrom, tran[index][-1], tran[index+1][0])
                if bam_splice_junction:
                    aver_frequency = sum(list(bam_splice_junction.values()))/len(bam_splice_junction)
                else:
                    aver_frequency = 0
                if aver_frequency > 2:
                    for junction, frequency in bam_splice_junction.items():
                        if frequency > aver_frequency:
                            candidate.append(junction)
                    for junction in candidate:
                        if junction[0] - tran[index][-1] == 0:#剪接位点左侧与候选剪接位点一致
                            confirm.append(junction); confirm.sort()
                        elif tran[index][0] < junction[0] < tran[index][-1]:#候选剪接位点在左侧
                            left.append(junction); left.sort()
                        elif tran[index][-1] < junction[0] < tran[index+1][0]:#候选剪接位点在右侧
                            right.append(junction); right.sort()
                    if confirm:
                        temp = confirm[0]
                        if temp[-1] < tran[index+1][-1]:#若候选剪接位点右侧在下一个exon右侧边界的左侧
                            tran[index+1] = (temp[-1],tran[index+1][-1])
                        elif temp[-1] >= tran[index+1][-1]:
                            if index+1 < len(tran) - 1:
                                tran.pop(index+1); tran[index+1] = (junction[-1],tran[index+1][-1])
                            else:
                                tran[index+1] = (junction[-1],junction[-1]+20)
                    elif left:
                        temp = left[0]
                        tran[index] = (curr_exon[0],temp[0]); tran[index+1] = (temp[-1],tran[index+1][-1])
                    elif right:
                        temp = right[0]
                        tran[index] = (curr_exon[0],temp[0]); tran[index+1] = (temp[-1],tran[index+1][-1])
                    else:#bam文件中不存在可以支撑或修改该剪接位点的信息
                        tran[index] = (curr_exon[0],tran[index+1][-1]); tran.pop(index+1)
                else:
                    tran[index] = (curr_exon[0],tran[index+1][-1]); tran.pop(index+1)
                index += 1
                            
        if match > mismatch:
            append = True
        
        if append == True:
            if novel_prime == True:
                de_nove_rectify.append(tran)
                pair_end[tran[0]] = tran[-1]
            elif novel_prime == False:
                de_nove_rectify.append(tran)
    return de_nove_rectify, pair_end
   
def assembly_sr(chrom, bam_file, de_nove_temp):
    #读取de nove transcript 表示成exon的形式
    de_nove_exon = []
    for tran in de_nove_temp:
        temp = []
        for index in range(0,len(tran)-2, 2):
            exon = ()
            exon += (tran[index],); exon += (tran[index + 1],)
            temp.append(exon)
        de_nove_exon.append(temp)
        
    de_nove_rectify = []; pair_end = {}
    for tran in deepcopy(de_nove_exon):
        index = 0; append = False; novel_prime = True; mismatch = 0; match = 0
        #将de nove transcript的3' 5'加入到剪接图中
        while index < len(tran) - 1:
            curr_junction = (tran[index][-1],tran[index+1][0]); curr_exon = tran[index]
            if read_bam(bam_file, chrom, tran[index][-1], tran[index+1][0]) == True:
                index += 1; match += 1
            else:#此时exon内部没有剪接位点且该剪接位点也并没有额外信息去支持
                mismatch += 1; candidate = []; confirm = []; left = []; right = []
                bam_splice_junction = out_put_bam(bam_file, chrom, tran[index][-1], tran[index+1][0])
                if bam_splice_junction:
                    aver_frequency = sum(list(bam_splice_junction.values()))/len(bam_splice_junction)
                else:
                    aver_frequency = 0
                if aver_frequency >= 1:
                    for junction, frequency in bam_splice_junction.items():
                        if frequency > aver_frequency:
                            candidate.append(junction)
                    for junction in candidate:
                        if junction[0] - tran[index][-1] == 0:#剪接位点左侧与候选剪接位点一致
                            confirm.append(junction); confirm.sort()
                        elif tran[index][0] < junction[0] < tran[index][-1]:#候选剪接位点在左侧
                            left.append(junction); left.sort()
                        elif tran[index][-1] < junction[0] < tran[index+1][0]:#候选剪接位点在右侧
                            right.append(junction); right.sort()
                    if confirm:
                        temp = confirm[0]
                        if temp[-1] < tran[index+1][-1]:#若候选剪接位点右侧在下一个exon右侧边界的左侧
                            tran[index+1] = (temp[-1],tran[index+1][-1])
                        elif temp[-1] >= tran[index+1][-1]:
                            if index+1 < len(tran) - 1:
                                tran.pop(index+1); tran[index+1] = (junction[-1],tran[index+1][-1])
                            else:
                                tran[index+1] = (junction[-1],junction[-1]+20)
                    elif left:
                        temp = left[0]
                        tran[index] = (curr_exon[0],temp[0]); tran[index+1] = (temp[-1],tran[index+1][-1])
                    elif right:
                        temp = right[0]
                        tran[index] = (curr_exon[0],temp[0]); tran[index+1] = (temp[-1],tran[index+1][-1])
                    else:#bam文件中不存在可以支撑或修改该剪接位点的信息
                        tran[index] = (curr_exon[0],tran[index+1][-1]); tran.pop(index+1)
                else:
                    tran[index] = (curr_exon[0],tran[index+1][-1]); tran.pop(index+1)
                index += 1
                            
        if match > mismatch:
            append = True
        
        if append == True:
            if novel_prime == True:
                de_nove_rectify.append(tran)
                pair_end[tran[0]] = tran[-1]
            elif novel_prime == False:
                de_nove_rectify.append(tran)
    return de_nove_rectify, pair_end

def test(de_nove_rectify):
    out_put_file = open(path + r'Binder.gtf','w')
    tran_index = 0
    for tran in de_nove_rectify:
        tran_index += 1
        out_put_file.write(chrom + '\t' +\
                          'Binder' + '\t' +\
                          'transcript' + '\t' +\
                          str(tran[0][0]) + '\t' +\
                          str(tran[-1][-1]) + '\t' +\
                          '1000' + '\t' +\
                          '-' + '\t' +\
                          '.' + '\t' +\
                          'gene_id ' + '"' + chrom + '.' + str(tran_index) + '"; ' +\
                          'transcript_id ' + '"' + chrom + '.' + str(tran_index) + '.' + str(tran_index) + '";' + '\n')
        exon_index = 0
        for index in range(0,len(tran)):
            exon_index += 1
            out_put_file.write(chrom + '\t' +\
                              'Binder' + '\t' +\
                              'exon' + '\t' +\
                              str(tran[index][0]) + '\t' +\
                              str(tran[index][-1]) + '\t' +\
                              '1000' + '\t' +\
                              '-' + '\t' +\
                              '.' + '\t' +\
                              'gene_id ' + '"' + chrom + '.' + str(tran_index) + '"; ' +\
                              'transcript_id ' + '"' + chrom + '.' + str(tran_index) + '.' + str(tran_index) + '";' +\
                              'exon_number ' + '"' + str(exon_index) + '";' + '\n')
    out_put_file.close()
    
# de_nove_rectify, pair_end = assembly(chrom, bam_file, de_nove_temp, genome_guided_temp); test(de_nove_rectify)
