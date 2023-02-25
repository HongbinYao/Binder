# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 22:15:04 2022

@author: yhb
"""
def write_gtf(assembly_trans, out_put_file):
    for chrom, trans in assembly_trans.items():
        for tran, exons in trans.items():
            out_put_file.write(chrom[0:-1] + '\t' +\
                              'Binder' + '\t' +\
                              'transcript' + '\t' +\
                              str(exons[0]) + '\t' +\
                              str(exons[-1]) + '\t' +\
                              '1000' + '\t' +\
                              chrom[-1] + '\t' +\
                              '.' + '\t' +\
                              'gene_id ' + '"' + tran[0:-2] + '"; ' +\
                              'transcript_id ' + '"' + tran + '";' + '\n')
            exon_index = 0
            for index in range(0,len(exons),2):
                exon_index += 1
                out_put_file.write(chrom[0:-1] + '\t' +\
                                  'Binder' + '\t' +\
                                  'exon' + '\t' +\
                                  str(exons[index]) + '\t' +\
                                  str(exons[index + 1]) + '\t' +\
                                  '1000' + '\t' +\
                                  chrom[-1] + '\t' +\
                                  '.' + '\t' +\
                                  'gene_id ' + '"' + tran[0:-2] + '"; ' +\
                                  'transcript_id ' + '"' + tran + '"; ' +\
                                  'exon_number ' + '"' + str(exon_index) + '";' + '\n')
    out_put_file.close()
    return None

def write_merge_ends(chrom, chain, output_gtf, merge_ends):
    for key, value in merge_ends.items():
        for exon in value:
            output_gtf.write('chr' + str(chrom) + '\t' +\
                              chain + '\t' +\
                              'start' + '\t' + str(key[0]) + '\t' +\
                              'end' + '\t' + str(exon[1]) + '\t')
            output_gtf.write('\n')
    
    return None

def write_prime(assembly_trans, output_gtf):
    for chrom, trans in assembly_trans.items():
        for tran_id, tran in trans.items():
            output_gtf.write(chrom[0:-1] + '\t' +\
                              'Binder' + '\t' +\
                              'transcript' + '\t' +\
                              str(tran[0]) + '\t' +\
                              str(tran[-1]) + '\t' +\
                              '1000' + '\t' +\
                              chrom[-1] + '\t' +\
                              '.' + '\t' +\
                              'gene_id ' + '"' + tran_id + '"; ' +\
                              'transcript_id ' + '"' + tran_id + '";' + '\n')
            exon_index = 0
            for index in range(0,len(tran),2):
                exon_index += 1
                output_gtf.write(chrom[0:-1] + '\t' +\
                                  'Binder' + '\t' +\
                                  'exon' + '\t' +\
                                  str(tran[index]) + '\t' +\
                                  str(tran[index + 1]) + '\t' +\
                                  '1000' + '\t' +\
                                  chrom[-1] + '\t' +\
                                  '.' + '\t' +\
                                  'gene_id ' + '"' + tran_id + '"; ' +\
                                  'transcript_id ' + '"' + tran_id + '"; ' +\
                                  'exon_number ' + '"' + str(exon_index) + '";')
                output_gtf.write('\n')
    return None

def write_merge(merge_write, out_put_file):
    last_chrom = next(iter(merge_write))[0:-1]; gene_index = 0
    for chrom, trans in merge_write.items():
        curr_chrom = chrom[0:-1]
        if curr_chrom == last_chrom:
            tran_index = 0
        else:#curr_chrom != last_chrom
            gene_index = 0; last_chrom = curr_chrom
            tran_index = 0
        for region in trans:
            gene_index += 1
            for tran in region:
                if not tran:
                    continue
                tran_index += 1
                
                out_put_file.write(chrom[0:-1] + '\t' +\
                                  'Binder' + '\t' +\
                                  'transcript' + '\t' +\
                                  str(tran[0]) + '\t' +\
                                  str(tran[-1]) + '\t' +\
                                  '1000' + '\t' +\
                                  chrom[-1] + '\t' +\
                                  '.' + '\t' +\
                                  'gene_id ' + '"' + chrom[0:-1] + '.' + str(gene_index) + '"; ' +\
                                  'transcript_id ' + '"' + chrom[0:-1] + '.' + str(gene_index) + '.' + str(tran_index) + '";' + '\n')
                exon_index = 0
                for index in range(0,len(tran),2):
                    exon_index += 1
                    out_put_file.write(chrom[0:-1] + '\t' +\
                                      'Binder' + '\t' +\
                                      'exon' + '\t' +\
                                      str(tran[index]) + '\t' +\
                                      str(tran[index + 1]) + '\t' +\
                                      '1000' + '\t' +\
                                      chrom[-1] + '\t' +\
                                      '.' + '\t' +\
                                      'gene_id ' + '"' + chrom[0:-1] + '.' + str(gene_index) + '"; ' +\
                                      'transcript_id ' + '"' + chrom[0:-1] + '.' + str(gene_index) + '.' + str(tran_index) + '"; ' +\
                                      'exon_number ' + '"' + str(exon_index) + '";' + '\n')
    out_put_file.close()
    return None


if __name__ == '__main__':
    path = r'C:\Users\yhb\Desktop\gtf文件\gtf\\'
    merge_file = open(path + 'new.gtf','w')
    a=[]
    a.append(de_nove_handle)
    write_rectify(14, "+", 0, merge_file, a)
    merge_file.close()
    write_merge_ends(chrom, chain, merge_file, merge_ends)
