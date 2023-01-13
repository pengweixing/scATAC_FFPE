#################################################
#  File Name:MarkOptDup.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Tue Nov 22 13:10:20 2022
#################################################

import sys
import pysam
import argparse
import copy
from orderedset import OrderedSet
import random
import math

def fargv():
    parser = argparse.ArgumentParser(usage="python ")
    parser.add_argument('-i',"--input",help="the sorted bam file ", required=True)
    parser.add_argument('-o',"--output",help="the sorted bam file ", required=True)
    parser.add_argument('-t',"--type",help="single or pair ", required=True)
    parser.add_argument('-n',"--Tn5name",help="the name of each Tn5 for same sample ", required=True)
    args = parser.parse_args()
    return args

class BamProcess:
    def __init__(self, bam_input = '',bam_output = '',name=''):
        self.bam_input = bam_input
        self.bam_output = bam_output
        self.cache = {'cellname':'','start':'','insert_len':'','R2_start':'','readsname':''}
        self.name = name
        
    def Classic_Opt_Mark_Dup_single(self):
        """for the reads within each cell, check if they have same [cellname, start ,end],similar with picardtools"""
        def assign():
            self.bam_output.write(a)
            self.cache['cellname'] = cell_name
            self.cache['start'] = mystart
    
        for a in self.bam_input:
            if not a.is_unmapped:
                query_name = a.query_name
                out_query_name = query_name+"_"+self.name
                a.query_name = out_query_name
                cell_name = query_name.split('_')[0]

            #    if a.is_reverse:
            #        mystart = a.reference_start + a.query_alignment_end
             #   else:
                mystart = a.reference_start
                
                if mystart == self.cache['start']:
                    
                    if cell_name != self.cache['cellname']:
                        assign()
                    else:
                        if a.flag < 1000 or a.flag > 2000:
                            a.flag = a.flag + 1024
                        assign()
                else:
                    assign()
            else:
                self.bam_output.write(a)
            
            
    def Harsh_Opt_Mark_Dup(self):
        """mark it as dup if R2 start position is same"""   
        
        def assign():
            self.cache['cellname'] = cell_name
            self.cache['start'] = a.reference_start
            self.cache['insert_len'] = a.template_length
            self.cache['R2_end'] = read2_end
            self.cache['readsname'] = a.query_name
            
            if read2_end in temp_corr_dict:
                temp_corr_dict[read2_end].add(cell_name)
            else:
                temp_corr_dict[read2_end] = set([cell_name])
        read2_start = ''
        read2_end = ''
        self.cache['R2_end'] = read2_end
        cell_linear_amp = {}
        cell_opt_amp = {}
        j = 1
        temp_corr_dict = {}
        temp_read1_start = 0
        for a in self.bam_input: 
            if j%100 == 0:
                for each in temp_corr_dict:
                    if each < read2_end - 3: 
                        del temp_corr_dict[each]
                        
            if not a.is_unmapped:
                query_name = a.query_name
                cell_name = query_name.split('_')[0]
                
                if a.is_read2:
                    if a.is_reverse:
                        read2_start = a.reference_start
                        read2_end = read2_start + a.query_alignment_end
                        self.cache['R2_end'] = read2_end
                    else:
                        read2_end = a.reference_start
                        self.cache['R2_end'] = read2_end
                    
                if read2_end == self.cache['R2_end']:
                    if read2_end in temp_corr_dict:
                        if query_name != self.cache['readsname'] and cell_name in temp_corr_dict[read2_end]:
                            if cell_name in cell_linear_amp:
                                cell_linear_amp[cell_name].add(query_name)
                            else:
                                cell_linear_amp[cell_name] = OrderedSet([query_name])
                            assign()
                        else:
                            assign()
                    else:
                        assign()
                else:
                    assign()

        self.bam_input.reset()

        for a in self.bam_input:
            query_name = a.query_name
            a.query_name = query_name+"_"+self.name
            cell_name = query_name.split('_')[0]
            if not a.is_unmapped:
                if cell_name in cell_linear_amp:
                    dup_name = cell_linear_amp[cell_name]
                    if query_name in dup_name:
                        if query_name != dup_name[-1]:
                            if a.flag < 1000 or a.flag > 2000:
                                a.flag = a.flag + 1024
                                self.bam_output.write(a)
                            else:
                                self.bam_output.write(a)
                        else:
                            self.bam_output.write(a)
                    else:
                        self.bam_output.write(a)
                else:
                    self.bam_output.write(a)
            else:
                self.bam_output.write(a)
    
def main(kwargs):
    f_bam_input = pysam.AlignmentFile(kwargs.input,mode='rb',check_header=True,threads=3)
    f_bam_output = pysam.AlignmentFile(kwargs.output,mode='wb',template=f_bam_input,threads=3)
    bamobject = BamProcess(bam_input = f_bam_input,bam_output = f_bam_output,name=kwargs.Tn5name)
  #  
    if kwargs.type == "single":
        bamobject.Classic_Opt_Mark_Dup_single()
    elif kwargs.type == "pair":
        bamobject.Harsh_Opt_Mark_Dup()
    else:
        bamobject.Harsh_Opt_Mark_Dup()

if __name__ == "__main__":
    kwargs = fargv()
    main(kwargs)
