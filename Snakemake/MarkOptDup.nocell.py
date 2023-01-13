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
import numpy

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
        self.cache = {'cellname':'','start':'','insert_len':0,'R2_start':0,'readsname':''}
        self.name = name

    """### this function has not been compeleted, don't use it!!"""
    def Classic_Opt_Mark_Dup_single(self):  
        pass            
    def Harsh_Opt_Mark_Dup(self):
        """mark it as dup if R2 start position is same"""   
        
        def assign():
            self.cache['cellname'] = cell_name
            self.cache['R2_start'] = R2_start
            self.cache['readsname'] = a.query_name
            
            if R2_start in temp_coor_dict:
                temp_coor_dict[R2_start].add(cell_name)
            else:
                temp_coor_dict[R2_start] = set([cell_name])
        
        ### Traversing the forward R2
        print("Do forward duplicates marking ")
        cell_linear_amp = {}
       # j = 1
        temp_coor_dict = {}  ## save the cell names for the reads which have the same coordinates
        i = 1
        max_len = 0
        max_reads_name = 'na'
        max_cell_name = 'na'

        for a in self.bam_input: 
            ## clear temp_coor_dict
            if i % 1000000==0:
                print("%s Million reads has been done" % str(i/1000000))
            if len(temp_coor_dict)>5:
                aa = list(temp_coor_dict.keys())
                for each in aa:
                    if each < R2_start: 
                        del temp_coor_dict[each]
            i += 1 
            if not a.is_unmapped:
                query_name = a.query_name
                cell_name = query_name.split('_')[0]
                if a.is_read2:
                    if not a.is_reverse:
                        R2_start = a.reference_start
                        if R2_start == self.cache['R2_start']:
                            ## check if they are belongs to same cell
                      #      if cell_name in temp_coor_dict[R2_start]:  
                            if a.template_length > max_len: ## if dna length > last one, the last one is labeled with duplicates      
                                if max_cell_name in cell_linear_amp:
                                    cell_linear_amp[max_cell_name].add(max_reads_name)
                                else:
                                    cell_linear_amp[max_cell_name] = OrderedSet([max_reads_name])
                                
                                max_cell_name = cell_name
                                max_reads_name = query_name
                                max_len = int(a.template_length)
                            else: ## if dna length < last one, the current one is labeled with duplicates
                                if cell_name in cell_linear_amp:
                                    cell_linear_amp[cell_name].add(query_name)
                                else:
                                    cell_linear_amp[cell_name] = OrderedSet([query_name])
                                
                        #        assign()
                        #    else:
                            assign()
                        else:
                            assign()
                            max_len = 0
                            max_reads_name = 'na'
                            max_cell_name = 'na'
                    

        self.bam_input.reset()
        ### Traversing the reverse Reads2
        print("Do reverse duplicates marking ")
        temp_coor_dict = {}  ## save the cell names for the reads which have the same coordinates
        max_len = 0
        max_reads_name = 'na'
        max_cell_name = 'na'
        iteri = 1
        for a in self.bam_input:
            if iteri % 1000000==0:
                print("%s Million reads has been done" % str(iteri/1000000))
            iteri += 1    
            if not a.is_unmapped:
                query_name = a.query_name    
                if a.is_read2:
                    if a.is_reverse:     
                        R2_start = a.reference_end
                        frag_len = abs(a.template_length)
                        if not temp_coor_dict:
                            if R2_start in temp_coor_dict:
                                temp_coor_dict[R2_start].add(query_name+"___"+str(frag_len))
                            else:
                                temp_coor_dict[R2_start] = set([query_name+"___"+str(frag_len)])
                        elif a.reference_start < min(temp_coor_dict):

                            if R2_start in temp_coor_dict:

                                temp_coor_dict[R2_start].add(query_name+"___"+str(frag_len))
                            else:
                                temp_coor_dict[R2_start] = set([query_name+"___"+str(frag_len)])
                        else:
                            all_value = list(temp_coor_dict[min(temp_coor_dict)])
              
                            if len(all_value)>1:  # all_reads_name and length
                                all_len = [each_reads.split("___")[-1] for each_reads in all_value]
                                max_index = numpy.argmax(all_len)
                                for ii in range(len(all_len)):
                                    if ii != max_index:
                                        readsname = all_value[ii].split("___")[0]
                                        cell_name = readsname.split("_")[0]
                                        if cell_name in cell_linear_amp:
                                            cell_linear_amp[cell_name].add(readsname)
                                        else:
                                            cell_linear_amp[cell_name] = OrderedSet([readsname])
                                del temp_coor_dict[min(temp_coor_dict)]   
                            else:
                                del temp_coor_dict[min(temp_coor_dict)]

        ## Traversing the end of the bam file
        temp_cell = {}
        keys = list(temp_coor_dict.keys())
        for bb in keys:
            all_value = list(temp_coor_dict[bb])
            if len(all_value)>1:
                all_len = [each_reads.split("___")[1] for each_reads in all_value]
                max_index = numpy.argmax(all_len)
                for ii in range(len(all_len)):
                    if ii != max_index:
                        readsname = all_value[ii].split("___")[0]
                        cell_name = readsname.split("_")[0]
                        if cell_name in cell_linear_amp:
                            cell_linear_amp[cell_name].add(readsname)
                        else:
                            cell_linear_amp[cell_name] = OrderedSet([readsname])
               
                del temp_coor_dict[bb]   
            else:
                del temp_coor_dict[bb]
        

        ### do duplicates marking 
        self.bam_input.reset()
        i = 0
        print("Do writing")
        for a in self.bam_input:
            if i % 1000000==0:
                print("%s Million reads has been written" % str(i/1000000))
            i += 1
            query_name = a.query_name
            a.query_name = query_name+"_"+self.name
            cell_name = query_name.split('_')[0]
            if not a.is_unmapped:
                if cell_name in cell_linear_amp:
                    dup_name = cell_linear_amp[cell_name]
                    if query_name in dup_name:
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
