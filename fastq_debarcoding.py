#################################################
#  File Name:fastq_to_10x.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com, pengweixing@igp.uu.se
#  Created Time: Thu Nov 25 17:31:28 2021
#################################################
HELP = """ Example:
---------------------------------------------------------------------------
R2 sequence with four structures: 5' -> 3'                                                                      <--- R1         
cell_bc1              cell_bc2             cell_bc3    linker   sample_index                   
ACGATTGNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNNNNNNNNNNNNNNNN   -> ME -> genomicDNA -> ME-rev-com -> 
ACGATTGNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNNNNNNNNNNNNNNNN   -> ME -> genomicDNA -> ME-rev-com ->
ACGATTGNNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNNNNNNNNNNNNNNNN   -> ME -> genomicDNA -> ME-rev-com ->
ACGATTGNNNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNNNNNNNNNNNNNNNN   -> ME -> genomicDNA -> ME-rev-com ->
---------------------------------------------------------------------------
linker ATCCACGAGCATTCG
---------------------------------------------------------------------------
Trimming:        |        |
R1 5' NNNNN->ME->genomicDNA->ME-rev_comp->NNNN
R2 5' NNNNN->ME->genomicDNA->ME-rev_comp->NNNN
                 |        |          

ME: AGATGTGTATAAGAGACAG
ME_rev_comp: CTGTCTCTTATACACATCT

---------------------------------------------------------------------------
Barcode.txt:
---------------------------------------------------------------------------
[Sample_index]
Brain   ATC
cere    TGA
...

[Cell_BC1]
R02_#01 AAACCGG
R02_#02 AAACGTC
R02_#03 AAAGATG
...

[Cell_BC2]
R03_#01 AAACCGG
R03_#02 AAACGTC
R03_#03 AAAGATG

[Cell_BC3]
R04_#01 AAACCGG
R04_#02 AAACGTC
R04_#03 AAAGATG
...
---------------------------------------------------------------------------
---------------------------------------------------------------------------

"""

import sys
import gzip
import argparse
import collections
from collections import Counter
import numpy as np
import os
import pandas as pd
import math
import Levenshtein
from fuzzywuzzy import fuzz
from matplotlib import pyplot as plt
from multiprocessing import Pool,Process
import threading
from itertools import islice

def fargv():
    parser = argparse.ArgumentParser(usage="python debarcode.py -r1 R1.fastq.gz -r2 R2.fastq.gz -b barcode.list -o output_name")
    parser.add_argument('-r1',"--R1",help="the R1 of file ", required=True)
    parser.add_argument('-r2',"--R2",help="the R2 of file ", required=True)
    parser.add_argument('-b',"--barcode",help="the barcode list ", required=True)
    parser.add_argument('-o',"--out_dir",help="the directory for output", default='./')
    parser.add_argument('-a',"--adaptor",help="the adaptor's sequence",type=str,default='GGAGAAGATGTGTATAAGAGACAG')
    parser.add_argument('-p',"--processor",help="the number of processors",default=10,type=int)
    parser.add_argument('-ml',"--min_length_seq",help="the minimum length of R1 and R2 after trimming ",default=19,type=int)
    
    args = parser.parse_args()
    return args

def process_R1_R2_wrapper(Fastq_transform,data_R1_each_chunk,data_R2_each_chunk):
    return Fastq_transform.process_R1_R2(data_R1_each_chunk,data_R2_each_chunk)

class Fastq_transform:
    names = locals()
    def __init__(self,R1_input = [], R2_input = [],output_dir = './', min_frags_cutoff_plot = 2,\
     bins=100, adaptor = '',min_length_seq = 21,cpu = 10):

        self.R1_input = R1_input
        self.R2_input = R2_input
        self.pos = [(0,7),(22,29),(44,51)]
        """linker sequence between sample index and cell BC1"""
        self.linker = "ATCCACGAGCATTCG"
        self.me = "AGATGTGTATAAGAGACAG"
        self.sample_index_map,self.bc1,self.bc2,self.bc3 = {},[],[],[]
        self.output_dir = output_dir
        self.adaptor = adaptor
        self.threhold_for_Levenshtein = 80
        self.adaptor_set = []
        self.adaptor_rev_comp_set = []
        self.mismatch = 3
        self.min_length_seq = min_length_seq
        self.cutoffR1 = 75
        self.cell_bc_total_len = 91
        self.mychunksize = 5000000  #### read 5M lines per time and distribute them to multi cpus
        self.Processor = cpu
        self.all_sc_barcodes = []
        self.reads_number = 1
        self.cutoffR2 = 59 ## keep 59bp of R2 from 3'
        self.length_for_complete = 15

    def init_adaptors(self):
        adaptor_comp = self.DNA_complement(self.adaptor)
        adaptor_rev_comp = self.DNA_reverse(adaptor_comp)
        self.adaptor_set = self.__gen_adaptor(self.adaptor,rev=False)
        adaptor_rev_comp_set = self.__gen_adaptor(adaptor_rev_comp,rev=True)
        self.adaptor_rev_comp_set = adaptor_rev_comp_set[::-1]

    def mkdir(self): 
        if not self.output_dir == "./":
            try:
                os.mkdir(self.output_dir)
            except OSError as error:
                print("Warning: The directory of %s has already exists\n" % self.output_dir)

    def DNA_complement(self,sequence):
        sequence = sequence.upper()
        sequence = sequence.replace('A', 't')
        sequence = sequence.replace('T', 'a')
        sequence = sequence.replace('C', 'g')
        sequence = sequence.replace('G', 'c')
        return sequence.upper()

    def DNA_reverse(self,sequence):
        sequence = sequence.upper()  
        return sequence[::-1]
    def print_error(self,value):
        print(value)

    def distribute_to_processor(self):
        
        data_R1 = gzip.open(self.R1_input,'rt')
        data_R2 = gzip.open(self.R2_input,'rt')
        i = 1
        Barcode_reads_total = {}
        cell_number_each_total = {}
        Total_reads = 0
        Barcode_reads_total_C =  Counter(Barcode_reads_total)
        cell_number_each_total_C =  Counter(cell_number_each_total)
        data_R1_chunk = []
        while 1:
            if not data_R1_chunk:
                data_R1_chunk = list(islice(data_R1, 0, self.mychunksize, None))
                data_R2_chunk = list(islice(data_R2, 0, self.mychunksize, None))
     
            if data_R1_chunk and len(data_R1_chunk)/4 > self.Processor:
                pool = Pool(self.Processor)
                reads = self.mychunksize/4        ## The total number of reads
                block = int(reads/self.Processor) ## The number of reads for each processor
                results = []
                ## split the total reads into different part and distribute them to each cpu
                for index in range(self.Processor): 
                    start = index*block*4
                    end = (index+1)*block*4
                    data_R1_each_chunk = data_R1_chunk[start:end]
                    data_R2_each_chunk = data_R2_chunk[start:end]
                    results.append(pool.apply_async(process_R1_R2_wrapper,args = (self,data_R1_each_chunk,data_R2_each_chunk,)\
                    ,error_callback=self.print_error))
                """pre-read the data for next loop while the process is running """
                data_R1_chunk = list(islice(data_R1, 0, self.mychunksize, None))
                data_R2_chunk = list(islice(data_R2, 0, self.mychunksize, None))
                                
                pool.close()
                pool.join()
                for result in results:
                    result = result.get()
                    """
                    result[0]: R1
                    result[1]: cell barcodes, R2
                    result[2]: R3 
                    """
                    #### write processed reads to different files with different threads
                    if not result[0]:
                        print("Warning: No reads found in this decoding loop")
                    x1 = threading.Thread(target=self.__write_to_file_R1_R3, args=(result[0],result[1],'R1',))# R1 barcode
                    x3 = threading.Thread(target=self.__write_to_file_R1_R3, args=(result[2],result[1],'R3',))# R2 barcode
                    x1.start(),x3.start()
                    x1.join(),x3.join()      
                    ### merge all dicts of cell * fragments and reads number * sample 
                    Barcode_reads_each_C = Counter(result[3])
                    cell_number_each_C = Counter(result[4])
                    Total_reads_each = result[5]
                    Total_reads = Total_reads + Total_reads_each
                    Barcode_reads_total_C = Counter(dict(Barcode_reads_total_C+Barcode_reads_each_C))
                    cell_number_each_total_C = Counter(dict(cell_number_each_total_C+cell_number_each_C))
                print("%s reads have been done" % str(i*self.mychunksize/4))
                i += 1
            else:
                break
        all_stat = self.__write_stat(Barcode_reads_total_C,cell_number_each_total_C,Total_reads)
   #     self.hist_plot(all_stat)
        
    def __write_stat(self,Barcode_reads_total_C,cell_number_each_total_C,Total_reads):

        """write out the fragments number * cells to file"""
        Barcode_reads_total = dict(Barcode_reads_total_C)
        cell_number_each_total = dict(cell_number_each_total_C)
        f = open(self.output_dir+"/Barcodes_fragments_qc.txt",'w')
        f.write('sample\tBC1\tBC2\tBC3\tFrags_number\n')
        all_stat = []
        for key,value in cell_number_each_total.items():
            all_stat.append([self.sample_index_map[key[0]],key[1],key[2],key[3],value])
            f.write("%s\t%s\t%s\t%s\t%s\n" % (self.sample_index_map[key[0]],key[1],key[2],key[3],value))
        f.close()

        """write the statistics of debarcoding rate to file """
        f = open(self.output_dir+"/Barcoding_rate_qc.txt",'w')
        f.write('Total_reads\t%s\n' % Total_reads)
        f.write('sample\tFinal_reads\n')
        for key in Barcode_reads_total.keys():
            barcoded_reads = Barcode_reads_total[key]
            f.write('%s\t%s\n' % (key,barcoded_reads))
        f.close()
        return all_stat
        
    def process_R1_R2(self,data_R1_each_chunk,data_R2_each_chunk):
        
        buffer_R1 = []
        buffer_R3 = []
        buffer_R2 = []
        Barcode_reads = {}
        cell_number_order_dict = {}
        f_R1 = iter(data_R1_each_chunk)
        f_R2 = iter(data_R2_each_chunk)
        Total_reads = 0

        for R1_field,R2_field in zip(self.__readfq(f_R1),self.__readfq(f_R2)):
            R1_name,R1_seq,R1_qual = self._trimming(R1_field,reads ='R1')
            R2_name,R2_seq,R2_qual = self._trimming(R2_field,reads = 'R2')
            sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3 = self.__R2_Demultiplexing(R2_field)    
            Total_reads += 1 
                

            if sample_index in self.sample_index_map and \
            (temp_cell_bc1,temp_cell_bc2,temp_cell_bc3) in self.all_sc_barcodes:
                ### count the reads which have correct barcodes
                if sample_index in self.sample_index_map.keys():
                    if self.sample_index_map[sample_index] in Barcode_reads.keys():  
                        Barcode_reads[self.sample_index_map[sample_index]] += 1
                    else:
                        Barcode_reads[self.sample_index_map[sample_index]] = 1
                #  count the fragments number for each cell
                if (sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3) in cell_number_order_dict.keys():
                    cell_number_order_dict[(sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3)] += 1
                else:
                    cell_number_order_dict[(sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3)] = 1
                    
                if len(R1_seq) >= self.min_length_seq and len(R2_seq) >= self.min_length_seq:
                    ### keep the trimed and decoded reads in the list
                    buffer_R1.append((R1_name+'\n',R1_seq,R1_qual))
                    buffer_R2.append((sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3))
                    buffer_R3.append((R2_name+'\n',R2_seq,R2_qual))

        return buffer_R1,buffer_R2,buffer_R3,Barcode_reads,cell_number_order_dict,Total_reads
       
    def init_write(self):
        """initialize the output file"""
        for outname in self.sample_index_map.values():
        #    Fastq_transform.names['fI'+outname] = gzip.open(self.output_dir+"/"+outname+'_S1_L001_I1_001.fastq.gz','wb')  ### sample index
            Fastq_transform.names['f1'+outname] = gzip.open(self.output_dir+"/"+outname+'_S1_L001_R1_001.fastq.gz','wb')   ### Read 1 of genomic DNA
            Fastq_transform.names['f3'+outname] = gzip.open(self.output_dir+"/"+outname+'_S1_L001_R2_001.fastq.gz','wb')   ### Read 2 of genomic DNA
            
    def close_file(self):

        for outname in self.sample_index_map.values():
         #   Fastq_transform.names['fI'+outname].close()
            Fastq_transform.names['f1'+outname].close()
            Fastq_transform.names['f3'+outname].close()

            print("the %s has been closed\n" % outname)

    def __write_to_file_R1_R3(self,buffer_R1_R3,buffer_R2,which_reads):
        for R1_record,R2_record in zip(buffer_R1_R3,buffer_R2):
            R1_name,R1_seq,R1_qual = R1_record
            sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3 = R2_record
            outname = ''

      #      if sample_index in self.sample_index_map:  ## get the sample name from each sample index
            outname = self.sample_index_map[sample_index] 
       #     if R2_record[1:4] in self.all_sc_map_10x:   ## get the 10x barcode from each triple sc barcode
            barcode_name = R2_record[1]+R2_record[2]+R2_record[3]
            R1_name = R1_name.strip().split()[0]
            read_header_1 = "@" + barcode_name+"_" + '_'.join(R1_name.split(':')[3:]) + " 1:N:0"
            read_header_2 = "@" + barcode_name+"_" + '_'.join(R1_name.split(':')[3:]) + " 2:N:0"
      #      read_header_2 = "@" + barcode_name+"_"+n.zfill(str(self.reads_number))+" 2:N:0"
            temp_for_write_R1 = {}
            temp_for_write_R2 = {}
            
            if which_reads == "R1":
                if outname not in temp_for_write_R1.keys():
                    temp_for_write_R1[outname]  = read_header_1 +"\n"+ R1_seq + '\n' + '+\n' + R1_qual + '\n'
                else:
                    temp_for_write_R1[outname]  = temp_for_write_R1[outname] + read_header_1 + "\n" + R1_seq + '\n' + '+\n' + R1_qual + '\n'
            
            elif which_reads == "R3":
                if outname not in temp_for_write_R1.keys():
                    temp_for_write_R1[outname]  = read_header_1 +"\n"+ R1_seq + '\n' + '+\n' + R1_qual + '\n'
                else:
                    temp_for_write_R1[outname]  = temp_for_write_R1[outname] + read_header_2 + "\n" + R1_seq + '\n' + '+\n' + R1_qual + '\n'
    
            for outname in temp_for_write_R1.keys():        
                if which_reads == "R1":
                    Fastq_transform.names['f1'+outname].write(temp_for_write_R1[outname].encode())
                    Fastq_transform.names['f1'+outname].flush()
                elif which_reads == "R3":
                    Fastq_transform.names['f3'+outname].write(temp_for_write_R1[outname].encode())
                    Fastq_transform.names['f3'+outname].flush()
                
    ## split the adaptor into many sub-sequences, because the adaptor in the reads might not be intact            
    def __gen_adaptor(self,adaptor,rev=False):
        adaptor_subset = []
        for i in range(len(adaptor)):
            if not rev:
                if len(adaptor[i:len(adaptor)]) > 8:
                    adaptor_subset.append(adaptor[i:len(adaptor)])
            else:
                if len(adaptor[0:i+1]) > 8:
                    adaptor_subset.append(adaptor[0:i+1])
        return adaptor_subset

    """modified from https://github.com/TheJacksonLaboratory/ATAC-seq/blob/master/auyar/pyadapter_trim.py """
    def __fuzz_align(self,reads,adaptor):
        if not isinstance(reads,str):
            reads = reads
        for idx, base in enumerate(reads):  # loop through equal size windows
            dist = 0
            reads_subset = reads[idx:idx+len(adaptor)]
            if len(reads_subset)<len(adaptor):
                break
            if reads_subset == adaptor:
                return idx,dist
                break
            else:
                dist = Levenshtein.distance(reads_subset,adaptor)
                if dist <= self.mismatch:  # find first then break
                    return idx,dist
                    break
                
    def _trimming(self,seq_field,reads =''):

        name, seq, qual = seq_field
        expand_for_overhang = 0
        minimal_adaptor = 15
        
        if reads == 'R1':
            seq = seq[0:self.cutoffR1]
            qual = qual[0:self.cutoffR1]
        else:
            seq = seq[-self.cutoffR2:]
            qual = qual[-self.cutoffR2:]

        for each_subset in self.adaptor_set:
            if len(each_subset) >= self.length_for_complete:
                if fuzz.partial_ratio(each_subset,seq)>self.threhold_for_Levenshtein:
                    hold = self.__fuzz_align(seq,each_subset)
                    if hold:
                        idx,dist = hold
                        if dist<=self.mismatch:
                            seq = seq[idx+len(each_subset):]
                            qual = qual[idx+len(each_subset):]
                            break
            elif len(each_subset) < self.length_for_complete and len(each_subset)>minimal_adaptor:
                line_front = seq[0:len(each_subset)+expand_for_overhang]
                #dist = Levenshtein.distance(line_front, each_subset)
                dist = fuzz.partial_ratio(each_subset,line_front)
                dist = len(each_subset)*(1-dist)
                if dist<=self.mismatch:
                    seq = seq[len(each_subset)-1+expand_for_overhang:]
                    qual = qual[len(each_subset)-1+expand_for_overhang:]
                    break
        for each_subset in self.adaptor_rev_comp_set:
            if len(each_subset) >= self.length_for_complete:
                hold = self.__fuzz_align(seq,each_subset)
                if fuzz.partial_ratio(each_subset,seq)>self.threhold_for_Levenshtein:
                    if hold:
                        idx,dist = hold
                        if dist <= self.mismatch:
                            seq = seq[0:idx]
                            qual = qual[0:idx]
                            break
            elif len(each_subset) < self.length_for_complete and len(each_subset)>minimal_adaptor:
                line_back = seq[-len(each_subset)-expand_for_overhang:]
                #dist = Levenshtein.distance(line_back, each_subset)
                dist = fuzz.partial_ratio(each_subset,line_front)
                dist = len(each_subset)*(1-dist)
                if dist<=self.mismatch:
                    seq = seq[0:-len(each_subset)-expand_for_overhang]
                    qual = qual[0:-len(each_subset)-expand_for_overhang]
        return name,seq,qual

    def __R2_Demultiplexing(self,R2_field):
        """demultiplexing R2 with specific location of barcodes"""
        R2_name,R2_seq,R2_qual = R2_field
        R2_seq = R2_seq[0:self.cell_bc_total_len]
        R2_qual = R2_qual[0:self.cell_bc_total_len]
        R2_seq_part = R2_seq.partition(self.linker)
        
        
        if len(R2_seq_part) == 3:
            if len(R2_seq_part[0]) == 51:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]:self.pos[1][1]]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]:self.pos[2][1]]
                me = R2_seq_part[2][3:22]

                if temp_cell_bc1 in self.bc1 and temp_cell_bc2 in self.bc2 and temp_cell_bc3 in self.bc3 and sample_index in self.sample_index_map and self.me == me:
                    return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
                else:
                    return None,None,None,None
            elif len(R2_seq_part[0]) == 52:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+1:self.pos[1][1]+1]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+1:self.pos[2][1]+1]
                me = R2_seq_part[2][3:22]

                if temp_cell_bc1 in self.bc1 and temp_cell_bc2 in self.bc2 and temp_cell_bc3 in self.bc3 and sample_index in self.sample_index_map and self.me == me:
                    return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
                else:
                    return None,None,None,None
            elif len(R2_seq_part[0]) == 53:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+2:self.pos[1][1]+2]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+2:self.pos[2][1]+2]
                me = R2_seq_part[2][3:22]
        
                if temp_cell_bc1 in self.bc1 and temp_cell_bc2 in self.bc2 and temp_cell_bc3 in self.bc3 and sample_index in self.sample_index_map and self.me == me:
                    return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
                else:
                    return None,None,None,None
            elif len(R2_seq_part[0]) == 54:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+3:self.pos[1][1]+3]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+3:self.pos[2][1]+3]
                me = R2_seq_part[2][3:22]
      
                if temp_cell_bc1 in self.bc1 and temp_cell_bc2 in self.bc2 and temp_cell_bc3 in self.bc3 and sample_index in self.sample_index_map and self.me == me:
                    return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
                else:
                    return None,None,None,None
            else:
                pass
            return None,None,None,None
        else:
            pass
        return None,None,None,None

    """Modified from https://github.com/lh3/readfq """
    def __readfq(self,fp): # this is a generator function
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq              
                for l in fp: # search for the start of the next record
                    l = l.rstrip()
                    if l[0] in '>@': # fasta/q header line
                        last = l # save this line, header line
                        break
            if not last: break
            name, seqs, last = last[1:], [], None
            for l in fp: # read the seq
                l = l.rstrip()
                if l[0] in '@+>':
                    last = l
                    break
                seq2 = l
            leng, seqs =  0, []
            for l in fp: # read the quality
                l = l.rstrip()
                last = None
                yield name, seq2, l; # yield a fastq record
                break
    def create_barcode_set(self):
        """create map betwen our barcode list"""
        all_barcode_list = [(temp_bc1,temp_bc2,temp_bc3) for temp_bc1 in self.bc1 for temp_bc2 in self.bc2 for temp_bc3 in self.bc3]
        self.all_sc_barcodes = set(all_barcode_list)
        
    def process_barcode(self,barcode_name):
        """get the barcode list from barcode file"""
        f_barcode = open(barcode_name,'r')
        all_barcode = np.array([line.strip() for line in f_barcode.readlines()])
        index_sample = np.where(all_barcode=="[Sample_index]")[0]
        index_bc1 = np.where(all_barcode=="[Cell_BC1]")[0]
        index_bc2 = np.where(all_barcode=="[Cell_BC2]")[0]
        index_bc3 = np.where(all_barcode=="[Cell_BC3]")[0]
     
        for i in range(len(all_barcode)):
            if i > index_sample and i < index_bc1:
                if all_barcode[i]:
                    temp_sample_bc = all_barcode[i].split()
                    self.sample_index_map[temp_sample_bc[1]] = temp_sample_bc[0]
            elif i > index_bc1 and i < index_bc2:
                if all_barcode[i]:
                    temp_sample_bc = all_barcode[i].split()
                    self.bc1.append(temp_sample_bc[1])        
            elif i > index_bc2 and i < index_bc3:
                if all_barcode[i]:
                    temp_sample_bc = all_barcode[i].split()
                    self.bc2.append(temp_sample_bc[1])
            elif i > index_bc3:
                if all_barcode[i]:
                    temp_sample_bc = all_barcode[i].split()
                    self.bc3.append(temp_sample_bc[1])
                    
        self.bc1 = set(self.bc1)
        self.bc2 = set(self.bc2)
        self.bc3 = set(self.bc3)
        
def main(kwargs):

    args = kwargs
    R1_input = args.R1
    R2_input = args.R2
    barcode_name = args.barcode
    output_dir = args.out_dir
    adaptor = args.adaptor
    cpu = args.processor
    min_length_seq = args.min_length_seq
    
    debarcoding_obj = Fastq_transform(R1_input = R1_input, R2_input = R2_input,output_dir = output_dir,
   adaptor = adaptor,cpu = cpu ,min_length_seq = min_length_seq)
    debarcoding_obj.mkdir()
    debarcoding_obj.process_barcode(barcode_name)
    debarcoding_obj.create_barcode_set()
    debarcoding_obj.init_write()
    debarcoding_obj.init_adaptors()
    debarcoding_obj.distribute_to_processor()
    debarcoding_obj.close_file()
    

if __name__ == "__main__":

    if sys.argv[-1].endswith(".py"):
        print(HELP)
    else:
        kwargs = fargv()
        main(kwargs)

