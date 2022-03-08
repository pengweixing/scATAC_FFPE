#################################################
#  File Name:fastq_to_10x.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com, pengweixing@igp.uu.se
#  Created Time: Thu Nov 25 17:31:28 2021
#################################################
HELP = """ Example:
---------------------------------------------------------------------------
R2 sequence with four structures: 5' -> 3'                                                      <--- R1         
cell_bc1              cell_bc2             cell_bc3    linker   sample_index                   
ACGATTGNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNNN   -> ME -> genomicDNA -> ME -> 
ACGATTGNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNN   -> ME -> genomicDNA -> ME ->
ACGATTGNNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNN   -> ME -> genomicDNA -> ME ->
ACGATTGNNNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNN   -> ME -> genomicDNA -> ME ->
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
usage: python debarcode.py -r1 R1.fastq.gz -r2 R2.fastq.gz -b barcode.list -o output_name -d 737K-cratac-v1.txt

optional arguments:
  -h, --help            show this help message and exit
  -r1 R1, --R1 R1       the R1 of file
  -r2 R2, --R2 R2       the R2 of file
  -l POS_OF_R1, --pos_of_R1 POS_OF_R1
                        the position of genomic DNA in R1, default is 95
  -b BARCODE, --barcode BARCODE
                        the barcode list
  -d A10X_BARCODE, --a10x_barcode A10X_BARCODE
                        the 10x barcode list
  -o OUT_DIR, --out_dir OUT_DIR
                        the directory for output
  -fi MIN_FRAGS_CUTOFF_FOR_PLOT, --min_frags_cutoff_for_plot MIN_FRAGS_CUTOFF_FOR_PLOT
                        the minimum number of frags for plot
  -fa MAX_FRAGS_CUTOFF_FOR_PLOT, --max_frags_cutoff_for_plot MAX_FRAGS_CUTOFF_FOR_PLOT
                        the maximum number of frags for plot
  -bins HIST_BINS, --hist_bins HIST_BINS
                        the number of bins for hist plot
  -wt WIDTH, --width WIDTH
                        the width of figure
  -ht HEIGHT, --height HEIGHT
                        the height of figure
  -a ADAPTOR, --adaptor ADAPTOR
                        the adaptor's sequence
  -p PROCESSOR, --processor PROCESSOR
                        the number of processors
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
import random
import time
import pickle
import Levenshtein
from matplotlib import pyplot as plt
from multiprocessing import Pool
import threading
from itertools import islice

def fargv():
    parser = argparse.ArgumentParser(usage="python debarcode.py -r1 R1.fastq.gz -r2 R2.fastq.gz -b barcode.list -d 737K-cratac-v1.txt -o output_name")
    parser.add_argument('-r1',"--R1",help="the R1 of file ", required=True)
    parser.add_argument('-r2',"--R2",help="the R2 of file ", required=True)
    parser.add_argument('-b',"--barcode",help="the barcode list ", required=True)
    parser.add_argument('-d',"--a10x_barcode",help="the 10x barcode list ", required=True)
    parser.add_argument('-o',"--out_dir",help="the directory for output", default='./')
    parser.add_argument('-fi',"--min_frags_cutoff_for_plot",help="the minimum number of frags for plot",default=2,type=int)
    parser.add_argument('-fa',"--max_frags_cutoff_for_plot",help="the maximum number of frags for plot",default=100000,type=int)
    parser.add_argument('-bins',"--hist_bins",help="the number of bins for hist plot",default=100,type=int)
    parser.add_argument('-wt',"--width",help="the width of figure",default=10,type=int)
    parser.add_argument('-ht',"--height",help="the height of figure",default=12,type=int)
    parser.add_argument('-a',"--adaptor",help="the adaptor's sequence",type=str,default='AGATGTGTATAAGAGACAG')
    parser.add_argument('-p',"--processor",help="the number of processors",default=10,type=int)
    args = parser.parse_args()
    return args

def process_R1_R2_wrapper(Fastq_transform,data_R1_each_chunk,data_R2_each_chunk):
    return Fastq_transform.process_R1_R2(data_R1_each_chunk,data_R2_each_chunk)

class Fastq_transform:
    names = locals()
    def __init__(self,R1_input = [], R2_input = [],barcode_10x_file = [],output_dir = './', min_frags_cutoff_plot = 2,\
    max_frags_cutoff_plot = 100000, bins=100, width = 10, height=12, adaptor = '',min_length_seq = 30,cpu = 10):

        self.R1_input = R1_input
        self.R2_input = R2_input
        self.pos = [(0,7),(22,29),(44,51)]
        self.barcode_database = barcode_10x_file
        """linker sequence between sample index and cell BC1"""
        self.linker = "ATCCACGAGCATTCG"
        self.sample_index_map,self.bc1,self.bc2,self.bc3 = {},[],[],[]
        self.all_sc_map_10x = []
        self.output_dir = output_dir
        self.min_frags_cutoff_plot = min_frags_cutoff_plot
        self.max_frags_cutoff_plot = max_frags_cutoff_plot
        self.plot_bins = bins
        self.width = width
        self.height = height
        self.adaptor = adaptor
        self.threhold_for_Levenshtein = 0.8
        self.adaptor_set = []
        self.adaptor_rev_comp_set = []
        self.mismatch = 3
        self.min_length_seq = min_length_seq
        self.cutoffR1 = 75
        self.cutoffR2 = 75
        self.mychunksize = 1000000  #### read 5M lines per time and distribute them to multi cpus
        self.Processor = cpu

    def init_adaptors(self):
        adaptor_comp = self.DNA_complement(self.adaptor)
        adaptor_rev_comp = self.DNA_reverse(adaptor_comp)
        self.adaptor_set = self.__gen_adaptor(self.adaptor,rev=False)
        self.adaptor_rev_comp_set = self.__gen_adaptor(adaptor_rev_comp,rev=True)

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
        
        data_R1 = pd.read_csv(self.R1_input,compression='gzip',sep="\t",chunksize=self.mychunksize,header=None) 
      #  data_R1 = gzip.open(self.R1_input,'rb')
        data_R2 = pd.read_csv(self.R2_input,compression='gzip',sep="\t",chunksize=self.mychunksize,header=None)
       # data_R2 = gzip.open(self.R2_input,'rb')
      #  lines_gen = islice(f1, 10)
      #  lines = [line.strip().decode() for line in lines_gen]

        Barcode_reads_total = {}
        cell_number_each_total = {}
        Total_reads = 0
        Barcode_reads_total_C =  Counter(Barcode_reads_total)
        cell_number_each_total_C =  Counter(cell_number_each_total)
        for data_R1_chunk,data_R2_chunk in zip(data_R1,data_R2):
            data_R1_chunk = data_R1_chunk.values.tolist()
            data_R2_chunk = data_R2_chunk.values.tolist()
            pool = Pool(self.Processor)
            reads = self.mychunksize/4        ## The total number of reads
            block = int(reads/self.Processor) ## The number of reads for each processor
            results = []
            ## split the total reads into different part and assign them to each cluster
            for index in range(self.Processor): 
                start = index*block*4
                end = (index+1)*block*4
                data_R1_each_chunk = data_R1_chunk[start:end]
                data_R2_each_chunk = data_R2_chunk[start:end]
                results.append(pool.apply_async(process_R1_R2_wrapper,args = (self,data_R1_each_chunk,data_R2_each_chunk,)\
                ,error_callback=self.print_error))    
            pool.close()
            pool.join()
            for result in results:
                result = result.get()
           #     print(result[0])
                """
                result[0]: R1
                result[1]: cell barcodes, R2
                result[2]: R3 
                """
                #### write processed reads to different files with different threads
                x1 = threading.Thread(target=self.__write_to_file_R1_R3, args=(result[0],result[1],'R1',))
                x2 = threading.Thread(target=self.__write_to_file_R2, args=(result[0],result[1],))
                x3 = threading.Thread(target=self.__write_to_file_R1_R3, args=(result[2],result[1],'R3',))
                x1.start(),x2.start(),x3.start()
                x1.join(),x2.join(),x3.join()        
                ### merge all dicts of cell * fragments and reads number * sample 
                Barcode_reads_each_C = Counter(result[3])
                cell_number_each_C = Counter(result[4])
                Total_reads_each = result[5]
                Total_reads = Total_reads + Total_reads_each
                Barcode_reads_total_C = Counter(dict(Barcode_reads_total_C+Barcode_reads_each_C))
                cell_number_each_total_C = Counter(dict(cell_number_each_total_C+cell_number_each_C))
            print("%s reads have been done" % str(i*self.mychunksize/4))
            i += 1
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
        buffer_R2 = []
        buffer_R3 = []
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
            if sample_index in self.sample_index_map and len(R1_seq) >= self.min_length_seq and len(R2_seq) >= self.min_length_seq and\
            (temp_cell_bc1,temp_cell_bc2,temp_cell_bc3) in self.all_sc_map_10x:
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
                
                ### keep the trimed and decoded reads in the list
                buffer_R1.append((R1_name+'\n',R1_seq,R1_qual))
                buffer_R2.append((sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3))
                buffer_R3.append((R2_name+'\n',R2_seq,R2_qual))
        return buffer_R1,buffer_R2,buffer_R3,Barcode_reads,cell_number_order_dict,Total_reads


    def hist_plot(self,all_line = []):
        data = pd.DataFrame(all_line)
        data.columns = ["sample","bc1","bc2","bc3","number"]
        data['number']=data['number'].astype(int)
        data = data.loc[data['number'] >= self.min_frags_cutoff_plot]
        data = data.loc[data['number'] <= self.max_frags_cutoff_plot]
        all_sampels = set(list(data['sample']))
        number = len(all_sampels)     
        i = 0
        j = 0
        ###  create the index of the subplot within compound graph
        comb = []
        for ii in range(number):
            if ii%2 == 0:
                comb.append((i,0))
            else:
                comb.append((i,1))
                i+=1
        ii = 0 
        fig, ax = plt.subplots(math.ceil(number/2),2,figsize=(self.width,self.height),frameon=False)   

        for each in all_sampels:
            each_data = pd.Series(data.loc[data['sample']==each]['number'])
            ax[comb[ii][0],comb[ii][1]].hist(each_data,bins=self.plot_bins)
            ax[comb[ii][0],comb[ii][1]].set_yscale("log")
            ax[comb[ii][0],comb[ii][1]].set_xlabel("The number of fragments for each cell")
            ax[comb[ii][0],comb[ii][1]].set_ylabel("The number of cells")
            ax[comb[ii][0],comb[ii][1]].set_title(each)
            ii += 1
            plt.tight_layout()
        fig.savefig(self.output_dir+"/Fragments_over_cells_distri.pdf",format='pdf',bbox_inches='tight')

    def create_barcode_map(self,barcode_database):
        """create map betwen our barcode list and 10x barcode list"""
        all_barcode_list = [(temp_bc1,temp_bc2,temp_bc3) for temp_bc1 in self.bc1 for temp_bc2 in self.bc2 for temp_bc3 in self.bc3]
        f_10x = open(barcode_database,'r')
        all_10x_barcode = [line.strip() for line in f_10x.readlines()]
        all_10x_barcode = all_10x_barcode[0:len(all_barcode_list)]
        self.all_sc_map_10x = dict(zip(all_barcode_list,all_10x_barcode))
       
    def init_write(self):
        """initialize the output file"""
        for outname in self.sample_index_map.values():
        #    Fastq_transform.names['fI'+outname] = gzip.open(self.output_dir+"/"+outname+'_S1_L001_I1_001.fastq.gz','wb')  ### sample index
            Fastq_transform.names['f1'+outname] = gzip.open(self.output_dir+"/"+outname+'_S1_L001_R1_001.fastq.gz','wb')   ### Read 1 of genomic DNA
            Fastq_transform.names['f2'+outname] = gzip.open(self.output_dir+"/"+outname+'_S1_L001_R2_001.fastq.gz','wb')   ### single cell barcode
            Fastq_transform.names['f3'+outname] = gzip.open(self.output_dir+"/"+outname+'_S1_L001_R3_001.fastq.gz','wb')   ### Read 2 of genomic DNA
    
    def close_file(self):
        for outname in self.sample_index_map.values():
         #   Fastq_transform.names['fI'+outname].close()
            Fastq_transform.names['f1'+outname].close()
            Fastq_transform.names['f2'+outname].close()
            Fastq_transform.names['f3'+outname].close()

    def __write_to_file_R1_R3(self,buffer_R1_R3,buffer_R2,which_reads):
        """ write the proccessed reads to fastq with 10x format"""
        for R1_record,R2_record in zip(buffer_R1_R3,buffer_R2):
            R1_name,R1_seq,R1_qual = R1_record
            sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3 = R2_record
            outname = ''
            ten_x_barcode = ''
            if sample_index in self.sample_index_map:  ## get the sample name from each sample index
                outname = self.sample_index_map[sample_index] 
            if R2_record[1:4] in self.all_sc_map_10x:   ## get the 10x barcode from each triple sc barcode
                ten_x_barcode = self.all_sc_map_10x[R2_record[1:4]]
            temp_for_write_I1 = {}
            temp_for_write_R1 = {}
                    

            if outname and ten_x_barcode: 
                if outname not in temp_for_write_I1.keys():
                    temp_for_write_I1[outname]  = '@' + R1_name + 'AAAGCATA\n' + '+\n' + 'FFFFFFFF\n'
                    temp_for_write_R1[outname]  = '@' + R1_name + R1_seq + '\n' + '+\n' + R1_qual + '\n'
                else:
                    temp_for_write_I1[outname]  = temp_for_write_I1[outname] + '@' + R1_name + 'AAAGCATA\n' + '+\n' + 'FFFFFFFF\n'
                    temp_for_write_R1[outname]  = temp_for_write_R1[outname] + '@' + R1_name + R1_seq + '\n' + '+\n' + R1_qual + '\n'
    
            for outname in temp_for_write_I1.keys():        
                if which_reads == "R1":
                    Fastq_transform.names['f1'+outname].write(temp_for_write_R1[outname].encode())
                elif which_reads == "R3":
                    Fastq_transform.names['f3'+outname].write(temp_for_write_R1[outname].encode())

    def __write_to_file_R2(self,buffer_R1,buffer_R2):
        """ write the proccessed reads to fastq with 10x format"""
        for R1_record,R2_record in zip(buffer_R1,buffer_R2):
            R1_name,R1_seq,R1_qual = R1_record
            sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3 = R2_record
            outname = ''
            ten_x_barcode = '' 
            if sample_index in self.sample_index_map:  ## get the sample name from each sample index
                outname = self.sample_index_map[sample_index] 
            if R2_record[1:4] in self.all_sc_map_10x:
               ## get the 10x barcode from each sc barcode
                ten_x_barcode = self.all_sc_map_10x[R2_record[1:4]]
            temp_for_write_I1 = {}
            temp_for_write_R2 = {}
            if outname and ten_x_barcode: 
                if outname not in temp_for_write_I1.keys():
                    temp_for_write_I1[outname]  = '@' + R1_name + 'AAAGCATA\n' + '+\n' + 'FFFFFFFF\n'
                    temp_for_write_R2[outname]  = '@' + R1_name + ten_x_barcode + '\n' + '+\n' + 'FFFFFFFFFFFFFFFF\n'
                else:
                    temp_for_write_I1[outname]  = temp_for_write_I1[outname] + '@' + R1_name + 'AAAGCATA\n' + '+\n' + 'FFFFFFFF\n'
                    temp_for_write_R2[outname]  = temp_for_write_R2[outname] + '@' + R1_name + ten_x_barcode + '\n' + '+\n' + 'FFFFFFFFFFFFFFFF\n'            
            for outname in temp_for_write_I1.keys():        
                Fastq_transform.names['f2'+outname].write(temp_for_write_R2[outname].encode())

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
      #  print(self.adaptor_set)
        name, seq, qual = seq_field
        if reads == 'R1':
            seq = seq[0:self.cutoffR1]
            qual = qual[0:self.cutoffR1]
        else:
            seq = seq[self.cutoffR2:]
            qual = qual[self.cutoffR2:]

        for each_ap in self.adaptor_set: ### if the direction of adaptor is 5'->3'
         #   print(len(each_ap),len(self.adaptor),sep="\t")
            if len(each_ap) == len(self.adaptor): #if the adaptor exit in the sequence with full length, it can appear at anywhere
                hold = self.__fuzz_align(seq,each_ap)
                if hold:
                    idx,dist = hold
                    seq = seq[idx+len(each_ap):]
                    qual = qual[idx+len(each_ap):]
            else:  ## if the adaptor is not intact, then it should be at 5' of the seq
                seq_5p = seq[0:len(each_ap)]
                dist = Levenshtein.distance(seq_5p,each_ap)
                if dist <= self.mismatch:
                    seq = seq[len(seq_5p)+1:] 
                    qual = qual[len(seq_5p)+1:] 
        """if the direction of adaptor is 3'->5', which means the genomic DNA is too short, the another side of DNA adaptor also is sequenced"""
        for each_ap_rc in self.adaptor_rev_comp_set:  

            if len(each_ap) == len(self.adaptor): #if the adaptor exit in the sequence with full length, it can appear at anywhere
                hold = self.__fuzz_align(seq,each_ap)
                if hold:
                    idx,dist = hold
                    seq = seq[:idx]
                    qual = qual[:idx]
            else:  ## if the adaptor is not intact, then it should be at 3' of the reads
                seq_5p = seq[-len(each_ap):]
                dist = Levenshtein.distance(seq_5p,each_ap)
                if dist <= self.mismatch:
                    seq = seq[0:-len(seq_5p)] 
                    qual = qual[0:-len(seq_5p)] 

        return name,seq,qual

    def __R2_Demultiplexing(self,R2_field):
        """demultiplexing R2 with specific location of barcodes"""
        R2_name,R2_seq,R2_qual = R2_field
        R2_seq = R2_seq[0:self.cutoffR2]
        R2_qual = R2_qual[0:self.cutoffR2]
        R2_seq_part = R2_seq.partition(self.linker)
        if len(R2_seq_part) == 3:
            if len(R2_seq_part[2]) == 9:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]:self.pos[1][1]]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]:self.pos[2][1]]
                return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
            elif len(R2_seq_part[2]) == 8:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+1:self.pos[1][1]+1]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+1:self.pos[2][1]+1]
                return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
            elif len(R2_seq_part[2]) == 7:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+2:self.pos[1][1]+2]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+2:self.pos[2][1]+2]
                return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
            elif len(R2_seq_part[2]) == 6:
                sample_index = R2_seq_part[2][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+3:self.pos[1][1]+3]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+3:self.pos[2][1]+3]
                return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
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
                    l = l[0]
                    if l[0] in '>@': # fasta/q header line
                        last = l # save this line, header line
                        break
            if not last: break
            name, seqs, last = last[1:], [], None
            for l in fp: # read the seq
                seq = l[0]
                if l[0] in '@+>':
                    last = l
                    break
                seq2 = seq
            leng, seqs =  0, []
            for l in fp: # read the quality
                seqs = l[0]
                last = None
                yield name, seq2, seqs; # yield a fastq record
                break

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
            
def main(kwargs):

    args = kwargs
    R1_input = args.R1
    R2_input = args.R2
    barcode_name = args.barcode
    barcode_10x_file = args.a10x_barcode
    output_dir = args.out_dir
    min_frags_cutoff_plot = args.min_frags_cutoff_for_plot
    max_frags_cutoff_plot = args.max_frags_cutoff_for_plot
    hist_bins = args.hist_bins
    width = args.width
    height = args.height
    adaptor = args.adaptor
    cpu = args.processor
    

    debarcoding_obj = Fastq_transform(R1_input = R1_input, R2_input = R2_input,output_dir = output_dir, min_frags_cutoff_plot = min_frags_cutoff_plot,\
    max_frags_cutoff_plot = max_frags_cutoff_plot, bins=hist_bins, width = width,height=height, adaptor = adaptor,cpu = cpu )
    debarcoding_obj.mkdir()
    debarcoding_obj.process_barcode(barcode_name)
    debarcoding_obj.create_barcode_map(barcode_10x_file)
    debarcoding_obj.init_write()
    debarcoding_obj.init_adaptors()
    debarcoding_obj.distribute_to_processor()
    debarcoding_obj.close_file()
  #  

if __name__ == "__main__":

    if sys.argv[-1].endswith(".py"):
        print(HELP)
    else:
        kwargs = fargv()
        main(kwargs)

