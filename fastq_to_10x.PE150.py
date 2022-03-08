#################################################
#  File Name:fastq_to_10x.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com, pengweixing@igp.uu.se
#  Created Time: Thu Nov 25 17:31:28 2021
#################################################
HELP = """ Example:
---------------------------------------------------------------------------
R2 sequence with four structures: 5' -> 3'
cell_bc1              cell_bc2             cell_bc3    linker   sample_index
ACGATTGNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNNN   -> ME -> genomicDNA
ACGATTGNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNNN   -> ME -> genomicDNA
ACGATTGNNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNNN   -> ME -> genomicDNA
ACGATTGNNNNNNNNNNNNNNNNNNAAACCGGNNNNNNNNNNNNNNNACCCTAANNNNNNNNNNNNNNNAGANNN   -> ME -> genomicDNA
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

[Cell_BC1]
R02_#01 AAACCGG
R02_#02 AAACGTC
R02_#03 AAAGATG

[Cell_BC2]
R03_#01 AAACCGG
R03_#02 AAACGTC
R03_#03 AAAGATG

[Cell_BC3]
R04_#01 AAACCGG
R04_#02 AAACGTC
R04_#03 AAAGATG
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
---------------------------------------------------------------------------

"""

import sys
import gzip
import argparse
import collections
import numpy as np
import os
import pandas as pd
import math
import Levenshtein
from matplotlib import pyplot as plt

def fargv():
    parser = argparse.ArgumentParser(usage="python debarcode.py -r1 R1.fastq.gz -r2 R2.fastq.gz -b barcode.list -d 737K-cratac-v1.txt -o output_name")
    parser.add_argument('-r1',"--R1",help="the R1 of file ", required=True)
    parser.add_argument('-r2',"--R2",help="the R2 of file ", required=True)
    parser.add_argument('-l',"--pos_of_R1",help="the position of genomic DNA in R1, default is 95", default="95",type=int)
    parser.add_argument('-b',"--barcode",help="the barcode list ", required=True)
    parser.add_argument('-d',"--a10x_barcode",help="the 10x barcode list ", required=True)
    parser.add_argument('-o',"--out_dir",help="the directory for output", default='./')
    parser.add_argument('-fi',"--min_frags_cutoff_for_plot",help="the minimum number of frags for plot",default=2,type=int)
    parser.add_argument('-fa',"--max_frags_cutoff_for_plot",help="the maximum number of frags for plot",default=100000,type=int)
    parser.add_argument('-bins',"--hist_bins",help="the number of bins for hist plot",default=100,type=int)
    parser.add_argument('-wt',"--width",help="the width of figure",default=10,type=int)
    parser.add_argument('-ht',"--height",help="the height of figure",default=12,type=int)
    parser.add_argument('-a',"--adaptor",help="the adaptor's sequence",type=str,default='AGATGTGTATAAGAGACAG')
    args = parser.parse_args()
    return args

class Fastq_transform:
    names = locals()
    def __init__(self,R1_input = [], R2_input = [],barcode_10x_file = [],output_dir = './', min_frags_cutoff_plot = 2,\
    max_frags_cutoff_plot = 100000, bins=100, width = 10, height=12, adaptor = '',min_length_seq = 30):

        self.R1_input = R1_input
        self.R2_input = R2_input
        self.pos = [(0,7),(22,29),(44,51)]
        self.barcode_database = barcode_10x_file
        """linker sequence between sample index and cell BC1"""
        self.linker = "ATCCACGAGCATTCG"
        self.sample_index_map,self.bc1,self.bc2,self.bc3 = {},[],[],[]
        self.all_sc_map_10x = []
        self.cache_records = 100000 ## fastqs for write
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

    def process_R1_R2(self):
        f_R1 = gzip.open(self.R1_input,'rt')
        f_R2 = gzip.open(self.R2_input,'rt')
        ii = 1
        Total_reads = 0
        Barcode_reads = {}
        buffer_R1 = []
        buffer_R2 = []
        buffer_R3 = []
        cell_number_order_dict = collections.OrderedDict()
        for R1_field,R2_field in zip(self.__readfq(f_R1),self.__readfq(f_R2)):
            R1_name,R1_seq,R1_qual = self._trimming(R1_field,reads ='R1')
            R2_name,R2_seq,R2_qual = self._trimming(R2_field,reads = 'R2')
            sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3 = self.__R2_Demultiplexing(R2_field)      
            Total_reads += 1
            if sample_index and len(R1_seq) >= self.min_length_seq and len(R2_seq) >= self.min_length_seq:
                ### count the reads which have correct barcodes
                if sample_index in self.sample_index_map.keys():
                    if self.sample_index_map[sample_index] in Barcode_reads.keys():  
                        Barcode_reads[self.sample_index_map[sample_index]] += 1
                    else:
                        Barcode_reads[self.sample_index_map[sample_index]] = 1
                """count the fragments number for each cell"""
                if (sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3) in cell_number_order_dict.keys():
                    if sample_index in self.sample_index_map.keys():
                        cell_number_order_dict[(sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3)] += 1
                else:
                    if sample_index in self.sample_index_map.keys():
                        cell_number_order_dict[(sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3)] = 1
                
                #### To reduce the time-consuming, we write 10000 reads to file per time
                if ii % self.cache_records == 0:
                    buffer_R1.append((R1_name+'\n',R1_seq,R1_qual))
                    buffer_R2.append((sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3))
                    buffer_R3.append((R2_name+'\n',R2_seq,R2_qual))
                    self.__write_to_file(buffer_R1, buffer_R2, buffer_R3)
                    ii = 1
                    buffer_R1 = []
                    buffer_R2 = []
                    buffer_R3 = []
                else:                
                    buffer_R1.append((R1_name+'\n',R1_seq,R1_qual))
                    buffer_R2.append((sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3))
                    buffer_R3.append((R2_name+'\n',R2_seq,R2_qual))
                    ii = ii + 1

        """write the fragments number for each cell to file"""
        f = open(self.output_dir+"/Barcodes_fragments_qc.txt",'w')
        f.write('sample\tBC1\tBC2\tBC3\tFrags_number\n')
        all_stat = []
        for key,value in cell_number_order_dict.items():
            all_stat.append([self.sample_index_map[key[0]],key[1],key[2],key[3],value])
            f.write("%s\t%s\t%s\t%s\t%s\n" % (self.sample_index_map[key[0]],key[1],key[2],key[3],value))
        f.close()

        """write the statistics of debarcoding rate to file """
        f = open(self.output_dir+"/Barcoding_rate_qc.txt",'w')
        f.write('Total_reads\t%s\n' % Total_reads)
        f.write('sample\tFinal_reads\n')
        for key in Barcode_reads.keys():
            barcoded_reads = Barcode_reads[key]
            f.write('%s\t%s\n' % (key,barcoded_reads))
        f.close()
        return all_stat

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

    def __write_to_file(self,buffer_R1, buffer_R2,buffer_R3):
        """ write the proccessed reads to fastq with 10x format"""
        for R1_record,R2_record,R3_record in zip(buffer_R1,buffer_R2,buffer_R3):
            R1_name,R1_seq,R1_qual = R1_record
            sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3 = R2_record
            R3_name,R3_seq,R3_qual = R3_record
            outname = ''
            ten_x_barcode = ''
            if sample_index in self.sample_index_map:  ## get the sample name from each sample index
                outname = self.sample_index_map[sample_index] 
            if R2_record[1:4] in self.all_sc_map_10x:   ## get the 10x barcode from each triple sc barcode
                ten_x_barcode = self.all_sc_map_10x[R2_record[1:4]]
            temp_for_write_I1 = {}
            temp_for_write_R1 = {}
            temp_for_write_R2 = {}
            temp_for_write_R3 = {}

            if outname and ten_x_barcode: 
                if outname not in temp_for_write_I1.keys():
                    temp_for_write_I1[outname]  = '@' + R1_name + 'AAAGCATA\n' + '+\n' + 'FFFFFFFF\n'
                    temp_for_write_R1[outname]  = '@' + R1_name + R1_seq + '\n' + '+\n' + R1_qual + '\n'
                    temp_for_write_R2[outname]  = '@' + R1_name + ten_x_barcode + '\n' + '+\n' + 'FFFFFFFFFFFFFFFF\n'
                    temp_for_write_R3[outname]  = '@' + R3_name + R3_seq + '\n' + '+\n' + R3_qual + '\n'
                else:
                    temp_for_write_I1[outname]  = temp_for_write_I1[outname] + '@' + R1_name + 'AAAGCATA\n' + '+\n' + 'FFFFFFFF\n'
                    temp_for_write_R1[outname]  = temp_for_write_R1[outname] + '@' + R1_name + R1_seq + '\n' + '+\n' + R1_qual + '\n'
                    temp_for_write_R2[outname]  = temp_for_write_R2[outname] + '@' + R1_name + ten_x_barcode + '\n' + '+\n' + 'FFFFFFFFFFFFFFFF\n'
                    temp_for_write_R3[outname]  = temp_for_write_R3[outname] + '@' + R3_name + R3_seq + '\n' + '+\n' + R3_qual + '\n'
            
            for outname in temp_for_write_I1.keys():        
            #    Fastq_transform.names['fI'+outname].write(temp_for_write_I1[outname].encode())
                Fastq_transform.names['f1'+outname].write(temp_for_write_R1[outname].encode())
                Fastq_transform.names['f2'+outname].write(temp_for_write_R2[outname].encode())
                Fastq_transform.names['f3'+outname].write(temp_for_write_R3[outname].encode())
    
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
    def fuzz_align(self,reads,adaptor):
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
                hold = self.fuzz_align(seq,each_ap)
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
                hold = self.fuzz_align(seq,each_ap)
                if hold:
                    idx,dist = hold
                    seq = seq[:idx]
                    qual = qual[:idx]
            else:  ## if the adaptor is not intact, then it should be at 3' of the seq
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
        R2_seq_part = R2_seq.split(self.linker)
        if len(R2_seq_part) == 2:
            if len(R2_seq_part[1]) == 9:
                sample_index = R2_seq_part[1][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]:self.pos[1][1]]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]:self.pos[2][1]]
                return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
            elif len(R2_seq_part[1]) == 8:
                sample_index = R2_seq_part[1][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+1:self.pos[1][1]+1]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+1:self.pos[2][1]+1]
                return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
            elif len(R2_seq_part[1]) == 7:
                sample_index = R2_seq_part[1][0:3]
                temp_cell_bc1 = R2_seq_part[0][self.pos[0][0]:self.pos[0][1]]
                temp_cell_bc2 = R2_seq_part[0][self.pos[1][0]+2:self.pos[1][1]+2]
                temp_cell_bc3 = R2_seq_part[0][self.pos[2][0]+2:self.pos[2][1]+2]
                return sample_index,temp_cell_bc1,temp_cell_bc2,temp_cell_bc3
            elif len(R2_seq_part[1]) == 6:
                sample_index = R2_seq_part[1][0:3]
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
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
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
    

    debarcoding_obj = Fastq_transform(R1_input = R1_input, R2_input = R2_input,output_dir = output_dir, min_frags_cutoff_plot = min_frags_cutoff_plot,\
    max_frags_cutoff_plot = max_frags_cutoff_plot, bins=hist_bins, width = width,height=height, adaptor = adaptor )
    debarcoding_obj.mkdir()
    debarcoding_obj.process_barcode(barcode_name)
    debarcoding_obj.create_barcode_map(barcode_10x_file)
    debarcoding_obj.init_write()
    debarcoding_obj.init_adaptors()
    all_stat = debarcoding_obj.process_R1_R2()
    debarcoding_obj.close_file()
    debarcoding_obj.hist_plot(all_stat)

if __name__ == "__main__":

    if sys.argv[-1].endswith(".py"):
        print(HELP)
    else:
        kwargs = fargv()
        main(kwargs)

