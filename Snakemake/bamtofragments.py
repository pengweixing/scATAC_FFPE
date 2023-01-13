#################################################
#  File Name:bamtofragments.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Tue Nov 29 15:46:09 2022
#################################################

import pysam
import sys
import argparse

def fargv():
    parser = argparse.ArgumentParser(usage="python ")
    parser.add_argument('-i',"--input",help="the name-sorted bam file ", required=True)
    parser.add_argument('-of',"--frag_output",help="the fragments file ", required=True)
    parser.add_argument('-ob',"--bam_output",help="the outputed bam file ", required=True)
    parser.add_argument('-m',"--mode",help="the outputed bam file ", required=True, choices=['SE','PE'])
    args = parser.parse_args()
    return args

class convertbamtofrags:
    def __init__(self,f_bam_input='',f_frag_output='',f_bam_output=''):
        self.bam_input = f_bam_input
        self.fragoutput = f_frag_output
        self.f_bam_output = f_bam_output
        
    def tofragments_se(self):

        def assign(header,query_name,mychr,left,right,cell_name):
            temp_bam = pysam.AlignedSegment(header)
            temp_bam.query_name = query_name
            temp_bam.flag = 0
            temp_bam.reference_name = mychr
            temp_bam.reference_start = left
            temp_bam.mapping_quality = 255
            length = right-left
            temp_bam.cigar = ((0,length),)
            temp_bam.next_reference_name = "="
            temp_bam.next_reference_start = 0
            temp_bam.template_length = 0
            temp_bam.tags = (("RG", cell_name),)
            return temp_bam

        for a in self.bam_input:
            query_name = a.query_name
            cell_name = "_".join([query_name.split('_')[0]]+query_name.split('_')[5:])
            read_start = a.reference_start
            read_end = read_start + a.query_alignment_end
            mychr = a.reference_name
            self.fragoutput.write("%s\t%s\t%s\t%s\n" % (mychr,read_start,read_end,cell_name))
            temp_bam = assign(self.bam_input.header,query_name,mychr,read_start,read_end,cell_name)
            self.f_bam_output.write(temp_bam)

    def tofragments_pe(self):
                
        def assign(header,query_name,mychr,left,right,cell_name):
            temp_bam = pysam.AlignedSegment(header)
            temp_bam.query_name = query_name
            temp_bam.flag = 0
            temp_bam.reference_name = mychr
            temp_bam.reference_start = left
            temp_bam.mapping_quality = 255
            length = right-left
            temp_bam.cigar = ((0,length),)
            temp_bam.next_reference_name = "="
            temp_bam.next_reference_start = 0
            temp_bam.template_length = 0
            temp_bam.tags = (("RG", cell_name),)
            return temp_bam
        index = -1
        temp = []
        temp_name = ''
        temp_chr = ''
        temp_cell = ''
        for a in self.bam_input:
            if index == -1:
                query_name = a.query_name
                cell_name = "_".join([query_name.split('_')[0]]+query_name.split('_')[5:])
                read_start = a.reference_start
                read_end = read_start + a.query_alignment_end
                temp.extend([read_start,read_end])
                index  = index * -1
                temp_name = query_name
                mychr = a.reference_name
                temp_chr = a.reference_name
                temp_cell = cell_name
                
            else:
                if a.query_name == temp_name:
                    query_name = a.query_name
                    mychr = a.reference_name
                    cell_name = "_".join([query_name.split('_')[0]]+query_name.split('_')[5:])
                    read_start = a.reference_start
                    read_end = read_start + a.query_alignment_end
                    temp.extend([read_start,read_end])
                    index = index * -1
                    left = min(temp)
                    right = max(temp)
                    if (right - left) < 1000:
                        self.fragoutput.write("%s\t%s\t%s\t%s\n" % (mychr,left,right,cell_name))
                        temp_bam = assign(self.bam_input.header,query_name,mychr,left,right,cell_name)
                        self.f_bam_output.write(temp_bam)
                        temp = []
                        
                else:
                    self.fragoutput.write("%s\t%s\t%s\t%s\n" % (temp_chr,temp[0],temp[1],temp_cell))
                    temp_bam = assign(self.bam_input.header,query_name,temp_chr,left,right,cell_name)
                    self.f_bam_output.write(temp_bam)
                    temp = []
                    temp_name = a.query_name
                    mychr = a.reference_name
                    temp_chr = a.reference_name
                    temp_cell = cell_name
                    read_start = a.reference_start
                    read_end = read_start + a.query_alignment_end
                    temp.extend([read_start,read_end])
                    index = 1

def main(kwargs):
    f_bam_input = pysam.AlignmentFile(kwargs.input,mode='rb',check_header=True,threads=3)
    f_bam_output = pysam.AlignmentFile(kwargs.bam_output,mode='wb',template=f_bam_input,threads=3)
    f_frag_output = open(kwargs.frag_output,'w')
    bamobject = convertbamtofrags(f_bam_input=f_bam_input,f_frag_output=f_frag_output,f_bam_output=f_bam_output)
    if kwargs.mode == "PE":
        bamobject.tofragments_pe()
    elif kwargs.mode == "SE":
        bamobject.tofragments_se()
        
if __name__ == "__main__":
    kwargs = fargv()
    main(kwargs)
