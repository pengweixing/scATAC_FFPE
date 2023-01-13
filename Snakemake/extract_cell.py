#################################################
#  File Name:extract_cell.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Wed Nov 30 19:36:03 2022
#################################################

import pysam
import sys
import argparse
import gzip

def fargv():
    parser = argparse.ArgumentParser(usage="python ")
    parser.add_argument('-n',"--names",help="the cellnames file ", required=True)
    parser.add_argument('-f',"--fragments",help="the fragments file ", required=True)
    parser.add_argument('-b',"--bam",help="the bam file ", required=True)
    parser.add_argument('-of',"--out_frags",help="the output of fragments ", required=True)
    parser.add_argument('-ob',"--out_bam",help="the output of bam ",required=True)

    args = parser.parse_args()
    return args

def get_frags(all_names,kwargs):
    f = gzip.open(kwargs.fragments,'rb')
    f_out = open(kwargs.out_frags,'w')
    for line in f:
        line = line.strip().decode()
        line1 = line.split()
        if line1[3] in all_names:
            f_out.write("%s\n" % line)
    f.close()
    f_out.close()
    
def get_bam(all_names,kwargs):
    f_bamin = pysam.AlignmentFile(kwargs.bam,mode='rb',check_header=True,threads=3)
    f_bamout = pysam.AlignmentFile(kwargs.out_bam,mode='wb',template=f_bamin,threads=3)
    for a in f_bamin:
        cellname = a.get_tag('RG')
        if cellname in all_names:
            f_bamout.write(a)
    f_bamin.close()
    f_bamout.close(
        
    )
    
    
def main(kwargs):

    f_names = open(kwargs.names,'r')
    all_names = set([line.strip() for line in f_names.readlines()])
    get_frags(all_names,kwargs)
    get_bam(all_names,kwargs)
    
if __name__ == "__main__":
    kwargs = fargv()
    main(kwargs)
