#################################################
#  File Name:filter.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Thu Dec 23 10:50:17 2021
#################################################

import gzip
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(usage="python filter_fragments.py  -f $i -q $csv -mi 500 -ma 100000 -o $out1 ")
parser.add_argument('-f',"--fragments",help="fragments.tsv.gz", required=True)
parser.add_argument('-q',"--qc",help="peak_annotation.tsv", required=True)
parser.add_argument('-o',"--output",help="fragments output", required=True)
parser.add_argument('-mi',"--min_frag",help="the minumum of fragments number", default=500,type=int)
parser.add_argument('-ma',"--max_frag",help="the maximum of fragments number", default=100000,type=int)
args = parser.parse_args()

f1 = gzip.open(args.fragments,'rb') ##gzip
f2 = open(args.qc,'r') ##singlecell.csv

all_frag = {}
for line in f1:
    line = line.decode().strip()
    line1 = line.split()
    if not line.startswith('#'):
        line1[1]=int(line1[1])
        line1[2]=int(line1[2])
        if line1[3] in all_frag.keys():
            all_frag[line1[3]].append(line1)
        else:
            all_frag[line1[3]] = [line1]

barcode_pass_filter = []
print("the total barcodes number: %s " % len(all_frag.keys()))

head = f2.readline().strip().split(',')
j = 0
for each in head:
    if each == "passed_filters":
        break
    j+=1
i = 0
for line in f2:
    line = line.strip()
    line1 = line.split(',')
    barcode = line1[0]
    pass_filter = int(line1[j])
    if pass_filter >= args.min_frag and pass_filter <= args.max_frag:
        barcode_pass_filter.append(barcode) 
        i+=1
print("Barcodes number with passed filter: %s" % i)

all_line = []
for each_bar in barcode_pass_filter:
    each_line = all_frag[each_bar]
    all_line.extend(each_line)
all_line2 = pd.DataFrame(all_line)
all_line2 = all_line2.sort_values([0, 1], ascending = (False, True))
all_line2.to_csv(args.output,sep='\t',header=False,index=False)
