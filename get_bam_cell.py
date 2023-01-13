#################################################
#  File Name:get_bam_cell.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Wed Feb  2 16:27:40 2022
#################################################

import pysam
import sys
import gzip

f_bed = gzip.open(sys.argv[1],'rt')

all_name = set([line.strip().split()[3] for line in f_bed.readlines()])

samfile = pysam.AlignmentFile(sys.argv[2], "rb")
out = sys.argv[2].replace('bam.bam','filtered.bam')
f_out = pysam.AlignmentFile(out,mode='wb',template=samfile)

for read in samfile:
    cell = read.get_tag("CB")
    if not read.is_duplicate:
        if cell in all_name:
            f_out.write(read)

f_out.close()
