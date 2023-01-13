#################################################
#  File Name:extract_bam_from_SCname.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Wed Aug 10 10:09:54 2022
#################################################


import sys
import pysam

f1 = open(sys.argv[1],'r') ##name
f_bam_input = pysam.AlignmentFile(sys.argv[2],mode='rb',check_header=True)
f_bam_output = pysam.AlignmentFile(sys.argv[3],mode='wb',template=f_bam_input)

all_cell = [line.strip() for line in f1.readlines()]
all_cell = set(all_cell)

for each in f_bam_input:
    c = each.get_tags()
    sc_name = dict(c)['CB']
    if sc_name in all_cell:
        f_bam_output.write(each)

f1.close()
f_bam_input.close()
f_bam_output.close()