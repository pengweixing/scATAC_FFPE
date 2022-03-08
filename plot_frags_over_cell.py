#################################################
#  File Name:temp.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Wed Dec  1 22:03:23 2021
#################################################

import sys
from fastq_to_10x import Fastq_transform
HELP = """\npython plot_frags_over_cell.py Barcodes_fragments_qc.txt\n"""

def main():
    f = open(sys.argv[1],'r')
    f.readline()
    allline = [line.strip().split() for line in f.readlines()]
    debarcoding_obj = Fastq_transform()
    debarcoding_obj.hist_plot(allline)

if sys.argv[-1].endswith(".py"):
    print(HELP)
else:
    main()