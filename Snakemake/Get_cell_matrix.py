#################################################
#  File Name:Get_cell_matrix.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Thu Dec  1 15:07:50 2022
#################################################


import pysam
import sys
import argparse

def fargv():
    parser = argparse.ArgumentParser(usage="python ")
    parser.add_argument('-p',"--peaks",help="the peaks file ", required=True)
    parser.add_argument('-n',"--cellnames",help="the selected cellnames file ", required=True)
    parser.add_argument('-ib',"--bam",help="the singlecell bam file ", required=True)
    parser.add_argument('-o',"--out",help="the output of matrix ",required=True)
    args = parser.parse_args()
    return args

def main(kwargs):

    f_peaks = open(kwargs.peaks,'r')
    f_names = open(kwargs.cellnames,'r')
    f_bamin = pysam.AlignmentFile(kwargs.bam,mode='rb',check_header=True,threads=3)
    
    f_out = open(kwargs.out,'w')
    allline = [line.strip() for line in f_names.readlines()]
    allcell_names = {line.strip():0 for line in allline}
    all_cell = allcell_names.keys()
    all_cell2 = '\t'.join(all_cell)
    all_cell3 = "barcode\t"+all_cell2
    f_out.write("%s\n" % all_cell3)
        
    for line in f_peaks:
        line = line.strip()
        line1 = line.split()
        achr = line1[0]
        start = int(line1[1])
        end = int(line1[2])
        peak = achr+"-"+str(start)+"-"+str(end)
        allcell_names = {line.strip():0 for line in allline}
        for a in f_bamin.fetch(achr, start, end):
            cellname = a.get_tag('RG')
            allcell_names[cellname] += 1
        values = []
        for key in allcell_names:
            value = allcell_names[key]
            values.append(str(value))
        values = '\t'.join(values)
        f_out.write("%s\t%s\n" % (peak,values))
        
if __name__ == "__main__":
    kwargs = fargv()
    main(kwargs)
