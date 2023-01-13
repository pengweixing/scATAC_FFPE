#################################################
#  File Name:get_singcell_qc.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Fri Nov 25 15:45:09 2022
#################################################


import pysam
import sys
import argparse
import copy
import math 

def fargv():
    parser = argparse.ArgumentParser(usage="python ")
    parser.add_argument('-b1',"--bam1",help="the sorted bam file ", required=True)
    parser.add_argument('-b2',"--bam2",help="the bam file within the peaks ", required=True)
    parser.add_argument('-b3',"--bam3",help="the bam file within the tss ", required=True)
    parser.add_argument('-o',"--output",help="the QC matrix file ", required=True)
    parser.add_argument('-q',"--mapQuality",help="the q value of mapping quality ", default=0,type=int)
    args = parser.parse_args()
    return args

def get_qc(f_bam1,f_bam2,f_bam3,f_output,mapq):

    f_out = f_output
    cell_stats = {}
    temp = {'Total':0,'Duplicates':0,'LowQ':0,'Unpaired':0,'Final':0}
    for a in f_bam1:
        if not a.is_unmapped:
            query_name = a.query_name
            cell_name = "_".join([query_name.split('_')[0]]+query_name.split('_')[5:])
            #cell_name = '_'.join([query_name.split('_')[i] for i in [0,5]])
            if not cell_name in cell_stats:
                cell_stats[cell_name] = copy.deepcopy(temp)
                cell_stats[cell_name]['Total'] +=1
            else:
                cell_stats[cell_name]['Total'] += 1
            if a.is_duplicate:
                cell_stats[cell_name]['Duplicates'] += 1
            else:
                if a.mapping_quality < mapq:
                    cell_stats[cell_name]['LowQ'] += 1
                else:
                    if not a.is_proper_pair:
                        cell_stats[cell_name]['Unpaired'] += 1
                    else:
                        cell_stats[cell_name]['Final'] += 1
    number_FRiP = {}
    for b in f_bam2:
        query_name = b.query_name
        cell_name = "_".join([query_name.split('_')[0]]+query_name.split('_')[5:])
        if cell_name in number_FRiP:
            number_FRiP[cell_name] += 1
        else:
            number_FRiP[cell_name] = 1
            
    number_FRiT = {}
    for b in f_bam3:
        query_name = b.query_name
        cell_name = "_".join([query_name.split('_')[0]]+query_name.split('_')[5:])
        if cell_name in number_FRiT:
            number_FRiT[cell_name] += 1
        else:
            number_FRiT[cell_name] = 1

    f_out.write('barcode\ttotal\tDups\tLowQ\tUnpaired\tFinal\tFRiP\tFRiT\n')             
    for i in cell_stats:
        total = math.ceil(cell_stats[i]['Total']/2)
        Duplicates = round(cell_stats[i]['Duplicates']/2/total,4)*100
        LowQ = int(cell_stats[i]['LowQ'])
        Unpaired = int(cell_stats[i]['Unpaired']/2)
        Final = math.ceil(cell_stats[i]['Final']/2)
        if i in number_FRiP:
            frags = number_FRiP[i]/2
        else:
            frags = 0
        
        if i in number_FRiT:
            fragst = number_FRiT[i]/2
        else:
            fragst = 0
            
        if Final>0:
            FRiP = frags/Final*100
            FRiT = fragst/Final*100
        else:
            FRiP = 0
            FRiT = 0
        out = "{cellname}\t{total}\t{Duplicates:.2f}\t{LowQ}\t{Unpaired}\t{Final}\t{FRiP:.2f}\t{FRiT:.2f}\n"
        f_out.write(out.format(cellname=i,total=total,Duplicates=Duplicates,LowQ=LowQ,Unpaired=Unpaired,Final=Final,FRiP=FRiP,FRiT=FRiT))

    f_bam1.close()
    f_out.close()
    f_bam2.close()
    f_bam3.close()
    
def main(kwargs):
    f_bam1 = pysam.AlignmentFile(kwargs.bam1,mode='rb',check_header=True,threads=3)
    f_bam2 = pysam.AlignmentFile(kwargs.bam2,mode='rb',check_header=True,threads=3)
    f_bam3 = pysam.AlignmentFile(kwargs.bam3,mode='rb',check_header=True,threads=3)
    
    f_output = open(kwargs.output,'w')
    mapq = int(kwargs.mapQuality)
    get_qc(f_bam1,f_bam2,f_bam3,f_output,mapq)
    
if __name__ == "__main__":
    kwargs = fargv()
    main(kwargs)
