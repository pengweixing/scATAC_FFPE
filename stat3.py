#################################################
#  File Name:stat.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 21 Nov 2020 10:49:34 PM UTC
#################################################

import sys
import re
import gzip

f1 = open(sys.argv[1],'r')
path = sys.argv[2]
all_file = [line.strip().split() for line in f1.readlines()]
raw_reads = {}
for each in all_file[0]:
    new = gzip.open(path+'/'+each,'rb')
    each2 = each.split('_L001')
    name=each2[0]
    num = len(new.readlines())
    num = num/4
    raw_reads[name]=[num]

for each in all_file[0]:
    each2 = each.split('_L001')
    myname=each2[0]
    myname2 = path+'/'+'Mapping/'+myname+'_bowtie2.txt'
    map_file = open(myname2,'r')
    all_line = [line.strip() for line in map_file.readlines()]
    after_trim = int(all_line[0].split()[0])
    map_rate = all_line[-1].split()[0]
    map_rate2 = float(map_rate.replace('%',''))
    map_reads = int(map_rate2*after_trim/100)
    
    stat_name = path+'/'+'Mapping/'+myname+'.q2.sort.rmdup.bam.stat'
    stat_file = open(stat_name,'r')
    stat_line = [line.strip() for line in stat_file.readlines()]
    final_reads = stat_line[0].split()[0]
    
    duplicate_name = path + '/' + 'Mapping' + '/' + myname + '.Picard.log'
    duplicate_f = open(duplicate_name,'r')
    all_dupliate_lines = duplicate_f.readlines()
    for each_dup_line in all_dupliate_lines:
        oneLine = each_dup_line.strip()
        words = oneLine.split()
        keyword = 'Marking'
        if keyword in words:
            duplicateCnt = int(words[5])
            dupPercent = duplicateCnt / float(raw_reads[myname][0])
            dupPercent = '{:.2%}'.format(dupPercent)
    raw_reads[myname].extend([after_trim,map_rate,map_reads,dupPercent,final_reads])

    ## add the single cell qc
    sc_name = path+'/'+'sc_Debarcode'+'/'+myname+'/'+'Barcoding_rate_qc.txt'
    sc_file = open(sc_name,'r')
    temp = sc_file.readline()
    temp = sc_file.readline()
    number = 0
    for line in sc_file:
        line = line.strip()
        line1 = line.split()
        number += int(line1[1])
    raw_reads[myname].extend([number])
    sc_rate = number/float(raw_reads[myname][0])   
    sc_rate = '{:.2%}'.format(sc_rate)
    raw_reads[myname].extend([sc_rate])

print('##final reads: removed chrM and with Q2 filtering\n')
print('name\tall_reads\tafter_trimmed\tmapping_rate\tmapped_reads\tdupPercent\tfinal\tsc_reads\tsc_rate')
for key,value in raw_reads.items():
    print(key,*value,sep="\t")

