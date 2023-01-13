#################################################
#  File Name:get_bulk_qc.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Thu Dec  1 16:08:42 2022
#################################################


import sys
outdir=sys.argv[1]

f_out = open(outdir+'/All.QC.txt','w')


with open(outdir+'/sample.group','r') as f:
    sampledict = {}
    all_line = [line.strip() for line in f.readlines()]
    for a in all_line:
        group = a.split()[0]
        if '_' in group:
            group = group.replace('_','.')
        tn5 = a.split()[1]
        if not group in sampledict:
            sampledict[group] = [tn5]
        else:
            sampledict[group].append(tn5)
            
f_barcoding = open(outdir+'/01.debarcoding/Barcoding_rate_qc.txt','r')
header = f_barcoding.readline().strip()
total = header.split()[1]
a=f_barcoding.readline()
all_rest = [line.strip() for line in f_barcoding.readlines()]

total_decoded = 0
f_out.write('\n#### Decoding Info ####\n')
f_out.write('%s\n' % header)
for group in sampledict:
    tn5s = sampledict[group]
    g_number = 0
    for tn5 in tn5s:
        for b in all_rest:
            b = b.split()
            btn5 = b[0]
            number = int(b[1])
            if btn5 == tn5:
                g_number += number
    total_decoded += g_number
    f_out.write('%s\t%s\n' % (group,g_number))

pct = round(total_decoded/int(total),3)
f_out.write('Total decoded\tDecoding rate\n')
f_out.write('%s\t%s\n' % (total_decoded,pct))
f_out.write('\n#### Mapping Info ####\n')
f_out.write('Name\tTrimmedReads\tMapping_rate\tProper_paired\tDup_Rate\tFinal\n')
for group in sampledict:
    with open(outdir+'/02.mapping/' + group + '.sort.markdup.merge.stat','r') as f:
        i = 0
        for c in f:
            c = c.strip()
            if i == 0:
                trimmed_number = int(int(c.split()[0])/2)
            elif i == 3:
                dups = int(c.split()[0])
            elif i == 4:
                map_rate = c.split('(')[1].split(':')[0]
            elif i == 8:
                proper = c.split('(')[1].split(':')[0]
            i += 1
        dup_rate = round(dups/2/trimmed_number,2)
        dup_rate2 = str(dup_rate*100)+"%"
    with open(outdir+'/02.mapping/' + group + '.sort.rmdup.merge.pair.bam.flagstat','r') as f1:
        for d in f1:
            d = d.strip()
            if 'read1' in d:
                number1 = int(d.split()[0])
            elif 'read2' in d:
                number2 = int(d.split()[0])
        number = max([number1,number2])
        f_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (group,trimmed_number,map_rate,proper,dup_rate2,number))

            
