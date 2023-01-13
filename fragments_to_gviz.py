#################################################
#  File Name:fragments_to_gviz.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Wed Sep 28 13:32:35 2022
#################################################
#cell_type_name.txt
"""
cell	type
spleen_FFPE#AAAAAACGCAGGGCGC-9	other
spleen_FFPE#AAAAAAGTCTGTAATA-7	other
spleen_FFPE#AAAAAAGTCTGTAGAT-4	other
spleen_FFPE#AAAAAACTCTGTCAGT-12	other
spleen_FFPE#AAAAAAGATTAGTGTG-10	other
"""

import gzip,sys
f1 = gzip.open(sys.argv[1],'rb') ##fragments.filter.sort.tsv.gz
f2 = open(sys.argv[2],'r')##cell_type_name.txt

f2.readline()
mydict = {}
allline = [line.strip() for line in f2.readlines()]
alltype = ['T_cells','B_cells','other']
    
tdict = {}
bdict = {}
odict = {}
t = 1
b = 1
o = 1
for line in allline:
    line1 =line.split()
    cell = line1[0]
    mytype = line1[1]
    cell2 = cell.split('#')[1]
    if mytype == 'T_cells':
        tdict[cell2] = [mytype,t]
        t += 1
    elif mytype == 'B_cells':
        bdict[cell2] = [mytype,b]
        b += 1
    elif  mytype == 'other':
        odict[cell2] = [mytype,o]
        o += 1

for line in f1:
    line = line.decode().strip()
    line1 = line.split()
    bc = line1[3]
    if bc in tdict:
        aa_type = tdict[bc][0]
        aa_id = tdict[bc][1]
        print(line1[0],line1[1],line1[2],aa_id,aa_type,sep="\t")

    elif bc in bdict:
        aa_type = bdict[bc][0]
        aa_id = bdict[bc][1]
        print(line1[0],line1[1],line1[2],aa_id,aa_type,sep="\t")

    elif bc in odict:
        aa_type = odict[bc][0]
        aa_id = odict[bc][1]
        print(line1[0],line1[1],line1[2],aa_id,aa_type,sep="\t")


f1.close()
f2.close()