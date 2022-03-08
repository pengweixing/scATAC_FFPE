#################################################
#  File Name:generate_database.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Mon Jan 17 09:29:47 2022
#################################################
import itertools as it
import sys
i = 1

for each in it.product('AGCT', repeat=16):
    print(''.join(each))
    if i>=int(sys.argv[1]):
        break
    i+=1
