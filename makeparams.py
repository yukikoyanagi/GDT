#!/usr/bin/env python3
#
# File: makeparams.py
#
# Time-stamp: <2016-10-24 16:50:24 yuki>
#
# A script to make params.txt file for parameter analysis.
#
# Author: Yuki Koyanagi
# History:
#  


from itertools import combinations

# Make list of avals
a = list(range(4, 13))
avs = [s for l in range(1, 10) for s in combinations(a, l)]
ns = 'nvals=6,8,10'
pc = [';'.join(['avals={}'.format(','.join(map(str, av))), ns])
      for av in avs]
print(len(pc))
# Write to a param file
f = '/home/yuki/QGM/gdt/data/linfit/params.txt'
with open(f, 'w') as outf:
    outf.write('\n'.join(pc))
