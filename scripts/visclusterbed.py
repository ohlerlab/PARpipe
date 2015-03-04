#!/usr/bin/env python
import re
import sys
#usage: PARclusters2bed $yourClusterFileName > $outFileName

headerfile1 = sys.stdin.readline()
headerfile1 = headerfile1.strip('\n')
field = headerfile1.split(',')
order = [-1, -1, -1, -1, -1, -1]
for i in range(len(field)):
  if re.search(r'^Chr', field[i]):
    order[0] = i
  elif re.search(r'Strand', field[i]):
    order[1] = i
  elif re.search(r'Start', field[i]):
    order[2] = i
  elif re.search(r'End', field[i]):
    order[3] = i
  elif re.search(r'ID$', field[i]):
    order[4] = i
  elif re.search(r'T2Cfraction', field[i]):
    order[5] = i
for line in sys.stdin:
  line = line.strip('\n')
  field = line.split(',')

  print(field[order[0]]+'\t'+str(int(field[order[2]])-1)+'\t'+field[order[3]]+'\t'+field[order[4]]+'\t'+field[order[5]]+'\t'+field[order[1]])
