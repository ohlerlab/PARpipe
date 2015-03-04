#!/usr/bin/env python
import sys
import math
import re
#usage: reads2bed $yourClusterFileName > $outFileName

def reads2bed(filename):
 file1 = open(filename,'r')
 headerfile1 = file1.readline()
 headerfile1 = headerfile1.strip('\n')
 field = headerfile1.split(',')
 order = [-1, -1, -1, -1, -1, -1, -1, -1]
 for i in range(len(field) - 1):
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
   elif re.search(r'mismatch', field[i]):
     order[5] = i
   elif re.search(r'count', field[i]):
     order[6] = i
   elif re.search(r'end_base', field[i]):
     order[7] = i
 for line in file1:
   line = line.strip('\n')
   field = line.split(',')
   print(field[order[0]]+'\t'+field[order[2]]+'\t'+field[order[3]]+'\t'+field[order[4]]+'\t'+field[order[5]]+','+field[order[6]]+','+field[order[7]]+'\t'+field[order[1]])

 file1.close()

if __name__ == "__main__":
    reads2bed(sys.argv[1])
