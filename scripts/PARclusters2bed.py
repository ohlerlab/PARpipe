#!/usr/bin/env python
import sys
import re
#usage: PARclusters2bed $yourClusterFileName > $outFileName

headerfile1 = sys.stdin.readline()
headerfile1 = headerfile1.strip('\n')
field = headerfile1.split(',')
order = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
for i in range(len(field)):
  if re.search(r'^Chr', field[i]):
    order[0] = i
  elif re.search(r'Strand', field[i]):
    order[1] = i
  elif re.search(r'Start', field[i]):
    order[2] = i
  elif re.search(r'End', field[i]):
    order[3] = i
  elif re.search(r'(Aligned to)|(TranscriptLocation)', field[i]):
    order[4] = i
  elif re.search(r'(GeneName)|(gene name)', field[i]):
    order[5] = i
  elif re.search(r'ID$', field[i]):
    order[6] = i
  elif re.search(r'Sequence$', field[i]):
    order[7] = i
  elif re.search(r'ReadCount', field[i]):
    order[8] = i
  elif re.search(r'ModeLocation', field[i]):
    order[9] = i
  elif re.search(r'ModeScore', field[i]):
    order[10] = i
  elif re.search(r'ConversionLocationCount', field[i]):
    order[11] = i
  elif re.search(r'^ConversionEventCount', field[i]):
    order[12] = i
  elif re.search(r'NonConversionEventCount', field[i]):
    order[13] = i
  elif re.search(r'AnnotationSource', field[i]):
    order[14] = i
vli = []
for i in range(4, 15):
  if order[i] >= 0:
    vli.append(order[i])
for line in sys.stdin:
  line = line.strip('\n')
  field = line.split(',')
  
#  if (location =='3\'utr' or location =='5\'utr' or location =='coding' or location == 'intron'):
#   print(chrom+'\t'+start+'\t'+end+'\t'+clusterId+'\t'+readCount+'\t'+strand)
#   print(chrom+'\t'+start+'\t'+end+'\t'+clusterId+'\t'+conversions+'\t'+strand)
  prli = []
  for i in range(len(vli)):
    prli.append(field[vli[i]])
  print(field[order[0]]+'\t'+field[order[2]]+'\t'+field[order[3]]+'\t'+field[order[6]]+'\t'+",".join(prli)+'\t'+field[order[1]])
