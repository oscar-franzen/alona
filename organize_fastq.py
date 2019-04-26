#!/usr/bin/env python3

"""
 Organizes barcodes and UMIs from a FASTQ file. 
 
 How to use:
 python3 organize_fastq.py -i sequences.fastq -b barcodes.bc -c 1000
 
 Details:
 https://github.com/oscar-franzen/alona/blob/master/preprocess_single_cell_rnaseq_data.md
 
 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import getopt, sys, re

target_bc = {}

def read_bc_file(input_bc_file, freq_min):
  with open(input_bc_file,'r') as f:
    for line in f:
      bc, freq = line.rstrip('\n').split(' ')
      
      if int(freq) < int(freq_min):
        continue
        
      if re.search('^A+$',line) or \
         re.search('^T+$',line) or \
         re.search('^G+$',line) or \
         re.search('^C+$',line):
        continue
      
      target_bc[bc] = 1

def prep_fq_file(input_fq_file):
  with open(input_fq_file,'r') as f:
    for line in f:
      line1 = line
      line2 = next(f)
      line3 = next(f)
      line4 = next(f)
      
      bc = line2[0:12]
      umi = line2[12:12+8]
      seq = line2[12+8:]
      quals = line4[12+8:]
      
      if bc in target_bc and seq.count('N') < 10:
        new_header = '%s_%s_%s' % (line1.split(' ')[0],bc,umi)
        print(new_header)
        print(seq)
        print('+%s' % new_header[1:])
        print(quals)
        
        d = sys.stdin.readline()

def main():
  try:
    args, values = getopt.getopt(sys.argv[1:], "i:b:c:", ["input=","barcode=","count"])
  except getopt.GetoptError as err:
    print(err)
    sys.exit(2)
  
  input_fq_file = ''
  input_bc_file = ''
  freq_min = 1000
  
  for o, a in args:
    if o == '-i':
      input_fq_file = a
    elif o == '-b':
      input_bc_file = a
    elif o == '-c':
      freq_min = a
      
  if input_fq_file == '' or input_bc_file == '':
    print("""Usage:
               -i <file>\tInput FASTQ file.
               -b <file>\tInput barcode file.
               -c <freq>\tMinimum number of detected barcodes.""")
    sys.exit(2)
  
  read_bc_file(input_bc_file,freq_min)
  prep_fq_file(input_fq_file)

if __name__ == "__main__":
  main()
