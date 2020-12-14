#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 20:24:33 2020

@author: bryan
"""

import sys

query_name = sys.argv[2]

with open(sys.argv[1],'r') as ali_fasta:
    #sys.stderr.write("got here\n")
    #sys.stderr.write(query_name + '\n')
    read_seq = False
    for line in ali_fasta:
          if line.strip().strip(">") == query_name:
               read_seq = True
               seq_list = []
               #sys.stderr.write("found query")
          elif read_seq == True:
               seq_list.append(line.strip())
          elif line.startswith(">"):
              #print(line)
              read_seq = False
    ali_fasta.close()
query_seq = ''.join(seq_list)

gap_dict = {}
for i in range(0,len(query_seq)):
     if query_seq[i] in ["-","."," "]:
          gap_dict[i] = True
     else:
          gap_dict[i] = False
          
with open(sys.argv[1],'r') as ali_fasta:
     with open("trunc_ali.fasta","w+") as fasta_out:
          seq_line = ["NULL"]
          for line in ali_fasta:
               if line.startswith(">"):
                    if seq_line != ["NULL"]:
                         fasta_out.write(''.join(seq_line)+'\n')
                    fasta_out.write(line.strip()+'\n')
                    seq_line = []
                    i = 0
               else:
                    for letter in line.strip():
                         if gap_dict[i] == False:
                              seq_line.append(letter)
                         i+=1
          fasta_out.write(''.join(seq_line)+'\n')

               