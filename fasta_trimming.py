#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 19:52:56 2020

@author: bryan
"""

def main(input_fasta, output, qname, min_cov):
    with open(input_fasta,'r') as fasta_in:
         query_list = []
         for line in fasta_in:
              if line.strip() == ">"+qname:
                   read_seq = True
              elif read_seq == True and line.startswith(">"):
                   break
              elif read_seq == True:
                   query_list.append(line.strip())
         query_seq = ''.join(query_list)
    
    min_seq_len = float(min_cov)*len(query_seq)
    
    seq_dict = {}
    dups_removed = 0
    too_short = 0
    seqs_written = 0
    
    with open(input_fasta,'r') as fasta_in:
         with open(output,'w+') as fasta_out:
              for line in fasta_in:
                   if line.startswith('>'):
                        seq_name = line.strip()
                   else:
                        seq = line.strip()
                        if seq_dict.get(seq) == None:
                             if len(seq) < min_seq_len:
                                  too_short += 1
                                  continue
                             else:
                                  fasta_out.write(seq_name + '\n')
                                  fasta_out.write(seq + '\n')
                                  seq_dict[seq] = seq_name
                                  seqs_written +=1             
                        else:
                             dups_removed +=1
    sys.stdout.write("%d duplicate sequences removed\n" % dups_removed)
    sys.stdout.write("%d short sequences removed\n" % too_short)
    sys.stdout.write("(min length = %d)\n" % int(min_seq_len))
    sys.stdout.write("%d sequences written to deduplicated_seqs.fasta\n" % seqs_written)

if __name__ == "__main__":

    import argparse
    import sys
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i",
                  "--input",
                  action = 'store',
                  dest="input_fasta",
                  help="input sequence alignment"
                  )
    parser.add_argument("-o",
                  "--output",
                  action='store',
                  dest="output",
                  default="Output.an",
                  help="Outputfile name. Default: Output.an"
                  )
    parser.add_argument("-t",
                  "--query_title",
                  action='store',
                  dest="qname",
                  default="NA",
                  help="Fasta containing only the query sequence"
                  )
    parser.add_argument("-m",
                  "--min_coverage",
                  action='store',
                  dest="min_cov",
                  default="0.8",
                  help="Fasta containing only the query sequence"
                  )

    options = parser.parse_args()

    main(options.input_fasta, options.output, options.qname, options.min_cov)