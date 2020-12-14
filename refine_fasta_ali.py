#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 12:59:06 2020

@author: bryan
"""

def main(input_fasta, output_fasta, exclude_fasta, query_name, maxGap, minID, maxID):

    #First, read the aligned fasta into memory
    seqs, seqs_added = read_input_fasta(input_fasta)
    sys.stderr.write("Read %d seqs from aligned fasta file\n" % seqs_added)

    #Second, read the names of sequences that you want to exclude (if none excluded, produce an empty list)
    if exclude_fasta != "NULL":
        exclude_seqs, exclude_added = read_input_fasta(exclude_fasta)
        exclude_headers = list(exclude_seqs.keys())
    else:
        exclude_seqs  = {}
        exclude_headers = []

    #Third, filter out sequences by the number of gaps they contain (bypass those in exclude_seqs)
    seqs, seqs_cov_filt = filter_by_coverage(seqs, max_gaps = maxGap, e_dict = exclude_seqs)
    sys.stderr.write("%d sequences removed for having too many gaps\n" % seqs_cov_filt)
    
    #Fourth, filter out sequences by %identity relative to the query (bypass those in exclude_seqs)
    seqs, seqs_id_filt = filter_by_pctid(seqs, query_name, min_id = minID, max_id = maxID, e_dict = exclude_seqs)
    sys.stderr.write("%d sequences removed for being to dissimilar or too similar to the query\n" % seqs_id_filt)
    
    #Fifth, if there are sequences to exclude, check if any sequences better match (%id) the exluded sequence than the query, and remove them if they do
    #This process will remove the excluded sequences themselves, of course
    if exclude_fasta != "NULL":
        sys.stderr.write("Comparing sequences against %d excluded sequences\n" % exclude_added)
        seqs = filter_by_undesired_match(seqs, query_name, exclude_headers)
    
    #Finally, print the sequences that passed all filters to a new fasta
    print_refined_ali(seqs, outfile = output_fasta)

def read_input_fasta(fasta_file):
     with open(fasta_file, 'r') as fasta_in:
          seqs_out = {}
          seqs_added = 0
          for line in fasta_in:
               if line.startswith(">"):
                    curr_seq = line.strip(">").strip()
                    seqs_out[curr_seq] = ''
                    seqs_added +=1
               else:
                    seqs_out[curr_seq] += line.strip()
          fasta_in.close()
     return (seqs_out, seqs_added)

def filter_by_coverage(seqs, max_gaps, e_dict = {}):
     seqs_out = {}
     seqs_removed = 0
     for name in seqs:
          if seqs[name].count("-") < (max_gaps * float(len(seqs[name]))) or e_dict.get(name) != None:
               seqs_out[name] = seqs[name]
          else:
               seqs_removed +=1
     return (seqs_out, seqs_removed)

def filter_by_pctid(seqs, query_name, min_id, max_id, e_dict = {}):
    with open("pct_ids_temp.tsv","w+") as pct_id_file: 
        seqs_out = {query_name: seqs[query_name]}
        seq_len = len(seqs[query_name])
        seqs_removed = 0
        for name in seqs:
             IDs = 0
             for i in range(0, seq_len):
                  if seqs[name][i] == seqs[query_name][i]:
                       IDs +=1
             pctid = float(IDs) / float(seq_len)
             pct_id_file.write(str(pctid)+'\n')
             if (pctid > min_id and pctid < max_id) or e_dict.get(name) != None:
                  seqs_out[name] = seqs[name]
             else:
                  seqs_removed +=1
        return (seqs_out, seqs_removed)
               
def filter_by_undesired_match(seqs, query_name, exclude_headers):
     seqs_out = {}
     seq_marks = {}
     for name in seqs:
          seq_marks[name] = True
     seq_len = len(seqs[query_name])
     for e_name in exclude_headers:
          seqs_excluded = 0
          for name in seqs:
               IDs_query = 0
               IDs_exclude = 0
               for i in range(0, seq_len):
                    if seqs[name][i] == seqs[query_name][i]:
                         IDs_query +=1
                    if seqs[name][i] == seqs[e_name][i]:
                         IDs_exclude +=1
               if IDs_query < IDs_exclude:
                    seq_marks[name] = False
                    seqs_excluded +=1
          sys.stderr.write("%d sequences flagged for removal for better matching " %seqs_excluded + e_name + " than the query\n")
     for name in seq_marks:
          if seq_marks[name] == True:
               seqs_out[name] = seqs[name]
     return seqs_out

def print_refined_ali(seqs, outfile = "refined_ali.fasta"):
    with open(outfile,'w+') as fasta:
          for name in seqs:
               fasta.write(">" + name + '\n')
               fasta.write(seqs[name]+'\n')
    sys.stderr.write("Wrote %d sequences to file\n" % len(seqs))

if __name__ == "__main__":

    from optparse import OptionParser
    import sys

    parser = OptionParser()
    parser.add_option('--in', action = 'store', type = 'string', dest = 'input_fasta', default = 'query.fasta', help = "sequence from which to construct alignment")
    parser.add_option('--out', action = 'store', type = 'string', dest = 'output_fasta', default = 'reined_ali.fasta', help = "sequence to write to")
    parser.add_option('--exclude', action = 'store', type = 'string', dest = 'exclude_fasta', default = 'NULL', help = "tsv list of undesired sequences. Sequences that match better to an undesired sequence than the query will be discarded")
    parser.add_option('--query', action = 'store', type = 'string', dest = 'query_name', default = 'Query', help = "name of the query sequence (no >)")
    parser.add_option('--maxGap', action = 'store', type = 'float', dest = 'maxGap', default = 0.25, help = "Maximum fraction of positions that can be gaps in retained sequences")
    parser.add_option('--minID', action = 'store', type = 'float', dest = 'minID', default = 0.15, help = "Minimum %identical positions relative to the query in retained sequences")
    parser.add_option('--maxID', action = 'store', type = 'float', dest = 'maxID', default = 0.99, help = "Maximum %identical positions relative to the query in retained sequences")
    (option, args) = parser.parse_args()

    main(option.input_fasta, option.output_fasta, option.exclude_fasta, option.query_name, option.maxGap, option.minID, option.maxID)