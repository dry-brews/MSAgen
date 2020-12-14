#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 10:50:20 2020

@author: bryan
"""

def main(InputMSA, output):
    ncbi = NCBITaxa()
    #ncbi.update_taxonomy_database()
    headers, seqs = readAlg(InputMSA)
    sys.stdout.write("Annotating headers for %d sequences..." % len(headers))
    
    for i in range(0,len(headers)):     
        head_terms = read_header(headers[i])
        lin = ncbi.get_lineage(head_terms["taxid"])
        #sp_name = ncbi.translate_to_names([tid])
        lin_name = ncbi.translate_to_names(lin)
        with open(output, 'w+') as output_fasta:
            output_fasta.write(">%s|%s|%s\n%s\n" % (head_terms["header"], lin_name[-1], ", ".join(lin_name[1:]), seqs[i]))
    sys.stdout.write("Done\n")
    
def readAlg(filename):
    """
    Read in a multiple sequence alignment in FASTA format, and return the
    headers and sequences.

    headers, sequences = readAlg(filename)
    """

    filelines = open(filename, "r").readlines()
    headers = list()
    sequences = list()
    notfirst = 0
    for line in filelines:
        if line[0] == ">":
            if notfirst > 0:
                sequences.append(seq.replace("\n", "").upper())
            headers.append(line[1:].replace("\n", ""))
            seq = ""
            notfirst = 1
        elif line != "\n":
            seq += line
    sequences.append(seq.replace("\n", "").upper())
    return headers, sequences

def read_header(header):
    terms = {"header":header}
    for t in header.split("&"):
        try:
            terms[t.split("=")[0]] = t.split("=")[1]
        except IndexError:
            sys.stderr.write("Error: failed to parse header: %s" % header)
            sys.exit()
    return terms



if __name__ == "__main__":

    from optparse import OptionParser
    import sys
    from ete3 import NCBITaxa
    parser = OptionParser()
    
    parser.add_option("-i",
                      "--input",
                      action = 'store',
                      dest="Input_MSA",
                      help="input sequence alignment"
                      )
    
    parser.add_option("-o",
                      "--output",
                      action='store',
                      dest="output",
                      default="Output.an",
                      help="Outputfile name. Default: Output.an"
                      )
    (option, args) = parser.parse_args()

    main(option.Input_MSA, option.output)
