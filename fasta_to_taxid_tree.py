#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 11:46:18 2020

@author: bryan
"""

from ete3 import NCBITaxa
import sys, os
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

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

headers, seqs = readAlg(sys.argv[1])
taxids = []
tax_to_head = {}
for i, h in enumerate(headers):
    head_terms = read_header(h)
    tid = int(head_terms["taxid"])
    taxids.append(tid)
    if tax_to_head.get(tid) == None:
        tax_to_head[tid] = []
    tax_to_head[tid].append(head_terms["header"])

for tid in tax_to_head:
    if len(tax_to_head[tid]) == 1:
        tax_to_head[tid] = tax_to_head[tid][0]
    else:
        tax_to_head[tid] = str(tuple(tax_to_head[tid])).replace("'",'')
        #print(tax_to_head[tid])
tree = ncbi.get_topology(taxids, intermediate_nodes = True)
leaves = []
for leaf in tree:
    leaves.append(leaf.name)
#print(leaves)
ancestor = tree.get_common_ancestor(leaves)
tree.set_outgroup('2')
#print(ancestor.name)
#print(tree.leafs())
#for leaf in tree:
#    print(leaf.name)
    #print(tax_to_head.get(int(leaf.name), "none"))   
    #leaf.add_features(name = tax_to_head.get(int(leaf.name), "none"))
#print(tree.get_ascii(attributes = ["name"]))
t = tree.write(format=6)
#tree.write(format=9, outfile="test.nw")
with open("test.nw",'w+') as newick_out:
    newick_out.write(t)   
#print(tax_to_head)   
#print(tax_to_head)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    