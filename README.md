# MSAgen
Scripts to generate annotated and trimmed multiple sequence alignments from single input sequences

Functions are provided by individual python programs, which can be called by the bash script construct_align.sh
You will need to edit variables near the top of script relevant to your protein. You may also want edit flags in the script to change various parameters that affect stringency of filtering, if the unedited script results in too few (or too many) sequences.
By default, sequence headers will be annotated with phylogenetic information in a way that can be used by pySCA (https://ranganathanlab.gitlab.io/pySCA/). The annotation process provided here is equivalent to, but much faster than, the annotateMSA function in pySCA when running multiple times. However, the first time it is run, it will need to set up a local database for NCBITaxa, which may be slow.

Dependencies:
Python 3
Biopython (https://biopython.org/docs/1.75/api/index.html)
Clustal Omega (http://www.clustal.org/omega/)
ete3 (http://etetoolkit.org)
FastTree (http://meta.microbesonline.org/fasttree/)
