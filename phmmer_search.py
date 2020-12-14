def main(input_fasta):
    ###Iteration 1: search through pdb sequences###
    ref_name, ref_seq = read_fasta(input_fasta)
    sys.stderr.write("Imported fasta successfully\nRunning first search on pdb database\n")
    get_hmmer_hits(sequence = ref_seq,
                    seqdb = 'pdb',
                    output_file = 'phmmer_pdb_results.txt')
    sys.stderr.write("Running second search on uniprot database\n")
    get_hmmer_hits(sequence = ref_seq,
                    seqdb = 'uniprotrefprot',
                    output_file = 'phmmer_uniprot_results.txt')
    sys.stderr.write("Running third search on ensembl database\n")
    get_hmmer_hits(sequence = ref_seq,
                    seqdb = 'ensembl',
                    output_file = 'phmmer_ensembl_results.txt')
    sys.stderr.write("Running fourth search on swissprot database\n")
    get_hmmer_hits(sequence = ref_seq,
                    seqdb = 'swissprot',
                    output_file = 'phmmer_swissprot_results.txt')
    sys.stderr.write("Running fifth search on Quest For Orthologs database\n")
    get_hmmer_hits(sequence = ref_seq,
                    seqdb = 'qfo',
                    output_file = 'phmmer_qfo_results.txt')
    sys.stderr.write("Running sixth search on rf15 reference proteome database\n")
    get_hmmer_hits(sequence = ref_seq,
                    seqdb = 'rp15',
                    output_file = 'phmmer_rp15_results.txt')


def read_fasta(fasta_file):
    with open(fasta_file,'r') as fasta:
        read_seq = False
        for line in fasta:
            if line.startswith(">") == True and read_seq == True:
                break
            elif line.startswith(">") == True:
                seq_name = line.strip()[1:]
                read_seq = True
                seq_list = []
            elif read_seq == True:
                seq_list.append(line.strip())
        fasta.close()
    return(seq_name, ''.join(seq_list))

def get_hmmer_hits(sequence, seqdb, output_file, max_hits = 100000):
    params = {'seqdb':seqdb, 'seq':'>Seq\n' + sequence}
    #post the search request to the server
    req = requests.post('https://www.ebi.ac.uk/Tools/hmmer/search/phmmer',data = params)
    results_url = str(req.url)[:-6]
    headers = {'Accept': 'text/xml'}
    params = (('output', 'text'), ('ali', '1'),('range', '1, max_hits'))
    results = requests.get(results_url, headers=headers, params=params)
    with open(output_file, 'w+') as file_out:
        file_out.write(results.text)

if __name__ == "__main__":

    from optparse import OptionParser
    import sys, requests

    parser = OptionParser()
    parser.add_option('--query', action = 'store', type = 'string', dest = 'input_fasta', default = 'query.fasta', help = "sequence from which to construct alignment")
    (option, args) = parser.parse_args()

    main(option.input_fasta)
