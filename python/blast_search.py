import time
import subprocess
import os
import sys
from Bio import SearchIO

class seqqer():
    def __init__(self):
        self.seq = None


def hmm_search_for_gene(fasta, gene, aa_dir_name, data_dir, gene_leng, tolerance, hmm_base):

    toc_aa_creator = time.perf_counter()
    print("Running HMM on aa")
    tic_hmm_run = time.perf_counter()
    # run HMM
    #gff_base = re.sub("\..*[a-z,A-Z,0-9].*$", "", gff_base)
    aa_base_name = aa_dir_name + "/" +  hmm_base + "_" + gene
    hmm_output_fn = aa_base_name
    print(hmm_output_fn)
    hmm_proc_out = 1
    while hmm_proc_out != 0:
        hmm_proc_out = subprocess.check_call('hmmsearch ' + data_dir  + gene + '.hmm ' + fasta + ' > ' + hmm_output_fn, shell=True)

    #subprocess.run("rm -r tmp_orfi_out", shell=True)

    # parse HMM
    print("Parsing HMM")
    hmm_output = list(SearchIO.parse(hmm_output_fn, 'hmmer3-text'))
    hmm_best_hsp_index = None
    hmm_best_hsp_bitscore = 0
    for i, hsp in enumerate(hmm_output[0]):
        if hsp.bitscore > hmm_best_hsp_bitscore:
            hmm_best_hsp_bitscore = hsp.bitscore
            hmm_best_hsp_index = i

    if hmm_best_hsp_index is None:
        return "Missing_HMM", hmm_output_fn

    hmm_best_hsp = hmm_output[0][hmm_best_hsp_index][0]
    #subprocess.call('rm ' + hmm_output_fn, shell=True)

    toc_hmm_run = time.perf_counter()
    print("Printing out the resistance")
    seq = None
    print(hmm_best_hsp)
    aa_seq = hmm_best_hsp.hit
    if len(aa_seq) < gene_leng - tolerance:
        print("Gene length not long enough has %s needs %s" % (len(aa_seq), gene_leng))
    else:
        seq = aa_seq
    status = None

    print("Took this long for %s HMM run: %s" % (gene, toc_hmm_run - tic_hmm_run))

    return seq


def blast_search_for_gene(fasta, gene, data_dir, search_dir, iso_name, threads):
    ## Blast search for particular MLST gene
    db_file = data_dir + "blast-db/" + gene + "_blastdb"
    out_csv = search_dir + iso_name + "_" + gene + ".csv"
    blast_cmd = "blastn -query " + fasta + " -db " + db_file + " -outfmt 10 -out " + out_csv + " -num_threads " + str(threads)

    try:
        subprocess.check_call(blast_cmd, shell=True)
    except subprocess.SubprocessError:
        sys.exit("Failed while running blast on isolate: %s for gene %s" % (iso_name, gene))

    return out_csv

def get_aln_pos_from_ref(hmm_aln, pos, offset):
    pos = pos - offset
    ref_pos = pos + 1
    upstream_length = 0
    while_counter = 0
    while upstream_length < pos:
        old_upstream = upstream_length
        upstream_frag = hmm_aln[0, 0:ref_pos].seq

        upstream_length = len(upstream_frag) - upstream_frag.count('.')
        ref_pos = ref_pos + 1
        while_counter += 1
        if while_counter > 10:
            if (upstream_length - old_upstream) == 0:
                return False

    return (ref_pos - 1)