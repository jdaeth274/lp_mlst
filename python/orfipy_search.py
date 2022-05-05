import subprocess
import time

def orfipy_search(fasta_loc, num_cores):
    tic_aa_creator = time.perf_counter()
    print("Creating aa")

    ## run orfipy ORF finder on the fasta file

    orfipy_run = "orfipy " + fasta_loc + " --dna tmp_out_dna.fa --outdir tmp_orfi_out --min 150 --procs " + str(num_cores)

    subprocess.run(orfipy_run, shell=True)

    aa_file = "./tmp_orfi_out/tmp_out_dna.fa"

    return aa_file