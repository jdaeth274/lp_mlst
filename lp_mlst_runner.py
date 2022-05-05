import argparse
from python.common import main

def parse_input():
    ## Just take in the fasta files and then the output name
    ## use the default path for the data

    purpose = ''' This is a scipt to use the MLST profiles from momPS https://www.sciencedirect.com/science/article/pii/S1198743X17300071 \n
    https://doi.org/10.1016/j.cmi.2017.01.002 \n
    and define these for a set of input 
    Legionella pneumophila isolates.\n Usage: \n
    python lp_mlst_runner.py --seqs <list_of_fastas>  --output <output_prefixes> --threads <num_cores_to_use>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='lp_mlst_runner.py')

    parser.add_argument('--seqs', required=True, help='List of seqeuence files (FASTA) (required)', type=str)
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)
    parser.add_argument('--threads', default=1, help='Number of threads to use for ORF finder', type = int)

    args = parser.parse_args()

    return args

def main_run():
    args = parse_input()
    main(args)


if __name__ == '__main__':
    main_run()

