import math
import pandas as pd


def process_blast(blast_csv, gene_length, gene):
    ## Import the blast csv and get the top hit from the alignment

    blast_cols = ['query','subject','pid','length','gap','mismatch',
                  'qstart','qend','sstart','send','eval','bitscore']
    gene_len = len(gene_length)
    blast_res = pd.read_csv(blast_csv, names=blast_cols, header = None)
    if blast_res.empty:
        print("No blast matches for this %s in isolate" % (gene))
        return math.nan
    else:
        hundred_matches = blast_res[blast_res['pid'] == 100]
        if hundred_matches.empty:
            print("No perfect matches for this isolate for gene: %s"  % gene)
            return int(blast_res.iloc[0,1])
        else:
            for k, length in enumerate(gene_length):
                perf_length = hundred_matches[hundred_matches['length'] == length + 1]
                if perf_length.empty:
                    perf_length = hundred_matches[hundred_matches['length'] == length]
                    if perf_length.empty:
                        if k == (gene_len - 1):
                            print("No perfect length matches for this isolate for gene: %s" % gene)
                            return int(hundred_matches.iloc[0, 1])
                        else:
                            continue
                    else:
                        return int(hundred_matches.iloc[0,1])
                else:
                    return int(hundred_matches.iloc[0, 1])


def ST_process(data_dir, isolate_row):
    ## left join the current isolates row to get the ST
    pub_mlst_sheet = data_dir + "lp_profiles.tsv"
    pub_mlst = pd.read_table(pub_mlst_sheet, sep="\t", header=0)

    updated_iso = isolate_row.merge(pub_mlst, how="left")

    return updated_iso
