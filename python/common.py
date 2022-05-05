import re
import sys
import pandas
import time
import os

import pandas as pd

from python import orfipy_search
from python import blast_search
from python import blast_processing


def main(input_args):
    tic_setup = time.perf_counter()
    print("Beginning setup")
    pandas.set_option('display.max_columns', 500)

    # parse command line


    ###############################################################################
    ## Lets go through the gff list in a for loop, we'll find the position of the #
    ## gene of interest and then from there extract this sequence and compare it ##
    ## to the gene alignment to test which is the right gene ######################
    ###############################################################################
    seq_files = open(input_args.seqs, "r")
    seq_lines = seq_files.read().splitlines()

    python_dir_name = os.path.dirname(os.path.realpath(__file__))
    data_dir = re.sub("python","data/", python_dir_name)
    print(data_dir)
    aa_dir_name = "./" +  input_args.output + "_aa_dir" + "/"
    if not os.path.exists(aa_dir_name):
        os.mkdir(aa_dir_name)

    df_names = ['id','ST','flaA','pilE','asd','mip','mompS','proA','neuA']
    out_df = pd.DataFrame(index=range(len(seq_lines)), columns =df_names)

    skip = False
    toc_setup = time.perf_counter()
    print("Took this long for initial set up: %s" % (toc_setup - tic_setup))
    print("Beginning iso run")
    for k, fasta_file in enumerate(seq_lines):
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("On isolate %s of %s" % (k + 1, len(seq_lines)))
        out_st = table_row_getter(fasta_file, data_dir, aa_dir_name, 1)
        out_df.loc[k] = out_st.iloc[0]
        # tic_iso_run = time.perf_counter()
        # bassio_nameo = os.path.basename(fasta_file)
        # print(bassio_nameo)

        # ## Perform the ORF finder for the isolate

        # #ORF_loc = orfipy_search.orfipy_search(fasta_file, input_args.threads)

        # ## Search for cpn60
        # flaA_seq = blast_search.blast_search_for_gene(fasta_file, "flaA", data_dir, aa_dir_name, bassio_nameo, input_args.threads)

        # # Search for fusA
        # pilE_seq = blast_search.blast_search_for_gene(fasta_file, "pilE", data_dir, aa_dir_name, bassio_nameo, input_args.threads)

        # # Search for gltA
        # asd_seq = blast_search.blast_search_for_gene(fasta_file, "asd", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # # Search for pyrG
        # mip_seq = blast_search.blast_search_for_gene(fasta_file, "mip", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # # Search for recA
        # mompS_seq = blast_search.blast_search_for_gene(fasta_file, "mompS", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # # Search for rplB
        # proA_seq = blast_search.blast_search_for_gene(fasta_file, "proA", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # # Search for rpoB
        # neuA_seq = blast_search.blast_search_for_gene(fasta_file, "neuA", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # toc_blast_run = time.perf_counter()
        # print("Took this long for isolate blast searching %s (s)" % (toc_blast_run - tic_iso_run))
        # tic_blast_process = time.perf_counter()

        # res_dict = { "id" : bassio_nameo,
        #            "flaA" : blast_processing.process_blast(flaA_seq, [182],"flaA"),
        #              "pilE" : blast_processing.process_blast(pilE_seq, [333],"pilE"),
        #             "asd" : blast_processing.process_blast(asd_seq, [473], "asd"),
        #              "mip" : blast_processing.process_blast(mip_seq, [402], "mip"),
        #              "mompS" : blast_processing.process_blast(mompS_seq, [352], "mompS"),
        #              "proA" : blast_processing.process_blast(proA_seq, [405], "proA"),
        #              "neuA" : blast_processing.process_blast(neuA_seq, [357,354,351], "neuA")}
        # res_df = pd.Series(res_dict).to_frame().transpose()
        # res_df.to_csv("./example_res_df.csv", index=None)

        # out_st = blast_processing.ST_process(data_dir, res_df)
        # out_df.loc[k] = out_st.iloc[0]
        # toc_blast_run = time.perf_counter()
        # print("Took this long for isolate blast processing %s (s)" % (toc_blast_run - tic_blast_process))

    toc_blast_run = time.perf_counter()
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print()
    print("Took this long for %s isolates %s (minutes)" % (len(seq_lines),round((toc_blast_run - tic_setup)/ 60, 3)))
    out_df.to_csv((input_args.output + "_ST.csv"), index=None)


def multiprocess_run(seq_lines):
    tot_list = mp_list(num_lines,input_args.threads)
    ## Set up the pool
    mp_pool = mp.Pool(processes=input_args.threads)
    rough_num = math.ceil(num_lines / input_args.threads)
    print_file = open(input_args.print_file, "a")
    print_file.write("Starting apply async runs " + str(datetime.datetime.now()) + "\n")
    print_file.write("Starting mem usage (MB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2) + "\n")
    print_file.close()
    df_list = [mp_pool.apply_async(table_getter, args=(input_args.embl[0], rough_num, rn, fs, fe,th)) for fs,fe,rn,th in tot_list]
    print_file = open(input_args.print_file, "a")
    print_file.write("Finishing apply async runs " + str(datetime.datetime.now()) + "\n")
    print_file.write("Starting mem usage (MB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2) + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()
    df_res = [p.get() for p in df_list]
    df_res = [p.dropna(how = "all") for p in df_res]
    print("")
    print("")
    newer_table = df_combiner(df_res)

def table_row_getter(fasta_file, data_dir, aa_dir_name, threads = 1):
    tic_iso_run = time.perf_counter()
    bassio_nameo = os.path.basename(fasta_file)
    print(bassio_nameo)

    ## Perform the ORF finder for the isolate

    #ORF_loc = orfipy_search.orfipy_search(fasta_file, input_args.threads)

    ## Search for cpn60
    flaA_seq = blast_search.blast_search_for_gene(fasta_file, "flaA", data_dir, aa_dir_name, bassio_nameo, threads)

    # Search for fusA
    pilE_seq = blast_search.blast_search_for_gene(fasta_file, "pilE", data_dir, aa_dir_name, bassio_nameo, threads)

    # Search for gltA
    asd_seq = blast_search.blast_search_for_gene(fasta_file, "asd", data_dir, aa_dir_name, bassio_nameo, threads)
    # Search for pyrG
    mip_seq = blast_search.blast_search_for_gene(fasta_file, "mip", data_dir, aa_dir_name, bassio_nameo, threads)
    # Search for recA
    mompS_seq = blast_search.blast_search_for_gene(fasta_file, "mompS", data_dir, aa_dir_name, bassio_nameo, threads)
    # Search for rplB
    proA_seq = blast_search.blast_search_for_gene(fasta_file, "proA", data_dir, aa_dir_name, bassio_nameo, threads)
    # Search for rpoB
    neuA_seq = blast_search.blast_search_for_gene(fasta_file, "neuA", data_dir, aa_dir_name, bassio_nameo, threads)
    toc_blast_run = time.perf_counter()
    print("Took this long for isolate blast searching %s (s)" % (toc_blast_run - tic_iso_run))
    tic_blast_process = time.perf_counter()

    res_dict = { "id" : bassio_nameo,
                "flaA" : blast_processing.process_blast(flaA_seq, [182],"flaA"),
                    "pilE" : blast_processing.process_blast(pilE_seq, [333],"pilE"),
                "asd" : blast_processing.process_blast(asd_seq, [473], "asd"),
                    "mip" : blast_processing.process_blast(mip_seq, [402], "mip"),
                    "mompS" : blast_processing.process_blast(mompS_seq, [352], "mompS"),
                    "proA" : blast_processing.process_blast(proA_seq, [405], "proA"),
                    "neuA" : blast_processing.process_blast(neuA_seq, [357,354,351], "neuA")}
    res_df = pd.Series(res_dict).to_frame().transpose()
    res_df.to_csv("./example_res_df.csv", index=None)

    out_st = blast_processing.ST_process(data_dir, res_df)
    #
    # 
    toc_blast_run = time.perf_counter()
    print("Took this long for isolate blast processing %s (s)" % (toc_blast_run - tic_blast_process))
    return out_st






