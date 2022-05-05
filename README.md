# lp_mlst
Seq types of Legionella pneumophila assemblies

## Installation ##

Install conda first install conda first [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install). <br />
Once install clone this repo:

	git clone https://github.com/jdaeth274/lp_mlst

Then use the environment.yml file to install the dependencies with conda

	cd ./lp_mlst
	conda env create --file=environment.yml

Then activate the environment with:

	conda  activate lp_mlst_env

## Usage ##
To find the Seq types of a collection of Legionella pneumophila assemblies, first concatenate the full paths to your <br />
assemblies into a single file:

	ls -d $PWD/*.fasta > list_of_fastas.txt

Then run the following command to output a `.csv` file containing the file name, ST and full profile of the isolate

	python lp_mlst_runner.py --seqs <list of fastas> --output <outfile_prefix> --threads <int number of threads to use>

Uses MLST data from [mompS](https://github.com/bioinfo-core-BGU/mompS)


