.ONESHELL:


create-rseqc-env:
	# You need to replace the following line with the rogiht path to your conda installation
	source /home/pkus/programs/miniconda/build/etc/profile.d/conda.sh
	conda create -y --name rseqc_env python=3.11
	conda activate rseqc_env
	pip install RSeQC


run-rnaseq-workflow:
	snakemake --profile config/hpc_slurm_profile


unlock:
	snakemake -j 1 --unlock
