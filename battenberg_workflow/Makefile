.ONESHELL:

SHELL = /bin/bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

dry-run:
	$(CONDA_ACTIVATE) snakemake
	snakemake -pn --profile config/hpc_slurm_profile --configfile config/config.yml

run:
	$(CONDA_ACTIVATE) snakemake
	snakemake --profile config/hpc_slurm_profile --configfile config/config.yml --use-singularity
