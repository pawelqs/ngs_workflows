# Snakemake workflow for clonality analysis of WES data

The snakemake is created to easily be run on computing clusters with SLURM.

**File under construction**

## Requirements

* conda
* python environment with snakemake and pandas
* ...

## Installation of snakemake & pandas

```
# Install Mamba forge, or install mamba in an existing conda base environment
conda install mamba -n base -c conda-forge

# Install snakemake
mamba create -c conda-forge -c bioconda -n snakemake snakemake pandas

# if you experience error: module 'lib' has no attribute 'OpenSSL_add_all_algorithms'
mamba create -c conda-forge -c bioconda -c anaconda -n snakemake snakemake cryptography==38.0.4 pandas
```

## Usage:

### HPC + SLURM

Open screen session and run

```
conda activate snakemake
snakemake --profile config/hpc_slurm_profile --configfile config/config.yml
# or
snakemake --profile config/hpc_slurm_profile --configfiles config/config.hg19.yml config/Shlush_AML.yml
```

## Limitations
