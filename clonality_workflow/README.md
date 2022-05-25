# Snakemake workflow for clonality analysis of WES data

The snakemake is created to easily be run on computing clusters with SLURM. For `qsub` `run_cnv_workflow.sh` will need to be adjusted.

**File under construction**

## Requirements

* conda
* python environment with snakemake and pandas
* ...

## Usage:

### HPC + SLURM

Open screen session and run

```
snakemake --profile config/hpc_slurm_profile --configfile config/config.yml
```

## Limitations
