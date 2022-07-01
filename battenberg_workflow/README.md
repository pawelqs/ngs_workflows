# Snakemake workflow for running battengerg CNV calling

The snakemake is created to easily be run on computing clusters with SLURM. 

**File under construction**

## Requirements

* apptainer
* python environment with snakemake and pandas

## Usage:

### HPC + SLURM

Open screen session and run

```
snakemake --profile config/hpc_slurm_profile --configfile config/config.yml
```

## Limitations

It is not a part of clonality workflow, since it requires apptainer, which can only be run on wn136 (phi)

As for now does not work - ends with an error:

> snakemake --profile config/hpc_slurm_profile --configfile config/config.yml --use-singularity
Building DAG of jobs...
WorkflowError:
The singularity command has to be available in order to use singularity integration.
