#!/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=survivability
#SBATCH --account=coen
#SBATCH --output=survivability-%j-%N.out
#SBATCH --error=survivability-%j-%N.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=END
#SBATCH --mail-user=strenge@control.tu-berlin.de

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.4.2
module load hpc/2015
julia  survivability_setup_error.jl $SLURM_NTASKS
