#!/bin/bash

# ------------------------------------------------------------------------------
##           Run 01_simulate_replicates_x_scenario across scenarios
# ------------------------------------------------------------------------------

## SLURM script for BSU icelake jobs:

##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Section 1: SLURM Commands 
##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#SBATCH --job-name=sim_array         # Job name
#SBATCH --array=1-3                # Job array
#SBATCH --output=logs/%x_%A_%a.out   # Standard output log (%x=job name, %A=job ID, %a=array task ID)
#SBATCH --error=logs/%x_%A_%a.err    # Error log
#SBATCH --time=15:00:00              # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=50           # Number of CPUs
#SBATCH -A mrc-bsu-sl2-cpu           # Project account name (mrc-bsu-sl2-cpu for icelake and mrc-bsu-sl2-gpu for ampere)
#SBATCH -p icelake                   # The partition. Use icelake for normal jobs, or icelake-himem if needed
#SBATCH --partition=standard         # Cluster partition or queue (depends on your system)

##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Section 2: Modules
##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
. /etc/profile.d/modules.sh               # Enables the module command
module purge                              # Removes all modules still loaded
module load rhel8/default-icl             # REQUIRED - loads the basic environment

# Load the latest R version
module load R/4.5.1-icelake

##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Section 3: Commands to execute
##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# All scenarios
PARAM_FILE="../../00_define_scenarios/outputs/scenario.00.csv"

# Row for each task ID
LINE=$(awk -v row=$SLURM_ARRAY_TASK_ID 'NR==row {print $0}' "$PARAM_FILE")
# Split line into args
read -r main_scenario main_scenario_val sub_sce_varying_par sub_sce_varying_par_value N r q1 q2 BGX1 diff_BGX BXY1 diff_BXY BUX BUY <<< "$LINE"
# Print scenario args
echo "Running task $SLURM_ARRAY_TASK_ID with parameters: main_scenario = $main_scenario, main_scenario_val = $main_scenario_val, sub_sce_varying_par = $sub_sce_varying_par sub_sce_varying_par_value = $sub_sce_varying_par_value" 

# Run the R script with index of row containing args x scenario
Rscript ../scripts/01_simulate_replicates_x_scenario.R "$SLURM_ARRAY_TASK_ID" 

##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Print relevant info:
JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then
        echo "Running on nodes: $SLURM_JOB_NODELIST"
else
        echo "Running on node: `hostname`"
fi

echo "Current directory: `pwd`"
echo -e "\nNum tasks = $SLURM_NTASKS, Num nodes = $SLURM_JOB_NUM_NODES, OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
