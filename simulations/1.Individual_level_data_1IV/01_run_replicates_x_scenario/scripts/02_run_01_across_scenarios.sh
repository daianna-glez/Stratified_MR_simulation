#!/bin/bash

# ------------------------------------------------------------------------------
#             Run 01_simulate_replicates_x_scenario across scenarios
# ------------------------------------------------------------------------------

# All scenarios (uncomment according to main scenario run)
# PARAM_FILE="../../00_define_scenarios/outputs/scenario.00.csv"
PARAM_FILE="../../00_define_scenarios/outputs/scenario.01.csv"

n=$(cat  $PARAM_FILE | wc -l)

for i in $(seq 1 $n)
do

# Row for scenario
LINE=$(awk -v row=$i 'NR==row {print $0}' "$PARAM_FILE")
# Split line into args
read -r main_scenario main_scenario_val sub_sce_varying_par sub_sce_varying_par_value N r q1 q2 BGX1 diff_BGX BXY1 diff_BXY BUX BUY <<< "$LINE"
# Print scenario args
echo "Running task $i with parameters: main_scenario = $main_scenario, main_scenario_val = $main_scenario_val, sub_sce_varying_par = $sub_sce_varying_par sub_sce_varying_par_value = $sub_sce_varying_par_value" 

# Run the R script with index of row containing args x scenario
Rscript ../scripts/01_simulate_replicates_x_scenario.R "$i" 

done

