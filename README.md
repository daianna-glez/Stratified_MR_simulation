# Simulation-based estimation of power and type 1 error rate for stratified Mendelian Randomization
Repository for a stratified Mendelian Randomization (MR) simulation study, aimed to evaluate the performance of stratified MR studies according to variable stratum parameters.
This project is being developed under the supervision of [Dr. Janne Pott](https://pottj.github.io) and [Dr. Stephen Burgess](https://www.mrc-bsu.cam.ac.uk/staff/stephen-burgess) in the Medical Research Council (MRC) - Biostatistics Unit (BSU) at the University of Cambridge. 

## Background

- Disease risk and therapeutic treatment effects can differ across strata of the population. The characterization of such differences is crucial to identify subgroups of patients at higher risk for a disease or that would benefit the most from a therapeutic intervention.
- Stratified Genome Wide Association Studies (GWAS) and MR studies can reveal novel stratum-specific genetic and causal exposure associations and refine findings based on combined strata. 
- Through stratified drug target MR approaches (where variants mimic a drug effect on an exposure) on diseases and clinical traits, we can evaluate the feasibility of drug targets, explain comorbidities, explore drug repurposing opportunities, and investigate drug side effects in each stratum. 
- **Stratum-specific** results will refine our understanding of **disease aetiology** and potentially lead to **targeted, more effective treatments** and **prevention programs**.
- In stratified MR analyses, causal effects are estimated within each stratum either by stratifying individual level data and regressing the outcome on the genetically predicted exposure, or by computing IV Wald ratios using summary statistics from stratified GWAS on the exposure and the outcome (two samples). Stratum-specific causal effect estimates are then compared to determine if a difference exists between strata.

What’s not yet established is how between-strata differences in sample size, minor allele frequency (MAF) and strength of selected genetic instruments, as well as the amount of unmeasured confounding impact on the detection of true differences in causal effects (power) and the rejection of non-existent differences (type 1 error rate). When using GWAS summary statistics from two different samples (one for the exposure and a second for the outcome) additional between-study differences may further impact on downstream estimates of causal effect differences.

## Aims
To investigate the performance of stratified MR studies in the detection of differences in genetic and exposure causal effects across two population strata, based on simulated individual level data and GWAS summary statistics, and using one and multiple genetic instrumental variables (IVs). 


## Simulations

Used notation and formulas are introduced in the below table.

| Variable | Description |
|-----------|--------------|
| $N$ | Total sample size |
| $r = \frac{N_1}{N_2}$ | Ratio of stratum sample sizes |
| $q_k$ | MAF of variant in stratum $k$ |
| $\beta_{GXk}$ | True effect of genetic variant on the exposure in stratum $k$ |
| $\Delta{\beta_{GX}}= \beta_{GX2} - \beta_{GX1}$ | True difference in genetic effect across strata |
| $\beta_{XYk}$ | True causal effect in stratum $k$ |
| $\Delta{\beta_{XY}}= \beta_{XY2} - \beta_{XY1}$ | True difference in causal effect across strata |
| $\beta_{GYk}$ | True effect of genetic variant on the outcome in stratum $k$ |
| $\Delta{\beta_{GY}}= \beta_{GY2} - \beta_{GY1}$ | True difference in genetic effect on outcome across strata |
| $\beta_{UX}, \beta_{UY}$ | Effect of unmeasured confounder $U$ on exposure and outcome |


### 1. Individual level data and 1 IV simulations 
[`simulations/1.Individual_level_data_1IV/`](/simulations/1.Individual_level_data_1IV/)

Four **main scenarios** were defined according genetic and causal effect differences across the two strata:

- Scenarios 0,0: $\Delta{\beta_{GX}} = 0, \Delta{\beta_{XY}} = 0$
- Scenarios 0,1: $\Delta{\beta_{GX}} = 0, \Delta{\beta_{XY}} ≠ 0$
- Scenarios 1,0: $\Delta{\beta_{GX}} ≠ 0, \Delta{\beta_{XY}} = 0$
- Scenarios 1,1: $\Delta{\beta_{GX}} ≠ 0, \Delta{\beta_{XY}} ≠ 0$

For main scenarios with nonzero differences (i.e., 0,1, 1,0, and 1,1), **main subscenarios** were defined according to the size of the differences. For each main subscenario, **varying parameter subscenarios** were defined according to one or two parameters varying (16 total cases), and definite subscenarios were 
defined according to the values those varying parameters take. See details in [`1.Individual_level_data_1IV/`](/simulations/1.Individual_level_data_1IV/).

Roughly, code is divided into: 0) definition of all subscenarios to simulate (from all main scenarios and main subscenarios), 1) data generation and statistical association testing across replicates for each subscenario, and 2) aggregation and plotting of results across replicates of a subscenario, and across subscenarios of a main subscenario.
Supplementary analyses are described in [`1.Individual_level_data_1IV/`](/simulations/1.Individual_level_data_1IV/).

### 2. Individual level data and multiple IVs simulations 
`TODO`

### 3. GWAS summary statistics and one IVs simulations 
`TODO`

### 4. GWAS summary statistics and multiple IVs simulations 
`TODO`


## File organization
Each simulation study has its own folder. Each analysis performed has its own subfolder: all code its under the `scripts/` directory, all used inputs in `inputs/`, generated data outputs in `outputs/`, and plots derived from such analysis in `plots/`. 
Simulation results are further categorized by main scenario (0,0, 0,1, 1,0, and 1,1), main subscenario (given by value(s) of difference(s) in 0,1, 1,0, and 1,1 cases), varying parameter(s) subscenario, and value(s) of varying parameter(s) subscenario. 


## Internal access
- BSU Server repository path: `/localhome/daianna/Stratified_MR_simulation/`
- HPC repository path: `TODO`

