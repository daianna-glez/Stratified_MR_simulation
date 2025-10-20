
# 1. Simulation study based on individual level data for 1 IV

## Aims

To evaluate the effect of between-strata varying parameters on power and type I error rate for the detection of causal effect differences across two population strata ($\Delta{\beta_{XY}}$). 

As sanity checks, to evaluate the effect of varying these parameters on the estimation of:

1) genetic effects on the exposure per strata ($\beta_{GXk}$), 
2) genetic effects on the outcome per strata ($\beta_{GYk}$), 
3) causal effects per strata ($\beta_{XYk}$), 
4) the difference in genetic effect on the exposure between strata ($\Delta{\beta_{GX}}$), and 
5) the difference in genetic effect on the outcome between strata ($\Delta{\beta_{GY}}$). 

The effect of these parameters was interrogated for each in isolation and in combination with others, exploring scenarios that can be expected in practice when performing stratified GWAS and MR studies across different ancestry, sex, or age groups. 

## Data-generating mechanisms

### 0. Definition of simulation scenarios
[`00_define_scenarios/`](/00_define_scenarios/) folder. 

---

All subscenarios to simulate were derived from four main scenarios defined according to zero and non-zero $\Delta\beta_{GX}$ and $\Delta\beta_{XY}$:

- **Scenarios 0,0**: no genetic effect difference $\Delta{\beta_{GX}} = 0$ and no causal effect difference $\Delta{\beta_{XY}} = 0$. 
- **Scenarios 0,1**: no genetic effect difference $\Delta{\beta_{GX}} = 0$ and causal effect difference $\Delta{\beta_{XY}} ≠ 0$. 
- **Scenarios 1,0**: genetic effect difference $\Delta{\beta_{GX}} ≠ 0$ and no causal effect difference $\Delta{\beta_{XY}} = 0$. 
- **Scenarios 1,1**: genetic effect difference $\Delta{\beta_{GX}} ≠ 0$ and causal effect difference $\Delta{\beta_{XY}} ≠ 0$. 

Each subscenario is described by $N, r, q_1, q_2, \beta_{GX1}, \Delta{\beta_{GX}}, \beta_{XY1}$, $\Delta\beta_{XY}$, $\beta_{UX}$, and $\beta_{UY}$. 

The **base case** used was: $N$ = 50,000, $r$ = $1$, $q_1$ = $q_2$ = 0.25,  $\beta_{GX1}$ = 0.5, $\beta_{XY1}$  = 0.7, $\beta_{UX}$ = $\beta_{UY}$ = 1.

The following were the values considered for each parameter, spanning ranges based on large-scale study observations or restricted by their own definition (MAF):

1. $N$ = {10k, 20k, 30k, 40k, 50k, 60k, 70k, 80k, 90k, 100k}
2. $r$ = {0.20,  0.25,  0.33,  0.50,  1.0,  2.0,  3.0,  4.0,  5.0} 
3. $q_1, q_2$ = {0.05,  0.15,  0.25,  0.35,  0.45}
4. $\beta_{GX1}$ = {0.09,  0.25,  0.41,  0.57,  0.73,  0.89,  1.05,  1.21,  1.37}
5. $\beta_{XY1}$ = {0.00,  0.15,  0.30,  0.45,  0.60,  0.75,  0.90,  1.05,  1.20,  1.35}
6. $\beta_{UX}$ = $\beta_{UY}$ = {0.00,  0.25,  0.50,  0.75,  1.00, 1.25,  1.50}
7. $N$ = {10k, 40k, 70k, 100k} and $r$ = {0.20,  0.50,  1.0,  2.0, 5.0} 
8. $N$ = {10k, 40k, 70k, 100k} and $q_1, q_2$ = {0.05,  0.25,  0.45}
9. $N$ = {10k, 40k, 70k, 100k} and $\beta_{GX1}$ = {0.1,  0.4,  0.7,  1.0,  1.3} 
10. $N$ = {10k, 40k, 70k, 100k} and $\beta_{XY1}$ = {0.0,  0.3,  0.6,  0.9,  1.2}
11. $r$ = {0.20,  0.50,  1.0,  2.0, 5.0} and $q_1, q_2$ = {0.05,  0.25,  0.45}
12. $r$ = {0.20,  0.50,  1.0,  2.0, 5.0} and $\beta_{GX1}$ = {0.1,  0.4,  0.7,  1.0,  1.3} 
13. $r$ = {0.20,  0.50,  1.0,  2.0, 5.0} and $\beta_{XY1}$ = {0.0,  0.3,  0.6,  0.9,  1.2}
14. $q_1, q_2$ = {0.05,  0.25,  0.45} and $\beta_{GX1}$ = {0.1,  0.4,  0.7,  1.0,  1.3} 
15. $q_1, q_2$ = {0.05,  0.25,  0.45} and $\beta_{XY1}$ = {0.0,  0.3,  0.6,  0.9,  1.2}
16. $\beta_{GX1}$ = {0.1,  0.4,  0.7,  1.0,  1.3} and $\beta_{XY1}$ = {0.0,  0.3,  0.6,  0.9,  1.2}

**Scenarios 0,0**: $\Delta{\beta_{GX}}= 0, \Delta{\beta_{XY}}= 0$

Defined in [`scripts/01_define_00_scenarios.R`](/00_define_scenarios/scripts/01_define_00_scenarios.R).

Assess the rate of falsely detected between-strata differences (i.e. type I error rate) in the genetic effect on the exposure and the exposure effect on the outcome, across the 16 varying parameter(s) subscenarios presented above.

**Scenario 0,1**: $\Delta{\beta_{GX}}= 0, \Delta{\beta_{XY}}≠ 0$

Defined in [`scripts/02_define_01_scenarios.R`](/00_define_scenarios/scripts/02_define_01_scenarios.R).

Assess the rate of falsely detected between-strata differences in the genetic effect on the exposure and the rate of truly detected between-strata differences (i.e. power) in the causal effect of the exposure on the outcome, for the following difference sizes:

1. $\Delta\beta_{XY}$ = -1.7
2. $\Delta\beta_{XY}$ = -1.2
3. $\Delta\beta_{XY}$ = -0.7
4. $\Delta\beta_{XY}$ = -0.2
5. $\Delta\beta_{XY}$ = 0.3
6. $\Delta\beta_{XY}$ = 0.8
7. $\Delta\beta_{XY}$ = 1.3, 

each across the 16 subscenarios presented above.

**Scenario 1,0**: $\Delta{\beta_{GX}}≠ 0, \Delta{\beta_{XY}}= 0$

Defined in [`scripts/03_define_10_scenarios.R`](/00_define_scenarios/scripts/03_define_10_scenarios.R).

Assess the rate of truly detected between-strata differences in the genetic effect on the exposure and the rate of falsely detected between-strata differences in the causal effect of the exposure on the outcome, for the following difference sizes:

1. $\Delta\beta_{GX}$ = -0.6
2. $\Delta\beta_{GX}$ = -0.2
3. $\Delta\beta_{GX}$ = 0.2
4. $\Delta\beta_{GX}$ = 0.6
5. $\Delta\beta_{GX}$ = 1.0
6. $\Delta\beta_{GX}$ = 1.4

each across the 16 subscenarios presented above.

**Scenario 1,1**: $\Delta{\beta_{GX}}≠ 0, \Delta{\beta_{XY}}≠ 0$

Defined in [`scripts/04_define_11_scenarios.R`](/00_define_scenarios/scripts/04_define_11_scenarios.R).

Assess the rate of truly detected between-strata differences in the genetic effect on the exposure and the causal effect of the exposure on the outcome, for each combination of $\Delta\beta_{GX}$  = {-0.6  -0.2  0.2  0.6  1.0  1.4} and $\Delta\beta_{XY}$  = {-1.7  -1.2  -0.7  -0.2  0.3  0.8  1.3}, each, across the 16 subscenarios presented above.

In total, 896 subscenarios were explored.



### 1. Run simulation replicates per subscenario 
[`01_run_replicates_x_scenario/`](/01_run_replicates_x_scenario/) folder. 

---

**Simulation framework** 

According to each subscenario parameters, we simulated additively modelled genotypes $G$ and a continuous exposure $X$ and outcome $Y$ for unrelated individuals across the strata of a dichotomous stratifying variable K $=k\in$ {1,2}, genetic associations with $X$ and $Y$ were tested and causal effects estimated in [`scripts/01_simulate_replicates_x_scenario.R`](/01_run_replicates_x_scenario/scripts/01_simulate_replicates_x_scenario.R) for 100 replicates. 
Per-subscenario replicates were run with [`scripts/02_run_01_across_scenarios.R`](/01_run_replicates_x_scenario/scripts/02_run_01_across_scenarios.R) taking each subscenario parameters from [`00_define_scenarios/outputs`](/00_define_scenarios/outputs/). This was run separately for all subscenarios of each main scenario (0,0, 0,1, 1,0, and 1,1).  

For each replicate $i = 1,...,100$ of a subscenario we:

1. Simulate $N$ individuals across the two strata: $K | N, r$
2. Simulate genotype $G$ for the individuals in each stratum $k$: $G_k |N_k, q_k$
    
    $G_k \sim Binom(2, q_k)$ 
    
    2.1 Perform chi-square test for HWE in each stratum and globally.  
    
3. Simulate unknown confounder $U$ for all $N$ individuals
    
    $U \sim Unif(0,1)$ 
    
4. Simulate error terms for X and Y for all $N$ individuals 
    
    $\epsilon_{X} \sim N(0,1)$ and $\epsilon_{Y} \sim N(0,1)$ 
    
5. Generate exposure per stratum $X|K, G, U, \epsilon_X$, $\beta_{GX1}$, $\Delta\beta_{GX}$,  $\beta_{UX}$

    - For stratum 1: $X ={\beta_{GX1}}G+\beta_{UX}U+\epsilon_{X}$
    - For stratum 2: $X=({\beta_{GX1}}+\Delta\beta_{GX})G+\beta_{UX}U+\epsilon_{X}$
    
6. Test association between G and X in each stratum
    
    6.0 Regress $X \sim G$ separately in each stratum: get $\hat\beta_{GXk}$, se($\hat\beta_{GXk}$), $t_{\hat\beta_{GXk}}$, and $p_{\hat\beta_{GXk}}$
    
    6.1 Confirm variant is strong IV in both strata: F-statistic > 10 rule for IV in each stratum.
    
    *Sanity checks:
      - Difference between true and estimated genetic effect on X per stratum: $|{\beta_{GXk}} - {\hat\beta_{GXk}}|$
      - Difference between true and estimated genetic effect difference across strata: $|{\Delta{\beta_{GX}}} - {\Delta\hat\beta_{GX}}|$
      - Test significance of estimated difference:  $Z_{\Delta\hat\beta_{GX}}=\frac{{\hat\beta_{GX2}} - {\hat\beta_{GX1}}}{\sqrt{se^2(\hat\beta_{GX2}) + se^2(\hat\beta_{GX1})}}$ with $Z \sim N(0,1)|H_0$

7. Generate outcome $Y|K, G, U, \epsilon_Y, {\beta_{GX1}}, \Delta\beta_{GX}, \beta_{XY1}, \Delta\beta_{XY}$, $\beta_{UX}$, $\beta_{UY}$

      - For stratum 1: $Y ={\beta_{XY1}}({\beta_{GX1}}G + \beta_{UX}U)+\beta_{UY}U+\epsilon_{Y}$
      - For stratum 2:  $Y =({\beta_{XY1}}+\Delta\beta_{XY})(({\beta_{GX1}}+\Delta\beta_{GX})G + \beta_{UX}U)+\beta_{UY}U+\epsilon_{Y}$
      
      
8.  Test association between G and Y in each stratum 
    
    8.0 Regress $Y \sim G$ separately in each stratum: get $\hat\beta_{GYk}$, se($\hat\beta_{GYk}$), $t_{\hat\beta_{GYk}}$, and $p_{\hat\beta_{GYk}}$
    
    *Sanity checks:
    
    - Difference between true and estimated genetic effect on Y per stratum: $|{\beta_{GYk}} - {\hat\beta_{GYk}}|$
    - Difference between true and estimated genetic effect difference across strata: $|{\Delta{\beta_{GY}}} - {\Delta\hat\beta_{GY}}|$
    - Test significance of estimated difference: $Z_{\Delta\hat\beta_{GY}}=\frac{{\hat\beta_{GY2}} - {\hat\beta_{GY1}}}{\sqrt{se^2(\hat\beta_{GY2}) + se^2(\hat\beta_{GY1})}}$ with $Z \sim N(0,1)|H_0$


9. Estimate causal effects per stratum: get $\hat\beta_{XYk}$ 

    - For stratum 1: $\hat\beta_{XY1} = \frac{\hat\beta_{GY1}}{\hat\beta_{GX1}}$
    - For stratum 2: $\hat\beta_{XY2} = \frac{\hat\beta_{GY2}}{\hat\beta_{GX2}}$
    
    9.1 Calculate metrics of interest:

      - Difference between true and estimated causal effect per stratum: $|{\beta_{XYk}} - {\hat\beta_{XYk}}|$
      - Difference between true vs estimated causal effect difference across strata:  $|{\Delta{\beta_{XY}}} - {\Delta\hat\beta_{XY}}|$
      - Standard error of estimated causal effect as the first order term from a delta method expansion: se($\hat\beta_{XYk}$) =$\frac{se(\hat \beta_{GY_k})}{\hat\beta_{GX_k}}$
      - Standard error of estimated causal effect using the second term of the delta expansion: se($\hat\beta_{XYk}$) = $\sqrt{}\frac{se(\hat\beta_{GY_k})^2}{\hat\beta_{GX_k}^2}+\frac{\hat\beta_{GY_k}^2se(\hat\beta_{GX_k})^2}{\hat\beta_{GX_k}^4}$
      - Test significance of estimated difference: $Z_{\Delta\hat\beta_{XY}}=\frac{{\hat\beta_{XY2}} - {\hat\beta_{XY1}}}{\sqrt{se^2(\hat\beta_{XY2}) + se^2(\hat\beta_{XY1})}}$ with $Z \sim N(0,1)|H_0$


## Estimands and performance assessment

### 2. Summarize and plot results per subscenario 
[`02_summarize_and_plot_results/`](/02_summarize_and_plot_results/) folder.

---

The primary estimand of interest is $\Delta\beta_{XY}$. Secondary estimands (as part of sanity checks) are $\beta_{GX_k}$, $\beta_{GY_k}$, $\beta_{XY_k}$, $\Delta\beta_{GX}$, and $\Delta\beta_{GY}$. 

Each simulation run for a subscenario replicate yielded the following outputs of interest:

| Output | Description |
|-----------|--------------|
| $q_{k_{ob}}$ | Observed MAF of variant in stratum $k$ (or globally) |
| $\chi^2_k$ | Chi-square statistic for HWE test in stratum $k$ (or globally) |
| $P_{\chi^2_k}$ | *p*-value of chi-square statistic for HWE test in stratum $k$ (or globally) |
| $\hat\beta_{GXk}$ | Estimated genetic effect on the exposure in stratum $k$ |
| se($\hat\beta_{GXk}$) | Standard error for $\hat\beta_{GXk}$ |
| $t_{\hat\beta_{GXk}}$ | *t*-statistic for $\hat\beta_{GXk}$ |
| $p_{\hat\beta_{GXk}}$ | *p*-value of $t_{\hat\beta_{GXk}}$ |
| $F_{\hat\beta_{GXk}}$ | F-statistic for $\hat\beta_{GXk}$ |
| $\hat\Delta{{\beta_{GX}}}= {\hat\beta_{GX2}} - {\hat\beta_{GX1}}$ | Estimated difference in genetic effect on exposure across strata |
| $Z_{\Delta\hat\beta_{GX}}$ | Z-statistic for $\hat\Delta{{\beta_{GX}}}$ | 
| $p_{Z_{\Delta\hat\beta_{GX}}}$ | *p*-value for $Z_{\Delta\hat\beta_{GX}}$ |
| $\hat\beta_{GYk}$ | Estimated genetic effect on the outcome in stratum $k$ |
| se($\hat\beta_{GYk}$) | Standard error for $\hat\beta_{GYk}$ |
| $t_{\hat\beta_{GYk}}$ | *t*-statistic for $\hat\beta_{GYk}$ |
| $p_{\hat\beta_{GYk}}$ | *p*-value of $t_{\hat\beta_{GYk}}$ |
| $\hat\Delta{{\beta_{GY}}}= {\hat\beta_{GY2}} - {\hat\beta_{GY1}}$ | Estimated difference in genetic effect on outcome across strata |
| $Z_{\Delta\hat\beta_{GY}}$ | Z-statistic for $\hat\Delta{{\beta_{GY}}}$ | 
| $p_{Z_{\Delta\hat\beta_{GY}}}$ | *p*-value for $Z_{\Delta\hat\beta_{GY}}$ | 
| $\hat\beta_{XYk}$ | Estimated causal effect in stratum $k$ |
| se($\hat\beta_{XYk}$) | Standard error for $\hat\beta_{XYk}$ |
| se(${\hat\beta_{XYk}}$)$^{2nd}$ | Standard error for $\hat\beta_{XYk}$ including second term from delta expansion method |
| $\hat\Delta{{\beta_{XY}}}= {\hat\beta_{XY2}} - {\hat\beta_{XY1}}$ | Estimated difference in causal effect across strata |
| $Z_{\Delta\hat\beta_{XY}}$ | Z-statistic for $\hat\Delta{{\beta_{XY}}}$ | 
| $p_{Z_{\Delta\hat\beta_{XY}}}$ | *p*-value for $Z_{\Delta\hat\beta_{XY}}$ |

The above per-replicate metrics were aggregated and summarized across all 100 replicates of a subscenario in [`scripts/01_summarize_scenario_results.R`](/02_summarize_and_plot_results/scripts/01_summarize_scenario_results.R). Performance metrics considered across replicates were:

- Bias (mean error)
- Mean absolute error (MAE)
- Type I error rate
- Power

Results across subscenarios of a varying parameter scenario were plotted in [`scripts/02_plot_scenario_results.R`](/02_summarize_and_plot_results/scripts/02_plot_scenario_results.R) and [`scripts/03_plot_main_scenario_results.R`](/02_summarize_and_plot_results/scripts/03_plot_main_scenario_results.R).


---

## Supplementary analyses
[`Supplementary_analyses/`](/Supplementary_analyses/) folder.

### 1.0  Simulations with varying βxʏ₁ and βux effect on Y.
[`01_BXY_U_effects/`](/Supplementary_analyses/01_BXY_U_effects/) subfolder.

Repeated simulations for 0,0 increasing $\beta_{XY1}$ subscenarios with zero and non-zero $\beta_{UX}$ on $Y$ to examine the role of the unmeasured confounder $U$'s indirect effect on $Y$ (via $X$) on genetic and causal effect estimates.


