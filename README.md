# CovSocial PGS Analysis Code 
Analysis code used to produce main and supplementary results of the 
manuscript: "Genetic Predisposition for Negative Affect Predicts Mental 
Health Burden during the Covid-19 Pandemic"

## 00_helper_functions.R
Formatting function to add asterix to multinomial results, Table 3. 

## 01_sample_characteristics.R
Code to produce Table 1 (sample characteristics), supplementary file 3, 
Table S1 (comparison to overal CovSocial sample), and the 
resilience-vulnerability plot of the analytic sample (Fig 1) 

## 02_PGS_cors_validity_specificity.R 
Computes Pearson's correlation between all PGS, and test for association 
between each individual PGS, and the target phenotype (e.g., neuroticism 
PGS - trait neuroticism). Output of this script are shown in supplementary 
file 3, Table S2 and Figures S1 to S6. 

## 03_PGS_aov_latent.R
Tests for PGS mean differences between the trajectories. Output of this 
one-way ANOVA is shown in Table 3.  

## 04_PGS_multinom_latent.R
For those PGS with mean differences between the trajectories, we applied 
multinomial regression to test which trajectories can be distinguished 
based on the PGS. Results are shown in Table 4 and Fig 2.  

## 05_PGS_lm_latent.R 
In this script, we run linear regression to test associations between 
isolated latent baseline and change scores that form the basis of the 
resilience-vulnerability trajectories. Output is shown in Table 5.
