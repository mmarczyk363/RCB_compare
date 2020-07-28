# RCB_compare
R code to compare distribution of RCB scores from patients treated with one type of therapy (e.g. standard of care) versus another type (e.g. experimental therapy) to quantify the benefit of the experimental treatment defined as treatment cohort-wide down shift in RCB. 

We propose a new metric, Treatment Efficacy Score (TES) that can capture the relative efficacy of one treatment over the other. We explored three different methods to calculate TES:

•	Weighted two-sample Kolmogorov-Smirnov test (wKS)

•	Density ratio of RCB scores from two treatments (DensRatio)

•	Density difference of RCB scores from two treatments (DensDiff)
