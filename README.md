# Optimal threshold
This R code makes an easier access to the supplementary material of the papier entitled **on evaluating how well a biomarker can predict treatment response with survival data**.

The function 

- **simul_data** is a function to simulate data using three scenarios: constant, increasing and decreasing  effect of treatment. 

- **get_estimates** is a function to get at horizon time (timepoint)

   * the optimal threshold (threshold)
   * the confidence interval of the threshold(confint)
   * the proportion of subject who avoid or benefict the treatment (Pneg*100)
   * the risk in treated(risk_1) subjects 
   * the risk in untreated(risk_0) subjects
   

- **Marginal_effect_timepoint** returns rho_0(t)-rho_1(t). 

- **plotime_predictiveness_curve** displays time-dependent marker-by-treatment predictiveness curves.

R code 

- **extention_dataset_code** is the script using in section 5.3 in the paper to extend the sample size of the real datasets. 


