# Optimal threshold
This R code makes an easier access to the supplementary material of the papier entitled **on evaluating how well a biomarker can predict treatment response with survival data**.

The function 
- **check_effect** for checking the treatment effect.

- **simul_data** is a function to simulated data using three scenarios: constant, increasing and decreasing  effect of treatment. 

- **get_risk** is a function to get risk at a prediction time (timepoint) with treated(risk_1) and untreated(risk_1) subjects..

- **get_censoring_weights** returns the inverse of the probability of censoring weighting. This function help to take the standard decision at a prediction time t. 

- **plotime_predictiveness_curve** displays time-dependent marker-by-treatment predictiveness curves.

R code 

- **extention_dataset_code** is the script using in section 5.3 in the paper to extend the sample size of the real datasets. 
