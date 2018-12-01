# Optimal threshold
R-code to display  the time-dependent marker by treatment predictiveness curve  using time-depenedent logistic model.

The function **check_effect** :

The function **simul_data** is a function to simulated data using three scenarios: constant, increasing and decreasing  effect of treatment. 

The function **get_risk** is a function to get risk at a prediction time (timepoint) with treated and untreated subjects..

The function **get_censoring_weights** returns the inverse of the probability of censoring weighting. This function help to take the standard decision at a prediction time t. 

The function **plotime_predictiveness_curve**  to display time-dependent marker-by-treatment predictiveness curves.

The **extention_dataset_code** is the script using the section 5.3 in the paper to extend the sample size of the real datasets. 



