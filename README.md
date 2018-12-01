# Optimal threshold
R-code to choose the optimal threshold of a biomarker using time-depenedent logistic model.

The function **check_effect** in R takes 4 inputs:
- time is the time of event
-event is the status of event for example died or not and takes 0 or 1.
- treatment is the indicator of the treatment and takes 0 or 1.
- Marker is the contnuous biomarker.
and outputs the plots of the effects of the covariables (intercepts,treatment,Marker and interaction between treatment and biomarker. 

The function **simul_data** is a function to simulated data using three scenarios: constant, increasing and decreasing  effect of treatment. This function takes 2 inputs.
- n: the sample size 
- scenario : scenario=1 if constant effect, scenario=2 if increasing effect and scenario=3 if decreasing effect.
and outputs a data set with 4 variables: time,event, treatment, biomarker.

The function **get_risk** is a function to get risk at a prediction time (timepoint) with treated and untreated subjects. The inputs are the same to **check_effect** also 




Inputs: values - This is a single column matrix and contains values between 0 and 1. classes - This single colomn matrix contains the labels for each of the values mentioned under "values". The labels should just have 2 classes i.e. 0 and 1
