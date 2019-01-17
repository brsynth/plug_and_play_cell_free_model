# The aim of this file is to generate the model heatmap data as presented in #add ref
# It focuses on the TF/reporter data generation.

#Author: Mathilde Koch, JL Faulon's research group, INRA

# Cleaning previous R workspace
rm(list = ls())

# Importing necessary packages: 

library(deSolve)  # For solving ODEs
library(ggplot2)  # For visualisation
library(reshape2)   # For visualisation (array reshaping)


solve_for_x <- function(X, epsilon, 
                        n_GFP_here = n_GFP, n_enz_here = n_enz, n_TF_here = n_TF,
                        kappa_GFP_here = kappa_GFP, kappa_enz_here = kappa_enz, kappa_TF_here = kappa_TF
                        ) {
  x_search  = c(0:X)
  # Total sum of free RNAP + genes has to be closest to allowed number of RNAP
  x_error = (x_search + epsilon *n_GFP_here * x_search/(x_search + kappa_GFP_here) + 
               n_enz_here * x_search/(x_search + kappa_enz_here) + 
               n_TF_here * x_search/(x_search + kappa_TF_here) - X)/X
  
  x = which(abs(x_error) == min(abs(x_error)))
  return(x)
}

solve_for_y <- function(Y,
                        mRNA_GFP_here = mRNA_GFP,  mRNA_enz_here =mRNA_enz,
                        mRNA_TF_here =mRNA_TF,
                        k_GFP_here = k_GFP, k_enz_here = k_enz, k_TF_here = k_TF 
) {
  y_search  = c(0:Y)
  # total sum of free ribosomes + transcripts has to be closest to allowed number of ribosomes
  y_error = (y_search + mRNA_GFP_here * y_search/(y_search + k_GFP_here) +
               mRNA_enz_here * y_search/(y_search + k_enz_here) +
               mRNA_TF_here * y_search/(y_search + k_TF_here) - Y)/Y
  
  y = which(abs(y_error) == min(abs(y_error)))
  return(y)
}

# Model parameters

global_parameters = c(
        kappa_GFP = 100,
        kappa_TF = 3000,
        k_GFP = 1e3,
        k_TF = 1e4,
        gamma_GFP = 0.42,
        gamma_TF = 0.335,
        delta_GFP = 0.05,
        delta_enz = 0.05,  # for both CocE and HipO
        delta_TF = 0.05,
        pi_GFP = 3.75,
        pi_TF = 3.679,
        lambda_GFP = 0.0016,
        lambda_enz = 0.0016,  # for both CocE and HipO
        lambda_TF = 0.0016,
        n_copy_GFP = 100,
        n_copy_TF = 30,
        omega_GFP = 4, 
        omega_TF = 5, 
        benzoic_kd = 100,
        activated_kd = 50,
        X_tot = 30,
        Y_tot = 30,
        K_R_mRNA = 10,
        K_tox = 100
)

pars_biosensor = c(
  kappa_enz = 10,  # To avoid division by 0 in resource function, not used
  k_enz = 10,  # To avoid division by 0 in resource function, not used
  gamma_enz = 0,
  pi_enz = 0,
  omega_enz = 0,
  Km_enz = 0,
  k_cat_enz = 0,
  mRNA_enz_length = 0,
  spontaneous_hydrolisation = 0,
  n_copy_enz = 0
)

# Defining the parameters as a mix of global parameters and parameters specific to the biosensor
# This allows for easy comparison with the enzyme assay script

all_parameters = c(global_parameters, pars_biosensor)

name = "Benzoic acid"
colnames_for_visu = c("Transcription factor", 'Reporter', 'GFP')

init_solver <- c(
  mRNA_GFP = 0,
  mRNA_TF = 0,
  mRNA_enz = 0,
  GFP = 0, 
  TF = 0,
  enz = 0,
  benzoic_acid = 0,
  inducer = 0,
  R_mRNA = 10000000,
  tox = 0
)

solver_all_resource <- function(t, state, pars) {
  with (as.list(c(state, pars)),{

    if (k_cat_enz == 0) {
      dbenzoic_acid = 0
      dinducer = 0
      benzoic_acid = inducer
    } else {
      rate_conversion = 0.001 * enz * k_cat_enz * inducer/(inducer + Km_enz) # Rate is expressed In µM. * 0.001 since enzyme is in nM.
      dbenzoic_acid = rate_conversion + spontaneous_hydrolisation * inducer
      dinducer = - rate_conversion - spontaneous_hydrolisation * inducer
    }

    activated_TF = TF * (benzoic_acid/(benzoic_acid + benzoic_kd) + 0.0005)  # Hill activation of the TF by benzoic acid
    epsilon = activated_TF/(activated_TF + activated_kd * 1000)  # Hill activation of the responsive promoter by activated TF.
    # See main text, modifying copy numbers to account for RNAP per DNA strand (see main text)
    n_GFP = n_copy_GFP * omega_GFP
    n_enz = n_copy_enz * omega_enz
    n_TF = n_copy_TF * omega_TF
    
    # Solving for RNAP at each time step
    x = solve_for_x(X = X_tot, epsilon = epsilon, n_GFP_here = n_GFP, 
                    n_enz_here = n_enz, n_TF_here = n_TF,
                    kappa_GFP_here = kappa_GFP,
                    kappa_enz_here = kappa_enz,
                    kappa_TF_here = kappa_TF)
    # Calculating mRNA variation rates
    mRNA_GFP_prod = gamma_GFP * n_GFP * epsilon * x/(x + kappa_GFP) * K_tox/(K_tox + tox) * (R_mRNA/R_mRNA + K_R_mRNA)
    mRNA_TF_prod = gamma_TF * n_TF * x/(x + kappa_TF) * K_tox/(K_tox + tox) * (R_mRNA/R_mRNA + K_R_mRNA) 
    mRNA_enz_prod = gamma_enz * n_enz * x/(x + kappa_enz) * K_tox/(K_tox + tox) * (R_mRNA/R_mRNA + K_R_mRNA) 
    
    dmRNA_GFP = mRNA_GFP_prod - delta_GFP * mRNA_GFP
    dmRNA_TF = mRNA_TF_prod - delta_TF * mRNA_TF
    dmRNA_enz = mRNA_enz_prod - delta_enz * mRNA_enz
    
    # Solving for ribosomes at each time step
    y = solve_for_y(Y = Y_tot,
                    mRNA_GFP_here = mRNA_GFP, mRNA_enz_here = mRNA_enz, mRNA_TF_here = mRNA_TF,
                    k_GFP_here = k_GFP, k_enz_here = k_enz, k_TF_here = k_TF)
    
    dGFP = pi_GFP * mRNA_GFP * y/(y + k_GFP) * K_tox/(K_tox + tox)  - lambda_GFP * GFP
    dTF = pi_TF * mRNA_TF * y/(y + k_TF) * K_tox/(K_tox + tox)  - lambda_TF * TF
    denz = pi_enz * mRNA_enz * y/(y + k_enz) * K_tox/(K_tox + tox)  - lambda_enz * enz
    
    dtox = dGFP + denz + dTF
    dR_mRNA = -(720 * dGFP + mRNA_enz_length * denz + 954 * dTF)
    
    return(list(c(dmRNA_GFP, dmRNA_TF, dmRNA_enz, dGFP, dTF, denz, dbenzoic_acid, dinducer, dR_mRNA, dtox)))
  })
}

times <- seq(0, 360, by = 0.1)  # Simulate for 6 hours

# Total biosensor
name_for_saving = "biosensor_plate"
# Generating and visualising time data

benzoic_acid_TF_GFP_range = c(0,10,100,1000)
TF_range = c(0, 0.1, 0.3, 1, 3, 10, 30, 100)
GFP_range = c(0, 0.1, 0.3, 1, 3, 10, 30, 100)

TF_GFP_adder_df = NULL
colnames_time_TF_GFP_adder_df = c("Time")
name_for_saving = "biosensor"

for (benzoic_acid_concentration in benzoic_acid_TF_GFP_range) {
  for (TF_concentration in TF_range) {
    for (GFP_concentration in GFP_range) {
      init_solver["inducer"] = benzoic_acid_concentration
      all_parameters["n_copy_GFP"] = GFP_concentration
      all_parameters["n_copy_TF"] = TF_concentration
      # Solve the system
      out_mRNA_GFP <- ode(y = init_solver, times = times, func = solver_all_resource, parms = all_parameters)
      df_time <- data.frame(out_mRNA_GFP)
      # Select 4 hours
      value_here = df_time[which(df_time[,c("time")] == 240), c("GFP")]
      TF_GFP_adder_df = rbind(TF_GFP_adder_df, c("benzoic_acid" = benzoic_acid_concentration, "TF" =TF_concentration, "GFP" = GFP_concentration, 'TF_GFP' = value_here + 0.05))
    }
  }
}
colnames(TF_GFP_adder_df) = c("benzoic_acid", "TF", "GFP", "TF_GFP")
rownames(TF_GFP_adder_df) = NULL
TF_GFP_adder_df = data.frame(TF_GFP_adder_df)

folder_for_saving = ''  # Define path in your computer
write.csv(TF_GFP_adder_df, paste(folder_for_saving, "biosensor_TF.csv", sep = '/'), row.names=FALSE)

# For visualisation, use the attached Jupyter notebook (biosensor_visualisation.ipynb)
