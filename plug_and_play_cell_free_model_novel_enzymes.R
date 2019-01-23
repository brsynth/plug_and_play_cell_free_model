# The aim of this file is to generate the model heatmap data as presented in #add ref
# It generates data of the 2 enzymatic assays

#Author: Mathilde Koch, JL Faulon's research group, INRA

# Cleaning previous R workspace
rm(list = ls())

# Importing necessary packages: 

library(deSolve)  # For solving ODEs
library(ggplot2)  # For visualisation
library('reshape2')   # For visualisation


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



pars_HipO = c(
  kappa_enz = 3000,
  k_enz = 1e4,
  gamma_enz = 0.25,
  pi_enz = 3.6,
  omega_enz = 6,
  Km_enz = 764,
  k_cat_enz = 5880,
  spontaneous_hydrolisation = 0,
  mRNA_enz_length = 1200
)

pars_CocE = c(
  kappa_enz = 4500,
  k_enz = 3e4,
  gamma_enz = 0.176,
  pi_enz = 3.542,
  omega_enz = 8,
  Km_enz = 5.7,
  k_cat_enz = 3060,
  mRNA_enz_length = 1700,
  spontaneous_hydrolisation = 0.0001
)


# Choose HipO or CocE: parameters and name

all_parameters = c(global_parameters, pars_HipO)
# all_parameters = c(global_parameters, pars_CocE)
name = "Hippurate"
# name = "Cocaine"
colnames_for_visu = c("HipO", 'Hippurate', 'GFP')
# colnames_for_visu = c("CocE", 'Cocaine', 'GFP')

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
    rate_conversion = 0.001 * enz * k_cat_enz * inducer/(inducer + Km_enz)  
    # Rate is expressed In µM. * 0.001 since enzyme is in nM.
    dbenzoic_acid = rate_conversion + spontaneous_hydrolisation * inducer
    dinducer = - rate_conversion - spontaneous_hydrolisation * inducer
    
    activated_TF = TF * (benzoic_acid/(benzoic_acid + benzoic_kd) + 0.0005)
    epsilon = activated_TF/(activated_TF + activated_kd * 1000)
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

times <- seq(0, 360, by = 0.1)

inducer_range = c(0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
enzyme_range = c(0, 0.1, 0.3, 1, 3, 10, 30, 100)

# Solving the model for each inducer/enzyme DNA concentration

inducer_enzyme_adder_df = NULL
inducer_enzyme_adder_df_time = times
colnames_time_inducer_enzyme = c("Time")
name_for_saving = "inducer_enzyme_adder"
for (enz in enzyme_range) {
  for (inducer_concentration in inducer_range) {
    all_parameters["n_copy_enz"] = enz
    init_solver["inducer"] = inducer_concentration
    # Solving the system
    out_mRNA_GFP <- ode(y = init_solver, times = times, func = solver_all_resource, parms = all_parameters)
    # Time evolution of all species can be plotted by uncommenting the next line
    # plot(out_mRNA_GFP, mfrow = c(4, 3))
    df_time <- data.frame(out_mRNA_GFP)
    colnames_time_inducer_enzyme = c(colnames_time_inducer_enzyme, paste("enz", enz, "inducer", inducer_concentration, sep = '_'))
    inducer_enzyme_adder_df_time = cbind(inducer_enzyme_adder_df_time, "GFP" = df_time[, c("GFP")])
    # Getting value at 4 hours
    value_here = df_time[which(df_time[,c("time")] == 240), c("GFP")]
    inducer_enzyme_adder_df = rbind(inducer_enzyme_adder_df, c("enz" = enz,  "inducer" = inducer_concentration, 'GFP' = value_here + 0.05))
  }
}

# The model is solved, now is time for visualisation of the heatmap

saving_graph <- function(plot = last_plot(),filename, path = getwd(), device = "jpeg", 
                         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                         dpi = 300, limitsize = TRUE, reminder = FALSE) {
  
  ggsave(filename = filename, plot = plot, device = device, path = path, scale = scale, width = width, height = height, units = units,
         dpi = dpi, limitsize = limitsize)
  if (reminder == TRUE) {
    button <- tkmessageBox(title='Warning',message='Do not forget to save the data that produced that graph',type='ok')
    button <- tclvalue(button)
  }
}


plot_2D_heatmap = function(df_means, df_sd = NULL, title_base = "Visualising adder in 2D heatmap", name_for_saving, folder_for_concentration_images_strategy, saving = TRUE, x_lab = "(nM)", y_lab = "(µM)", subtitle_here = "GFP") {
  experiment_1 = colnames(df_means)[1]
  experiment_2 = colnames(df_means)[2]
  experiments = c(experiment_1, experiment_2)
  
  GFP_by_OD_factor <- melt(df_means, id.vars = experiments)
  GFP_by_OD_factor[,c(experiment_1)] <- as.factor(GFP_by_OD_factor[,c(experiment_1)])
  GFP_by_OD_factor[,c(experiment_2)] <- as.factor(GFP_by_OD_factor[,c(experiment_2)])
  
  xlab = paste(experiment_1, x_lab, sep = ' ')
  ylab = paste(experiment_2, y_lab, sep = ' ')
  
  print_heatmap(data_as_factor = GFP_by_OD_factor, x_lab = xlab, y_lab = ylab, x_axis = experiment_1, y_axis = experiment_2, title = title_base, subtitle = subtitle_here)
  if (saving) {  
    # folder_to_save <- paste(folder_for_concentration_images_strategy, "2D_heatmaps", sep = '/')
    # if (!file.exists(folder_to_save)) {
    #   dir.create(folder_to_save, showWarnings = FALSE)
    # }
    folder_to_save = folder_for_concentration_images_strategy
    saving_graph(filename = paste(name_for_saving, sep = '_'), reminder = FALSE, path = folder_to_save)
  }
}

print_heatmap <- function(data_as_factor, x_axis = "coc", y_axis = "hip", x_lab = "Compound 1", y_lab ="Compound 2", title = "", subtitle = "", 
                          print = TRUE) {
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  p <- ggplot(data = data_as_factor, aes(x= data_as_factor[,x_axis], y= data_as_factor[,y_axis], fill=value))
  p <- p + geom_tile()
  p <- p + scale_fill_gradientn("RFU", colours = jet.colors(7))

  p <- p + labs(x = x_lab, y = y_lab, title = title, subtitle = subtitle)
  p <- p + theme_bw() #theme black and white
  if (print) {
    print(p)
  } else {
    return(p)
  }
}

inducer_enzyme_adder_df = data.frame(inducer_enzyme_adder_df)

colnames(inducer_enzyme_adder_df) = colnames_for_visu
# Heatmap values can be saved by uncommenting the following line
write.csv(inducer_enzyme_adder_df, paste(name, '.csv', sep = ''))
plot_2D_heatmap(df_means= inducer_enzyme_adder_df, df_sd = NULL, subtitle_here = '',title_base = "", saving = TRUE, name_for_saving = "hip_hipO_adder", folder_for_concentration_images_strategy = "")
