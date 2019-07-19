

####################################################################################
##                                                                                ##
##  TITLE:    SH Broiler Transmission Model Visualization                         ##
##  AUTHOR:   Brennan Chapman                                                     ##
##  CONTACT:  @chapb                                                              ##
##  DATE:     2018 - 2019                                                         ##
##                                                                                ##
####################################################################################


compt_list      <- list("S", "E", "I", "R", "C", "D", "WFP", "Conc")
compt_list_bio  <- compt_list[1:4]
compt_list_env  <- compt_list[5:8]
compt_list_init <- list("E0", "I0", "R0", "C0", "D0")
init_list       <- tibble::lst(init_base, init_s1, init_s2, init_s3, init_s4, init_s5, init_s6)
res_list        <- tibble::lst(res_base, res_s1, res_s2, res_s3, res_s4, res_s5, res_s6)



# Parameter Plots ------------------------------------------------------------------
# __________________________________________________________________________________

param_plot(param_gamma, nsim)
param_plot(param_alpha, nsim)
param_plot(param_lambda, nsim)
param_plot(param_delta, nsim)
param_plot(param_eta, nsim)

daily_val_plot(param_omega, NA, nday, nsim, parameter = TRUE)
daily_val_plot(param_rho, NA, nday, nsim, parameter = TRUE)



# Scenario Plots -------------------------------------------------------------------
# __________________________________________________________________________________


dynamics(res_burn, compt_list_bio, nday, nsim)
dynamics_stat(res_burn, mean, compt_list_bio, nday)
dynamics_stat(res_burn, mean, compt_list_env, nday, logT = TRUE)

dynamics(res_base, compt_list_bio, nday, nsim)
dynamics_stat(res_base, mean, compt_list_bio, nday)
dynamics_stat(res_base, mean, compt_list_env, nday, logT = TRUE)

dynamics(res_s1, compt_list_bio, nday, nsim)
dynamics_stat(res_s1, mean, compt_list_bio, nday)
dynamics_stat(res_s1, mean, compt_list_env, nday, logT = TRUE)

dynamics(res_s2, compt_list_bio, nday, nsim)
dynamics_stat(res_s2, mean, compt_list_bio, nday)
dynamics_stat(res_s2, mean, compt_list_env, nday, logT = TRUE)

dynamics(res_s3, compt_list_bio, nday, nsim)
dynamics_stat(res_s3, mean, compt_list_bio, nday)
dynamics_stat(res_s3, mean, compt_list_env, nday, logT = TRUE)

dynamics(res_s4, compt_list_bio, nday, nsim)
dynamics_stat(res_s4, mean, compt_list_bio, nday)
dynamics_stat(res_s4, mean, compt_list_env, nday, logT = TRUE)

dynamics(res_s5, compt_list_bio, nday, nsim)
dynamics_stat(res_s5, mean, compt_list_bio, nday)
dynamics_stat(res_s5, mean, compt_list_env, nday, logT = TRUE)

dynamics(res_s6, compt_list_bio, nday, nsim)
dynamics_stat(res_s6, mean, compt_list_bio, nday)
dynamics_stat(res_s6, mean, compt_list_env, nday, logT = TRUE)



# Scenario Comparison Plots --------------------------------------------------------
# __________________________________________________________________________________

compare_scenario(res_list, compt_list, mean)
compare_init(init_list, compt_list_init)


