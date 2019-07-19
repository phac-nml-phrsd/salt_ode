

####################################################################################
##                                                                                ##
##  TITLE:    SH Broiler Transmission Model Sensitivity Analysis                  ##
##  AUTHOR:   Brennan Chapman                                                     ##
##  CONTACT:  @chapb                                                              ##
##  DATE:     2018 - 2019                                                         ##
##                                                                                ##
####################################################################################


set.seed(130315)          # John Snow's Birthday
source("salt_ode_fun.R")


usePackage("sensitivity")
usePackage("ODEsensitivity")
usePackage("lhs")


nbrdsa <- 20000
ndaysa <- 36
nsimsa <- 25000
NSA    <- rep(nbrdsa, nsimsa)



# PCC Sensitivity Analysis ---------------------------------------------------------
# __________________________________________________________________________________

useLHS <- FALSE

## SPECIFY BOUNDS
sa_bounds <- tibble::lst(lambda = c(0.05  , 0.15  ),
                         alpha  = c(1/28  , 1/21  ),
                         gamma  = c(0.5   , 1     ),
                         delta  = c(1     , 10^7  ),
                         omega  = c(0     , 40    ),
                         eta    = c(10^-5 , 10^-3 ),
                         rho    = c(0     , 0.4   ), 
                         I0     = c(0     , nbrdsa), 
                         C0     = c(0     , 10^10 ))

## GENERATE SAMPLE
if (useLHS) {
  sa_params <- as.data.frame(maximinLHS(n = nsimsa, k = length(sa_bounds)))
  
  for (i in 1:length(sa_bounds)) {
    sa_params[,i] <- sa_params[,i] * (sa_bounds[[i]][2] - sa_bounds[[i]][1]) + sa_bounds[[i]][1]
  }
  
  colnames(sa_params) <- names(sa_bounds)
  
} else {
  sa_params <- as.data.frame(do.call(cbind, lapply(sa_bounds, function(x) runif(nsimsa, x[1], x[2]))))
}

## EXPAND TEMPORAL PARAMETERS
sa_omega         <- vector("list", nsimsa)
sa_rho           <- vector("list", nsimsa)

for (i in 1:nsimsa) {
  sa_omega[[i]] <- (rep(sa_params$omega[[i]], ndaysa))
}

for (i in 1:nsimsa) {
  sa_rho_basic <- sa_params$rho[[i]]
  sa_rho_A     <- function (x) sa_rho_basic
  sa_rho_B     <- function (x) sa_rho_basic * 10^-1
  sa_rho_C     <- function (x) sa_rho_basic * 10^-6
  sa_day_rho   <- piecewise(breakpoints =  c(2, 7, ndaysa),
                            functions =  c(sa_rho_A, sa_rho_B, sa_rho_C),
                            nday = ndaysa)
  sa_rho[[i]]  <- sa_day_rho
}

## GENERATE INITIAL CONDITIONS
init_sa <- data.frame(E0 = rep(0, nsimsa),
                      I0 = sa_params$I0,
                      R0 = rep(0, nsimsa),
                      C0 = sa_params$C0,
                      D0 = rep(1, nsimsa))


## RUN MODEL
res_sa  <- run_simulations(init_conditions =  init_sa,
                           N               = NSA,
                           nday            = ndaysa,
                           param_rho       = sa_rho,
                           param_eta       = sa_params$eta,
                           param_gamma     = sa_params$gamma,
                           param_alpha     = sa_params$alpha,
                           param_delta     = sa_params$delta,
                           param_lambda    = sa_params$lambda,
                           param_omega     = sa_omega) 


## VISUALIZE RESULTS
#dynamics(res_sa, bio = TRUE)

compt_list      <- tibble::lst("S", "E", "I", "R", "C", "D", "WFP", "Conc")
compt_list_bio  <- compt_list[1:4]
compt_list_env  <- compt_list[5:8]


## SENSITIVITY ANALYSIS
sa_wfp    <- result_on_day(res_sa, "WFP", ndaysa)
sa_conc   <- result_on_day(res_sa, "Conc", ndaysa)
pcc_wfp  <- pcc(sa_params, sa_wfp, nboot = 100, rank = TRUE)
pcc_conc <- pcc(sa_params, sa_conc, nboot = 100, rank = TRUE)

#plot(pcc_wfp, ylim = c(-1, 1))
#plot(pcc_conc, ylim = c(-1, 1))



# SOBOL Sensitivity Analysis -------------------------------------------------------
# __________________________________________________________________________________

# A new copy of the compartmental is specified here to easily pass "N" (as
# nbrdsa) in the calculation of dD.

sa_compartment_model    <- function (t, x, params) {
  
  S <- x[1]    # Number Susceptible
  E <- x[2]    # Number Exposed
  I <- x[3]    # Number Infected
  R <- x[4]    # Number Recovered
  C <- x[5]    # Number CFU Environment
  D <- x[6]    # Number Grams Feces Environment
  
  with(as.list(params), {
    dS <- -(1 - exp(-rho * C * eta / D)) * S
    dE <-  (1 - exp(-rho * C * eta / D)) * S  -  gamma * E
    dI <-                                        gamma * E  -  alpha * I
    dR <-                                                      alpha * I
    dC <-                                                                 delta * omega * I  -  lambda * C
    dD <-                                                                                                   omega * nbrdsa
    dx <- c(dS, dE, dI, dR, dC, dD)
    list(dx)})
}

sas_mins <- sapply(sa_bounds[1:7], `[[`, 1)
sas_maxs <- sapply(sa_bounds[1:7], `[[`, 2)
sas_pars <- names(sa_bounds[1:7])

sas_I0 <- 1
sas_C0 <- mean(sa_bounds$C0)

sas_init <- c(S = nbrdsa - sas_I0, E = 0, I = sas_I0, R = 0, C = sas_C0, D = 1)

sa_sobol <- ODEsobol(mod        = sa_compartment_model,
                     pars       = sas_pars, 
                     state_init = sas_init,
                     times      = ndaysa,
                     n          = 25000,
                     rargs      = paste0("min = ", sas_mins, ", max = ", sas_maxs))

sobol_res <- data.frame(variable = sas_pars, value = sa_sobol$I$T[2:(length(sas_pars) + 1)])

#ggplot(data = sobol_res, aes(variable, value)) + geom_bar(stat = "identity")



