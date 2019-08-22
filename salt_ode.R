

####################################################################################
##                                                                                ##
##  TITLE:    SH Broiler Transmission Model                                       ##
##  AUTHOR:   Brennan Chapman                                                     ##
##  CONTACT:  @chapb                                                              ##
##  DATE:     2018 - 2019                                                         ##
##                                                                                ##
####################################################################################


set.seed(130315)          # John Snow's Birthday
source("salt_ode_fun.R")  # Include Functions



# Compartmental Model Pseudocode ---------------------------------------------------
# __________________________________________________________________________________

  # S   # Number Susceptible
  # E   # Number Exposed
  # I   # Number Infected
  # R   # Number Recovered
  # C   # Number CFU Environment
  # D   # Number Grams Feces Environment
  # 
  # dS <- -(1 - exp(-rho * C * eta / D)) * S
  # dE <-  (1 - exp(-rho * C * eta / D)) * S  -  gamma * E
  # dI <-                                        gamma * E  -  alpha * I
  # dR <-                                                      alpha * I
  # dC <-                                                                 delta * omega * I  -  lambda * C
  # dD <-                                                                                                   omega * N
    


# Model Parameters - General -------------------------------------------------------
# __________________________________________________________________________________

nsim    <- 25000            # Number of Iterations 
nday    <- 36               # Number of Days to Simulate
nbrd    <- 20000            # Flock Size
N       <- rep(nbrd, nsim)  # Flock Size Vector



# Model Parameters - Hatchery ------------------------------------------------------
# __________________________________________________________________________________

cap_hatch_opt   <- c(7040, 12960, 15120, 16800, 19200, 21120, 28800, 42240)    # Hatcher Capacities of 3 Main Brands Used in Canada
cap_dolly_opt   <- c(2160, 2520, 3640, 4160, 6000, 7040)                       # Dolly Capacities of 3 Main Brands Used in Canada
cap_tray_opt    <- c(  77,   88,  130,  144,  168)                             # Tray Capacities of 3 Main Brands Used in Canada

cap_hatch       <- sample(cap_hatch_opt, nsim, replace = TRUE)                 # Hatcher Capacity
cap_dolly       <- sample(cap_dolly_opt, nsim, replace = TRUE)                 # Dolly Capacity
cap_tray        <- sample(cap_tray_opt,  nsim, replace = TRUE)                 # Tray Capacity

n_hatch         <- ceiling(N / cap_hatch)                                      # Number of Hatchers Required
n_dolly         <- ceiling(N / cap_dolly)                                      # Number of Dollys Required

p_hatch_base    <- rbeta(nsim,  31, 1665)   
p_dolly_base    <- rbeta(nsim,  36,   88)        
p_chick_base    <- rbeta(nsim, 129, 3468)

hc_base         <- hatch_contam(p_hatch_base, p_dolly_base, p_chick_base, sims = nsim, flks = nbrd)



# Model Parameters - Stochastic - Define -------------------------------------------
# __________________________________________________________________________________

# GAMMA
# The rate at which broilers become infected, after being exposed (day^-1)
avg_latent  <- runif(nsim, 1, 2)             # Duration of Latent Period
param_gamma <- rexp(nsim, avg_latent)                                

for (i in 1:nsim) {
  while (param_gamma[i] > 1 || param_gamma[i] < 0.5) {
    param_gamma[i] <- rexp(1, avg_latent[i])
  }
}

# ALPHA
# The rate at which broilers recover from infection (day^-1)
avg_infect  <- rtriangle(nsim, 21, 28, 23)   # Duration of Infectious Period
param_alpha <- rexp(nsim, avg_infect)

for (i in 1:nsim) {
  while (param_alpha[i] > (1/21) || param_alpha[i] < (1/28)) {       
    param_alpha[i] <- rexp(1, avg_infect[i])  
  }
}

# OMEGA
# The rate of fecal excretion (grams/day)
dex_fec_A     <- function (x) 0.3057 *     x  + 0.7738
dex_fec_B     <- function (x) 9.4464 * log(x) - 4.1504
dex_fec       <- piecewise(breakpoints = c(3, nday),
                           functions = c(dex_fec_A, dex_fec_B),
                           nday = nday)
param_omega   <- vector("list", nsim)
param_omega[] <- list(dex_fec)

# DELTA
# The excretion rate of Salmonella per gram of feces (CFU/gram)
param_delta <- 10^(rpert(nsim, 1.7762, 3.8715, 5.9669))

# ETA
# The daily ingestion rate of feces per bird (grams/day)
param_eta   <- rpert(nsim, 0, 0.050, 2)

# RHO
# The probability a single organism causes illness (unitless dose-response parameter)
#rho_basic <- 0.21
rho_A     <- function (x) 0.00110976
rho_B     <- function (x) 0.000112252
rho_C     <- function (x) 1.07098*10^-09
day_rho   <- piecewise(breakpoints =  c(2, 7, nday),
                       functions =  c(rho_A, rho_B, rho_C),
                       nday = nday)
param_rho     <- vector("list", nsim)
param_rho[]   <- list(day_rho)

# LAMBDA
# The rate of bacterial decay in the environment (per day)
param_lambda <- rep(0.1, nsim)



# Model Parameters - Farm ----------------------------------------------------------
# __________________________________________________________________________________

base_downt <- rlogis(nsim, 16.4, 3.7)
base_downt <- ifelse(base_downt < 0, 0, base_downt)

base_eff_downt <- (1-param_lambda)^base_downt
base_eff_clean <- 10^(-(rpert(nsim,  0.8,  1.3,  1.8)))
base_eff_disin <- 10^(-(rpert(nsim, 2.22, 3.67, 4.26)))


delayedAssign("init_base",                      # Delay assignment until burn-in run is
              data.frame(E0 = rep(0, nsim),     # complete.
                         I0 = hc_base$n_chick,
                         R0 = rep(0, nsim),
                         C0 = residual_load(res_burn, nday, nsim) * base_eff_downt * base_eff_clean,
                         D0 = rep(1, nsim)))



# RESULTS --------------------------------------------------------------------------
# __________________________________________________________________________________

## Burn-in -----------------------------------------------------------------

p_hatch_burn    <- rep(1, nsim)                       # Force p_hatch = 1    
hc_burn         <- hatch_contam(p_hatch_burn, p_dolly_base, p_chick_base, sims = nsim, flks = nbrd)

init_burn    <- data.frame(E0 = rep(0, nsim), 
                           I0 = hc_burn$n_chick,
                           R0 = rep(0, nsim), 
                           C0 = rep(0, nsim), 
                           D0 = rep(1, nsim))

res_burn     <- run_simulations(init_burn)


## Baseline ----------------------------------------------------------------
## No Vaccination, No Disinfection

res_base     <- run_simulations(init_base)


## S1 - Vaccination --------------------------------------------------------

init_s1      <- init_base

# MODIFY
p_hatch_s1    <- rbeta(nsim,   2, 1409) 
p_dolly_s1    <- rbeta(nsim,  16,  113)       
p_chick_s1    <- rbeta(nsim,  16, 3724)
hc_s1         <- hatch_contam(p_hatch_s1, p_dolly_s1, p_chick_s1, sims = nsim, flks = nbrd)

init_s1$I0   <- hc_s1$n_chick

# RUN
res_s1       <- run_simulations(init_s1)


## S2 - Disinfection -------------------------------------------------------

init_s2      <- init_base

# MODIFY
init_s2$C0   <- init_s2$C0 * base_eff_disin

# RUN
res_s2       <- run_simulations(init_s2)


## S3 - Artificial C0 Reduction 10-5 ---------------------------------------

init_s3      <- init_base

# MODIFY
eff_reduc    <- 10^(-5)
init_s3$C0   <- init_s3$C0 * eff_reduc

# RUN
res_s3       <- run_simulations(init_s3)


## S4 - Artificial C0 Reduction 10-7 ---------------------------------------

init_s4      <- init_base

# MODIFY
eff_reduc    <- 10^(-7)
init_s4$C0   <- init_s4$C0 * eff_reduc

# RUN
res_s4       <- run_simulations(init_s4)


## S5 - Extended Down Time -------------------------------------------------

init_s5      <- init_base

# MODIFY
downt_s5     <- rlogis(nsim, 23.4, 3.7)
downt_s5     <- ifelse(downt_s5 < 0, 0, downt_s5)
eff_downt_s5 <- (1-param_lambda)^downt_s5
init_s5$C0   <- residual_load(res_burn, nday, nsim) * eff_downt_s5 * base_eff_clean

# RUN
res_s5       <- run_simulations(init_s5)


## S6 - Reduced Excretion of Salmonella ------------------------------------

init_burn_s6    <- init_burn
init_s6         <- init_base

# MODIFY
param_delta_s6  <- 10^(rpert(nsim, 0.7762, 2.8715, 4.9669))

# RUN BURN
res_burn_s6     <- run_simulations(init_burn_s6, param_delta = param_delta_s6)
  
# MODIFY
init_s6$C0      <- residual_load(res_burn_s6, nday, nsim) * base_eff_downt * base_eff_clean

# RUN
res_s6          <- run_simulations(init_s6, param_delta = param_delta_s6)



# Scenario Comparison --------------------------------------------------------------
# __________________________________________________________________________________

res_list  <- tibble::lst(res_base, res_s1, res_s2, res_s3, res_s4, res_s5, res_s6)
WFP36     <- generate_intervals(res_list,  "WFP", nday, 0.95, TRUE)
CONC36    <- generate_intervals(res_list, "Conc", nday, 0.95, TRUE)

