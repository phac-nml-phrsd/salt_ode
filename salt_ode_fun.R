

####################################################################################
##                                                                                ##
##  TITLE:    SH Broiler Transmission Model Functions                             ##
##  AUTHOR:   Brennan Chapman                                                     ##
##  CONTACT:  @chapb                                                              ##
##  DATE:     2018 - 2019                                                         ##
##                                                                                ##
####################################################################################


usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) install.packages(p, dep = TRUE)
  library(p, character.only = TRUE)
}


# Load Packages --------------------------------------------------------------------
# __________________________________________________________________________________

usePackage("deSolve")      # ODE Solver
usePackage("triangle")     # Probability Distributions (Triangular)
usePackage("mc2d")         # Probability Distributions (Pert)
usePackage("dplyr")        # Data Frame Manipulation
usePackage("tolerance")    # Tolerance Interval Creation
usePackage("magrittr")     # Pipe Operators
usePackage("tidyverse")    # Tidyverse
usePackage("reshape2")     # Manupulating Dataframes for ggplot2
usePackage("stringr")      # Manupulating Strings
usePackage("ggthemes")     # Adding Themes for ggplot2
usePackage("parallel")     # Use Parallel Computing



# Parallel Options -----------------------------------------------------------------
# __________________________________________________________________________________

useParallel     <- TRUE    # Use Parallel Computation to Speed Up Model

autoDetectCores <- TRUE    # Automatically Detect Cores on Machine and Use All Cores
setManualCores  <- 2       # Set Core Number if Auto-Detect Fails/Incorrect

detectedCores   <- parallel::detectCores()

if (autoDetectCores & !is.na(detectedCores)) {
  cores <- detectedCores
} else {
  cores <- setManualCores
}



# Model Functions ------------------------------------------------------------------
# __________________________________________________________________________________

compartment_model       <- function (t, x, params) {
  
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
    dD <-                                                                                                   omega * N
    dx <- c(dS, dE, dI, dR, dC, dD)
    list(dx)})
}

simulate_production     <- function(N, nday, E0, I0, R0, C0, D0, rho, eta, gamma, alpha, delta, lambda, omega) {
  
  # USE:
  #   Simulate the compartmental model.
  # ARG:
  #   Lists, of length nsim, containing initial conditions and parameters.
  #   Alternatively, if the initial conditions or parameters do not vary
  #   stochastically, a single value to be recycled.
  # OUT:
  #   A list, of length nsim, containing the results of each simulation of the
  #   model. The results themselves are the daily estimates of the counts in
  #   each compartment, for days 1:nday. These results are not easily
  #   human-readable - use process_results() to reformat the results as a data
  #   frame, and add additional computed columns.
  # SEE:
  #   Simulation of the compartmental model is encapsulated in a function to benefit
  #   from the speed and parallelization offered by the apply() family. The number
  #   of simulations is determined by the number of objects in the passed parameter
  #   lists; the number of simulations is equal to max(length(parameter1,
  #   parameter2, ...)), and arguments are recycled where necessary.
  #
  #   Unlike a traditional implementation of a compartmental model (wherein the
  #   solver function is called once to evaluate the system of ODEs at each time
  #   point), here the solver is called to evaluate the ODEs in a step-wise fashion,
  #   as to allow the parameters to vary on a daily basis. In this configuration,
  #   the results of the solver for timestep = b - 1, b are passed as the initial
  #   conditions for timestep = b, b + 1.
  
  # Create an Empty Data Frame for The Results of Each Simulation
  sim_result <- data.frame(time = 1:nday, S = 0, E = 0, I = 0, R = 0, C = 0, D = 0)
  
  # Create a Vector of Start Values for Each Simulation
  xstart     <- c(S = (N - E0 - I0 - R0), E = E0, I = I0, R = R0, C = C0, D = D0)
  
  for (b in 1:nday) {
    
    # Create a Two-Day Time Step (Previous Day, and Current)
    timestep <- c(b - 1, b)                    
    params   <- list(rho    = rho[b],
                     alpha  = alpha,
                     gamma  = gamma,
                     delta  = delta,
                     eta    = eta,
                     lambda = lambda,
                     omega  = omega[b],
                     N      = N)
    
    # Run the Solver, and Return a Data Frame with the Previous Days Current Day's Values
    dayout <- as.data.frame(deSolve::lsoda(xstart, timestep, compartment_model, params))
    
    # Fill the Data Frame with the Current Day's Values
    sim_result[b, ] <- dayout[2, ]
    
    # Set the Current Day's Values as the Starting Conditions for the Next Day's Calculation
    xstart <- as.numeric(dayout[2, 2:7])    # Drop 'time' in the First Position
    
  }
  
  round(sim_result)
  
}

run_simulations         <- function(init_conditions, ...) {
  
  # USE:
  #   Run nsim simulations of the model.
  # ARG:
  #   init_conditions:  data frame; the initial conditions to be used
  #                     for each simulation.
  # OUT:
  #   A list, of length nsim, of data frames containing simulation results.
  
  if (!missing(...)) {
    arglist <- list(...)
    for (i in 1:length(arglist)) {
      assign(names(arglist)[i], arglist[[i]])
    }
  }
  
  
  E0 <- init_conditions$E0
  I0 <- init_conditions$I0
  R0 <- init_conditions$R0
  C0 <- init_conditions$C0
  D0 <- init_conditions$D0
  
  if (useParallel) {
    
    cl  <- parallel::makeCluster(cores)
    clusterExport(cl, c("compartment_model"))
    
    res <- clusterMap(cl, simulate_production, N, nday, E0, I0, R0, C0, D0, 
                      param_rho, param_eta, param_gamma, param_alpha, 
                      param_delta, param_lambda, param_omega, SIMPLIFY = FALSE)
  
    stopCluster(cl)
    
    lapply(res, function (x) mutate(x, WFP = I/(S+E+I+R), Conc = ifelse(D > 0, C/D, 0)))
    
  } else {
    
    res <- mapply(simulate_production, N, nday, E0, I0, R0, C0, D0, 
                param_rho, param_eta, param_gamma, param_alpha, 
                param_delta, param_lambda, param_omega, SIMPLIFY = FALSE)
    
    lapply(res, function (x) mutate(x, WFP = I/(S+E+I+R), Conc = ifelse(D > 0, C/D, 0)))

  }
  
}



# Basic Functions ------------------------------------------------------------------
# __________________________________________________________________________________

hatch_contam            <- function(p_hatch, p_dolly, p_chick, sims, flks) {
  
  # USE:
  #   Tabulate the probability and number of contaminated/infected
  #   objects/birds at the hatchery for each simulation.
  # ARG:
  #   p_hatch:      vect; the probability a hatcher is contaminated.
  #   p_dolly:      vect; the probability a dolly is contaminated.
  #   p_chick:      vect; the probability a chick is contaminated/infected.
  #   sims:         int ; number of simulations.
  #   flks:         int ; flock size.
  # OUT:
  #   A data frame containing the probability and number of contaminated/infected
  #   objects/birds at the hatchery for each simulation.
  
  n_hatch       <- rbinom(sims, n_hatch, p_hatch)
  n_dolly       <- rbinom(sims, n_dolly, p_dolly)
  
  n_chick_dolly <- rpert(sims,
                         min   = cap_tray,
                         mode  = (p_chick * cap_dolly),
                         max   = cap_dolly,
                         shape = 4)
  
  n_chick       <- ceiling(n_hatch * n_dolly * n_chick_dolly)
  n_chick       <- ifelse(n_chick > flks, flks, n_chick)
  
  data.frame(p_hatch, p_dolly, p_chick,
             n_hatch, n_dolly, n_chick_dolly, n_chick)
  
} # UP

residual_load           <- function(prev_batch, days, sims) {
  
  # USE:
  #   Extract the load of Salmonella in the environment at day days (last day of
  #   rearing) from a burn-in model run.
  # ARG:
  #   prev_batch:   list; output of run_simulations() for a burn-in run.
  #   days:         int ; number of days simulated.
  #   sims:         int ; number of simulations.
  # OUT:
  #   A vector, of length sims, containing random samples from a pert
  #   distribution, defined by the mean and 95th percentiles of the load in the
  #   environment ("C"), on the final day of rearing in the passed results.
  
  cvals <- sapply(prev_batch, '[[', "C")
  csumm <- summary(cvals[days, ])
  rpert(sims, csumm[1], csumm[3], csumm[6])
} # UP

result_on_day           <- function(results, compartment, day) {
  
  # USE:
  #   Get the value of the compartment for each simulation on the given day.
  # ARG:
  #   results:      list; output of run_simulations().
  #   compartment:  str : the compartment to extract.
  #   day:          int ; the day to extract.
  # OUT:
  #   A vector of length nday containing the computed statistic.
  # EXM:
  #   result_on_day(res_s1, "S", 36)  
  
  sapply(results, '[[', day, compartment)
} # UP

compare_scenario        <- function(results, compartments, statistic, logT = FALSE, jitter = FALSE) {

  # USE:
  #   Compare the statistic for the compartment, across scenarios.
  # ARG:
  #   results:      list; a tibble::lst or named list of outputs of run_simulations().
  #   compartments: list; compartments to be plotted.
  #   statistic:    fun ; a statistical function.
  #   logT:         logi; log transform the results.
  #   jitter:       logi; include jitter in the plot.
  # OUT:
  #   Plot.

  if (is.null(names(results))) {
    print("WARNING: Non-named list of results. Plot legend will not include series names.")
  }
  
  subFun <- function(compartment) {
    
    x   <- as.data.frame(sapply(results, 
                                daily_stat, 
                                compartment = compartment, 
                                statistic = statistic))
    day   <- nrow(x)
    x     <- melt(x, variable.name = "Scenario", id.vars = NULL)
    x$day <- 1:day
    
    ti   <- compt_key(compartment, description = FALSE)
    ylab <- compt_key(compartment, description = TRUE)
    
    if(logT == TRUE) {
      x$value  %<>% log10()
      ylab     %<>% paste("(Log10)")
    }
    
    if (jitter == TRUE) {
      position <- "jitter"
    } else {
      position <- "identity"
    }
    
    print(ggplot(data = x, aes(x = day, y = value)) + 
            geom_line(aes(color = Scenario, group = Scenario), position = position) +
            labs(title = ti, x = "Time (Days)", y = ylab))
  }
  
  invisible(lapply(compartments, subFun))

} # UP

compare_init            <- function(inits, compartments, logT = FALSE, jitter = FALSE) {
  # USE:
  #   Compare the initial values for the compartments, between the scenarios.
  # ARG:
  #   inits  :      list; a tibble::lst or named list of initial conditions for scenarios.
  #   compartments: list; compartments to be plotted.
  #   logT:         logi; log transform the results.
  #   jitter:       logi; include jitter in the plot.
  # OUT:
  #   Plot printed to viewer.
  # EXM:
  #   compare_init("S", logT = TRUE, jitter = FALSE, init_s1, init_s2, init_s3)
  
  subFun <- function(compartment) {
    x     <- do.call(cbind, lapply(inits, '[[', compartment))
    x     <- melt(x)
    names(x) <- c("sim", "Scenario", "value")
    
    cname <- substr(compartment, 1, 1)
    
    xlab <- "Simulation"
    ylab <- compt_key(cname, char = TRUE)
    ti   <- paste(compt_key(cname, char = TRUE, description = FALSE), "at T = 0")
    
    if(logT == TRUE) {
      x$value  %<>% log10()
      ylab     %<>% paste("(Log10)")
    }
    
    if (jitter == TRUE) {
      position <- "jitter"
    } else {
      position <- "identity"
    }
    
    print(ggplot(data = x, aes(x = sim, y = value)) + 
          geom_point(aes(color = Scenario, group = Scenario), position = position) +
          labs(title = ti, x = xlab, y = ylab))
  
  }
  
  invisible(lapply(compartments, subFun))
  
} # UP

generate_intervals      <- function(results, compartment, day, interval, export = FALSE) {
  
  # USE:
  #   Generate the tolerance and confidence intervals for the 
  #   specified compartment, on the specified day. The interval 
  #   represents both the level of the confidence interval, and 
  #   the proportion of the simulations covered by the tolerance
  #   interval.
  # ARG:
  #   results:     a list containing the results of multiple scenarios, 
  #                processed by process_results().
  #   compartment: string; the compartment for which the values should be
  #                retrieved.
  #   day:         the desired day.
  #   interval:    the interval used for the confidence and tolerance 
  #                interval. The inverval must be specified as a decimal 
  #                e.g. 0.95, 0.97, 0.99.
  #   export:      logical; a logical, indicating whether the result should
  #                be exported as a CSV in the working directory.
  # OUT:
  #   A data frame with the mean, intervals, and standard deviation. Where
  #   export = TRUE, a CSV in the working directory.
  
  cinterval <- (1 - interval)/2 + interval
  
  res <- as.data.frame(lapply(results, result_on_day, compartment, day))
  
  out <- data.frame(scenario = character(),
                    ti.alpha = double(),
                    ti.prop  = double(),
                    mean     = double(),
                    ti.lower = double(),
                    ti.upper = double(),
                    ci.lower = double(),
                    ci.upper = double(),
                    sd       = double())
  
  for (i in 1:length(results)) {
    
    out[i,2:6] <- normtol.int(res[,i], alpha = (1 - interval), P = interval, side = 2)
    
    n   <- length(res[,i])
    avg <- mean(res[,i])
    std <- sd(res[,i])
    err <- qt(p = cinterval, df = (n - 1)) * std / sqrt(n)
    
    out[i, "ci.lower"] <- avg - err
    out[i, "ci.upper"] <- avg + err
    out[i, "sd"]       <- std
    
  }
  
  out$scenario <- names(results)
  
  if (export == TRUE) {
    write_csv(out, paste(compartment, deparse(day), "INT.csv", sep = ""))
  }
  
  return(out)
  
}



# Plotting Functions ---------------------------------------------------------------
# __________________________________________________________________________________

case_convert            <- function(text, case = "NULL") {
  
  out <- switch(case,
                title    = str_to_title(text),
                lower    = str_to_lower(text),
                upper    = str_to_upper(text),
                sentence = {text <- str_to_lower(text)
                            l    <- str_to_upper(substr(text, 1, 1))
                            sub("^[[:alpha:]]", replacement = l, text)},
                text)
  
  out <- str_replace(out, "[Cc][Ff][Uu]", "CFU")
  out <- str_replace(out, "salmonella", "Salmonella")
  
  return(out)
  
}

compt_key               <- function(compartment, case = "NULL", description = TRUE, 
                                    char = TRUE, includeLog = FALSE) {
  
  # USE:
  #   Provide the title or description of the compartment.
  # ARG:
  #   compartment:      the compartment of interest.
  #   case:             string; case-modification of the returned string,
  #                     see case_convert() for details.
  #   description:      logical; return a description where true, a title 
  #                     where false.
  #   char:             logical; indicates whether compartment is passed as
  #                     a character string (true) or an object (false).
  #   includeLog:       logical; adds "(Log10)" to the returned string.
  # OUT:
  #   A string containing the compartment description or title.
  
  if (char == TRUE) {
    v <- compartment
  } else {  # Ensures function operates correctly where nested
    c <- quote(substitute(compartment))
    v <- eval(c)
    for (i in rev(head(sys.frames(), -1L))) {
      c[[2]] <- v
      v      <- eval(c, i)
    }
    v <- deparse(substitute(v))
  }
  
  if (description == TRUE) {
    out <- switch(v, 
                  S    = "The Number of Susceptible Birds",
                  E    = "The Number of Exposed Birds",
                  I    = "The Number of Infected Birds",
                  R    = "The Number of Recovered Birds",
                  C    = "The Number of CFUs in the Environment",
                  D    = "The Number of Grams of Feces in the Environment",
                  WFP  = "The Within-flock Prevalence of Salmonellosis",
                  Conc = "The Concentration of Salmonella in Feces in the Environment (CFU/gram)",
                  v)
  } else {
    out <- switch(v, 
                  S    = "Susceptible",
                  E    = "Exposed",
                  I    = "Infected",
                  R    = "Recovered",
                  C    = "CFUs Environment",
                  D    = "Feces Environment",
                  WFP  = "Within-flock Prevalence",
                  Conc = "Concentration of Salmonella in Feces",
                  v)
  }
  
  if (includeLog == TRUE && v %in% c("C", "D", "Conc")) {
    out <- paste(out, "(Log10)")
  }
  
  case_convert(out, case)
  
}

param_key               <- function(parameter,   case = "NULL", description = TRUE, 
                                    char = FALSE, includeLog = FALSE) {
  
  # USE:
  #   Provide the title or description of the parameter.
  # ARG:
  #   compartment:      the parameter of interest.
  #   case:             string; case-modification of the returned string,
  #                     see case_convert() for details.
  #   description:      logical; return a description where true, a title 
  #                     where false.
  #   char:             logical; indicates whether parameter is passed as
  #                     a character string (true) or an object (false).
  #   includeLog:       logical; adds "(Log10)" to the returned string.
  # OUT:
  #   A string containing the parameter description or title.
  
  if (char == TRUE) {
    v <- parameter
  } else {  # Ensures function operates correctly where nested
    c <- quote(substitute(parameter))
    v <- eval(c)
    for (i in rev(head(sys.frames(), -1L))) {
      c[[2]] <- v
      v      <- eval(c, i)
    }
    v <- deparse(substitute(v))
  } 
  
  if (description == TRUE) {
    out <- switch(v, 
                  param_alpha  = "The rate at which broilers recover from infection (day^-1)",
                  param_gamma  = "The rate at which broilers become infected, after being exposed (day^-1)",
                  param_delta  = "The excretion rate of Salmonella per gram of feces (CFU/gram)",
                  param_eta    = "The volume of feces ingested per day (grams/bird-day)",
                  param_lambda = "The rate of bacterial decay in the environment (day^-1)",
                  param_rho    = "The probability a single organism causes illness",
                  param_omega  = "The rate of fecal excretion (grams/bird-day)",
                  v)
  } else {
    out <- switch(v, 
                  param_alpha  = "Alpha",
                  param_gamma  = "Gamma",
                  param_delta  = "Delta",
                  param_eta    = "Eta",
                  param_lambda = "Lambda",
                  param_rho    = "Rho",
                  param_omega  = "Omega",
                  v)
  }
  
  if (includeLog == TRUE) {
    out <- paste(out, "(Log10)")
  }
  
  case_convert(out, case)
  
}

scen_key                <- function(scenario,    case = "NULL", description = TRUE, 
                                    char = FALSE) {
  
  # USE:
  #   Provide the title or description of the scenario.
  # ARG:
  #   scenario:         the scenario of interest.
  #   case:             string; case-modification of the returned string,
  #                     see case_convert() for details.
  #   description:      logical; return a description where true, a title 
  #                     where false.
  #   char:             logical; indicates whether parameter is passed as
  #                     a character string (true) or an object (false).
  # OUT:
  #   A string containing the scenario description or title.
  
  
  if (char == TRUE) {
    v <- scenario
  } else {  # Ensures function operates correctly where nested
    c <- quote(substitute(scenario))
    v <- eval(c)
    for (i in rev(head(sys.frames(), -1L))) {
      c[[2]] <- v
      v      <- eval(c, i)
    }
    v <- deparse(substitute(v))
  }
  
  if (description == TRUE) {
    out <- switch(v,
                  res_burn = "SB: Model Burn-in",
                  res_base = "S0: Baseline",
                  res_s1   = "S1: Vaccination",
                  res_s2   = "S2: Disinfection",
                  res_s3   = "S3: Reduction of C0 by 10^-5",
                  res_s4   = "S4: Reduction of C0 by 10^-7",
                  res_s5   = "S5: Extended Down Time",
                  res_s6   = "S6: Reduced Salmonella Excretion",
                  v)
  } else {
    out <- switch(v,
                  res_burn = "Model Burn-in",
                  res_base = "Baseline",
                  res_s1   = "Scenario 1",
                  res_s2   = "Scenario 2",
                  res_s3   = "Scenario 3",
                  res_s4   = "Scenario 4",
                  res_s5   = "Scenario 5",
                  res_s6   = "Scenario 6",
                  v)
  }
  
  case_convert(out, case)
  
}



# Parameter Functions --------------------------------------------------------------
# __________________________________________________________________________________

piecewise               <- function(breakpoints, functions, nday) {
  
  # USE:
  #   Apply multiple functions at different intervals along the series 1:nday,
  #   guarding against nday << breakpoint.
  # ARG:
  #   breakpoints:  a vector containing the breakpoints (intervals) over which 
  #                 the corresponding functions should be applied.
  #   functions:    a vector containing the functions to be applied. By default, 
  #                 piecewise() passes a vector 1:nday as the arguement to the 
  #                 supplied functions.
  # OUT:
  #   A vector of length nday.
  # EXM:
  #   peicewise(c(3,7,nday), c(function(x) x, function(x) x^2, function(x) x^3))
  
  dmax    <- max(nday, breakpoints)
  days    <- 1:dmax
  output  <- vector("list", dmax)
  r_break <- rev(breakpoints)
  r_funct <- rev(functions)
  output  <- mapply(r_funct[[1]], days)
  
  for (i in 2:length(r_break)) {
    output[1:r_break[[i]]] <- mapply(r_funct[[i]], days[1:r_break[[i]]])
  }
  output[1:nday]
}

param_plot              <- function(parameter, sims) {
  
  # USE:
  #   Plot the value of the passed parameter for each simulation.
  # ARG:
  #   parameter:        parameter values to be plotted.
  #   sims:             the number of simulations.
  # OUT:
  #   A scatter plot of the passed values.
  
  ggplot(mapping = aes(x = 1:sims, y = parameter)) + 
    geom_point() +
    labs(title = param_key(parameter, description = FALSE), x = "Simulation", y = param_key(parameter))
  
}



# Temporal Object Functions --------------------------------------------------------
# __________________________________________________________________________________

daily_val               <- function(results, compartment, days, sims) {
  
  # USE:
  #   Reshape the list - data frame structure of model results into a data frame.
  # ARG:
  #   results:      list; output of run_simulations().
  #   compartment:  str ; the compartment of interest.
  #   days:         int ; number of days simulated.
  #   sims:         int ; number of simulations.
  # OUT:
  #   A data frame of model results, of days columns, bt sims rows.
  # EXM:
  #   daily_val(res_baseline, "S", 36, 10000)
  
  values <- data.frame(row.names = 1:sims)
  
  for (i in 1:days) {
    values <- cbind(values, as.vector(sapply(results, '[[', i, compartment)))
  }
  
  colnames(values) <- c(1:days)
  return(values)

} # UP

daily_val_param         <- function(parameter, days, sims) {
  
  # USE:
  #   Reshape the list-list structure of a parameter into a data frame.
  # ARG:
  #   parameter:    list: parameter that varies over time.
  #   days:         int ; number of days simulated.
  #   sims:         int ; number of simulations.
  # OUT:
  #   A data frame of parameter values, of days columns, by sims rows.
  # EXM:
  #   daily_val_param(param_alpha, 36, 10000)
  
  values <- data.frame()
  
  for (i in 1:sims) {
    values <- rbind(values, parameter[[i]])
  }
  
  colnames(values) <- c(1:days)
  return(values)
  
} # UP

daily_val_plot          <- function(object, compartment, days, sims, parameter = FALSE) {
  
  # USE:
  #   Plot the value of the passed obect over time (1:days) for each simulation.
  # ARG:
  #   object:       list; output of run_simulations() or a parameter object.
  #   compartment:  str ; the compartment of interest.
  #   days:         int ; number of days simulated.
  #   sims:         int ; number of simulations.
  #   parameter:    logi; a parameter object.
  # OUT:
  #   A line plot.
  
  if (parameter == TRUE) {
    yval        <- daily_val_param(object, days, sims)
    ylab        <- param_key(object)
    ti          <- param_key(object, description = FALSE)
  } else {
    yval        <- daily_val(object, compartment, days, sims)
    ylab        <- compt_key(compartment)
    ti          <- deparse(substitute(object))
  }
  
  yval_long     <- melt(yval, variable.name = "day")
  yval_long$sim <- 1:sims
  
  ggplot(data = yval_long, aes(x = day, y = value)) + 
    geom_line(aes(color = sim, group = sim)) +
    labs(title = ti, x = "Time (Days)", y = ylab)

} # UP

dynamics                <- function(results, compartments, days, sims, logT = FALSE, jitter = FALSE) {
  
  # USE:
  #   Plot the results for each of the compartments over time.
  # ARG:
  #   results:      list; output of run_simulations().
  #   compartments: list; compartments to be plotted.
  #   days:         int ; number of days simulated.
  #   sims:         int ; number of simulations.
  #   logT:         logi; log transform the results.
  #   jitter:       logi; include jitter in the plot.
  # OUT:
  #   Plot.
  
  clist_bio <- list("S", "E", "I", "R")
  clist_env <- list("C", "D", "WFP", "Conc")
  
  ylab <- "Birds/Value"
  xlab <- "Time (Days)"
  ti   <- scen_key(results)
  
  if (all(compartments %in% clist_bio)) {
    ylab <- "Birds"
  } else if (all(compartments %in% clist_env)) {
    ylab <- "Value"
  }
  
  res      <- do.call(rbind, results)
  res$sim  <- rep(1:sims, each = days)
  res_long <- melt(res, variable.name = "compartment", measure.vars = compartments)
  
  if(logT == TRUE) {
    res_long$value %<>% log10()
    ylab           %<>% paste("(Log10)")
  }
  
  if (jitter == TRUE) {
    position <- "jitter"
  } else {
    position <- "identity"
  }
  
  ggplot(data = res_long, aes(x = time, y = value)) +
    geom_line(aes(color = compartment, group = interaction(compartment, sim)), position = position) +
    labs(title = ti, x = xlab, y = ylab)

} # UP


daily_stat              <- function(results, compartment, statistic) {
  
  # USE:
  #   Compute the statistic for the compartment for each day in 1:nday.
  # ARG:
  #   results:      list; output of run_simulations().
  #   compartment:  str ; the compartment of interest.
  #   statistic:    fun ; a statistical function.
  # OUT:
  #   An unnamed vector, of length nday, of the computed statistic.
  # EXM:
  #   daily_stat(res_base, "S", mean)
  
  unname(apply(sapply(results, '[[', compartment), 1, statistic))
  
} # UP

daily_stat_param        <- function(parameter, statistic, days, sims) {
  
  # USE:
  #   Compute the statistic for the parameter for each day in 1:days.
  # ARG:
  #   parameter:    list: parameter that varies over time.
  #   statistic:    fun ; a statistical function.
  #   days:         int ; number of days simulated.
  #   sims:         int ; number of simulations.
  # OUT:
  #   An unnamed vector, of length days, of the computed statistic.
  # EXM:
  #   daily_stat_param(param_alpha, mean, 36, 10000)
  
  unname(apply(daily_val_param(parameter, days, sims), 2, statistic))
  
} # UP

daily_stat_plot         <- function(object, compartment, statistic, days, sims, parameter = FALSE) {
  
  # USE:
  #   Plot the statistic of the passed obect over time (1:days).
  # ARG:
  #   object:       list; output of run_simulations() or a parameter object.
  #   compartment:  str ; the compartment of interest.
  #   days:         int ; number of days simulated.
  #   sims:         int ; number of simulations.
  #   parameter:    logi; a parameter object.
  # OUT:
  #   A line plot.

  if (parameter == TRUE) {
    yval        <- daily_stat_param(object, statistic, days, sims)
    ylab        <- param_key(object)
    ti          <- paste("The", str_to_title(deparse(substitute(statistic))), "of", param_key(object, description = FALSE))
  } else {
    yval        <- daily_stat(object, compartment, statistic)
    ylab        <- compt_key(compartment)
    ti          <- paste("The", str_to_title(deparse(substitute(statistic))), "of", deparse(substitute(object)))
  }
  
  yval_long     <- data.frame("value" = yval, "day" = c(1:days))
  
  ggplot(data = yval_long, aes(x = day, y = value)) +
    geom_line() + 
    labs(title = ti, x = "Time (Days)", y = ylab)
  
} # UP

dynamics_stat           <- function(results, statistic, compartments, days, logT = FALSE, jitter = FALSE) {
  
  # USE:
  #   Plot the statistic for each of the compartments over time.
  # ARG:
  #   results:      list; output of run_simulations().
  #   statistic:    fun ; a statistical function.
  #   compartments: list; compartments to be plotted.
  #   days:         int ; number of days simulated.
  #   logT:         logi; log transform the results.
  #   jitter:       logi; include jitter in the plot.
  # OUT:
  #   Plot.
  
  clist_bio <- list("S", "E", "I", "R")
  clist_env <- list("C", "D", "WFP", "Conc")
  
  ylab <- "Birds/Value"
  xlab <- "Time (Days)"
  ti   <- scen_key(results)
  
  if (all(compartments %in% clist_bio)) {
    ylab <- "Birds"
  } else if (all(compartments %in% clist_env)) {
    ylab <- "Value"
  }
  
  compt <- mapply(daily_stat, 
                  results     = list(results), 
                  compartment = compartments, 
                  statistic   = list(statistic))
  
  colnames(compt) <- compartments
  compt         %<>% as.data.frame()
  compt_long      <- melt(compt, variable.name = "compartment")
  compt_long$day  <- 1:days
  
  if (logT == TRUE) {
    compt_long$value %<>% log10()
    ylab             %<>% paste("(Log10)", sep = " ")
  }
  
  if (jitter == TRUE) {
    position <- "jitter"
  } else {
    position <- "identity"
  }
  
  ggplot(data = compt_long, aes(x = day, y = value)) +
    geom_line(aes(color = compartment, group = compartment), position = position) +
    labs(title = ti, x = xlab, y = ylab)

} # UP






