

####################################################################################
##                                                                                ##
##  TITLE:    SH Broiler Transmission Model Publication Figures                   ##
##  AUTHOR:   Brennan Chapman                                                     ##
##  CONTACT:  @chapb                                                              ##
##  DATE:     2018 - 2019                                                         ##
##                                                                                ##
####################################################################################

####################################################################################
##                                                                                ##
##  The following code is used to generate figures for publication. It may        ##
##  include less robust error handling or validation than other sections.         ##
##                                                                                ##
####################################################################################

source("salt_ode_fun.R")  # Include Functions


# Figure 1 -------------------------------------------------------------------------
## Biological Compartments Over Time in Baseline, w. Variability  

makeFigure1 <- function(results) {
  
  n <- length(results)
  day <- nrow(results[[1]])
  c_list <- c("S", "E", "I", "R")
  
  out <- data.frame()
  
  for (compt in c_list) {
  
    raw <- lapply(results, '[[', compt)
    raw <- do.call(rbind, raw)
    avg <- apply(raw, 2, mean)
    std <- apply(raw, 2, sd)
    err <- qt(p = 0.975, df = (n - 1)) * std / sqrt(n)
    cil <- avg - err
    ciu <- avg + err
    max <- apply(raw, 2, max)
    min <- apply(raw, 2, min)
    tol <- apply(raw, 2, tolerance::normtol.int, P = 0.95, side = 2)
    til <- unlist(lapply(tol, '[[', "2-sided.lower"))
    tiu <- unlist(lapply(tol, '[[', "2-sided.upper"))
    cpt <- as.data.frame(cbind(day = 1:36, avg, cil, ciu, min, max, til, tiu))
    cpt$Compartment <- compt
    out <- rbind(out, cpt)
  
  }
  
  return(out)
  
}

fig1data <- makeFigure1(res_base)
fig1data$Compartment <- factor(fig1data$Compartment, levels = c("S", "E", "I", "R"))

ggplot(data = fig1data, aes(x = day, y = avg, colour = Compartment)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = til, ymax = tiu, fill = Compartment), alpha = 0.1, linetype = 2) +
  labs(x = "Time (day)", y = "Number of birds") + 
  theme_classic() +
  scale_color_manual(values=c("#669900", "orange", "red", "blue")) + 
  scale_fill_manual(values=c("#669900", "orange", "red", "blue"))



# Figure 2 and 3 -------------------------------------------------------------------
## WFP and Conc Over Time, Across Scenarios

res_list      <- tibble::lst(res_base, res_s1, res_s2, res_s3, res_s4, res_s5, res_s6)
comp_scen_fig <- function(results, compartments, statistic, ylab, logT = FALSE, jitter = FALSE) {

  names(results) <- c("Baseline", "Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6")
  
  subFun <- function(compartment) {
    
    x   <- as.data.frame(sapply(results, 
                                daily_stat, 
                                compartment = compartment, 
                                statistic = statistic))
    day   <- nrow(x)
    x     <- melt(x, variable.name = "Scenario", id.vars = NULL)
    x$day <- 1:day
    
    if(logT == TRUE) {
      x$value  %<>% log10()
      ylab     %<>% paste("(Log10)")
    }
    
    if (jitter == TRUE) {
      position <- "jitter"
    } else {
      position <- "identity"
    }
    
    cbp1 <- c("#000000", "#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#0072B2")
    
    print(ggplot(data = x, aes(x = day, y = value, color = Scenario)) + 
            geom_line(position = position, size = 1, linejoin = "round") +
            labs(x = "Time (day)", y = ylab) +
            theme_classic() +
            scale_color_manual(values = cbp1))
    
  }
  
  invisible(lapply(compartments, subFun))
  
}

comp_scen_fig(res_list, tibble::lst("WFP"), mean, "Within-flock prevalence", logT = FALSE, jitter = FALSE)
comp_scen_fig(res_list, tibble::lst("Conc"), mean, "Bacterial concentration (CFU/g of feces)", logT = TRUE, jitter = FALSE)



# Figure 4 and 5 -------------------------------------------------------------------
## Sensitivity Analysis of WFP and Conc

sar_wfp <- pcc_wfp[[7]]
names(sar_wfp) <- c("original", "bias", "ste", "min", "max")

ggplot(data = sar_wfp) +
  geom_pointrange(aes(x = names(sa_bounds), y = original, ymin = min, ymax = max), size = 0.25) +
  ylim(-1, 1) +
  labs(y = "Partial correlation coefficient") +
  scale_x_discrete(labels = c(expression(alpha), expression(C[0]), expression(delta), expression(eta), 
                              expression(gamma), expression(I[0]), expression(lambda), expression(omega), 
                              expression(rho))) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "italic"))


sar_conc <- pcc_conc[[7]]
names(sar_conc) <- c("original", "bias", "ste", "min", "max")

ggplot(data = sar_conc) +
  geom_pointrange(aes(x = names(sa_bounds), y = original, ymin = min, ymax = max), size = 0.25) +
  ylim(-1, 1) +
  labs(y = "Partial correlation coefficient") +
  scale_x_discrete(labels = c(expression(alpha), expression(C[0]), expression(delta), expression(eta), 
                              expression(gamma), expression(I[0]), expression(lambda), expression(omega), 
                              expression(rho))) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "italic"))


