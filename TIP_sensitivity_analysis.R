# ------------------------------------------------------------------------------
# Script Name: TIP_sensitivity_analysis.R
# Description: This script reproduces Supplementary Figure 3 from the manuscript
#              "Immune Responses May Make HIV-1 Therapeutic Interfering Particles Less Effective"
#              It generates parameter sets using Latin hypercube sampling and
#              performs global sensitivity analysis on the basic model using partial
#              rank correlation coefficients.
#
# Authors: Griffin Kutler Dodd, Rob J. de Boer
#
# Date Created: 2025-08-24
# Last Modified: 2025-08-24
#
# Requirements:
#   - Packages: ggplot2, lhs, sensitivity
#   - The script grind.R should be in the same directory
#
# Usage:
#   Run the entire script to generate Supplementary Figures 3 from the paper.
#   Outputs will be saved in the '../Figures/' directory.
#
#
# ------------------------------------------------------------------------------

dir.create("../Figures", showWarnings = FALSE)

###############################
#Generate Supplementary Figure 3
###############################

library(lhs)
library(sensitivity)
library(ggplot2)

source("grind.R")

#Define the model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dTt = lambda - d*Tt - b*V*Tt - b_prime*V_M*Tt - b_prime*V_T*Tt
    dI = b*V*Tt + b_prime*V_M*Tt - delta*I
    dI_T = b_prime*V_T*Tt - d*I_T - b*V*I_T - b_prime*V_M*I_T
    dI_D = b*V*I_T + b_prime*V_M*I_T - delta_prime*I_D
    dV = n*delta*I  - c*V
    dV_M = (f^2)*n_prime*delta_prime*I_D - c*V_M
    dV_T = ((1-f)^2)*n_prime*delta_prime*I_D - c*V_T
    dV_H = 2*f*(1-f)*n_prime*delta_prime*I_D - c*V_H
    
    return(list(c(dTt, dI, dI_T, dI_D, dV, dV_M, dV_T, dV_H)))  
  }) 
}

#Check constraints: R0>1; delta, delta' > d
constraint_fn <- function(params) {
  with(as.list(params), {
    R0 <- (b*n*lambda)/(20*d)
    return(R0>1 && delta>d && delta_prime>d)
  })
}

#Function to draw a desired number of LHS samples where R0>1
get_valid_lhs_samples <- function(n_req, param_ranges, oversample_factor = 5) {
  n_total <- n_req*oversample_factor
  repeat {
    #Draw samples and scale them to parameter ranges
    lhs_raw <- randomLHS(n_total, ncol(param_ranges))
    scaled <- lhs_raw
    for (i in 1:ncol(param_ranges)) {
      min_val <- param_ranges[1, i]
      max_val <- param_ranges[2, i]
      scaled[,i] <- min_val + lhs_raw[,i] * (max_val - min_val)
    }
    colnames(scaled) <- colnames(param_ranges)
    valid_rows <- apply(scaled, 1, constraint_fn)
    valid_samples <- scaled[valid_rows, , drop = FALSE]
    if (nrow(valid_samples) >= n_req) {
      return(as.data.frame(valid_samples[1:n_req, , drop = FALSE]))
    }
    
    #Double the amount of samples if we do not draw enough
    n_total <- n_total*2
  }
}

#Generate samples with R0>1
param_ranges <- data.frame(
  lambda = c(200, 20000),
  d = c(0.02, 2),
  b = c(4e-7, 4e-5),
  delta = c(0.1, 10),
  n = c(100, 10000),
  f = c(0, 0.5),
  delta_prime = c(0.1, 10),
  n_prime_n_ratio = c(1, 10),
  b_prime_b_ratio = c(0.1, 10)
)
num_samples <- 50000

samples <- get_valid_lhs_samples(num_samples, param_ranges)
samples$n_prime <- samples$n*samples$n_prime_n_ratio
samples$b_prime <- samples$b*samples$b_prime_b_ratio
samples$c <- 20

change_in_VL <- numeric(num_samples)

#For each parameter set, find the % change in VL
for (i in 1:num_samples) {
  pars <- samples[i,]
  p <<- as.list(pars)
  
  #Compute the VL in the absence of TIPs
  s <- c(Tt = as.numeric(pars['lambda'])/as.numeric(pars['d']),I=0,I_T=0,I_D=0,V=5,V_M=0,V_T=0,V_H=0)
  s <- run(150, timeplot=FALSE)
  s <- newton(silent=TRUE)$state
  baseline_VL <- s['V']
  
  #Compute the VL with TIPs
  s['V_T'] <- 5
  s <- run(150, timeplot=FALSE)
  s <- newton(silent=TRUE)$state
  treatment_VL <- s['V']+s['V_M']
  change_in_VL[i] <- ((treatment_VL-baseline_VL)/baseline_VL)*100
}

samples$change_in_VL <- change_in_VL

#Perform global sensitivity analysis
prcc_data <- pcc(X=samples[, names(param_ranges)], y=samples$change_in_VL, rank=TRUE)$PRCC
prcc_df <- data.frame(names(param_ranges), prcc_data)
colnames(prcc_df) <- c('Parameter', 'PRCC')

#Plot partial rank correlation coefficients
plot_prcc <- ggplot(prcc_df, aes(x = Parameter, y = PRCC, fill=PRCC)) +
  coord_flip() +
  scale_x_discrete(labels=c("lambda" = expression(lambda), 
                            "d" = expression(italic(d)), 
                            "b" = expression(beta),
                            "delta"=expression(delta),
                            "n"=expression(italic(n)),
                            "f"=expression(italic(f)),
                            "delta_prime"=expression(italic(delta * "'")),
                            "n_prime_n_ratio"=expression(italic(n * "'"  / n)),
                            "b_prime_b_ratio"=expression(italic(beta * "'"/ beta)))) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = round(PRCC, 3), hjust = ifelse(PRCC > 0, -0.2, 1.2)), size = 6, color = "black") +
  scale_fill_gradient2(low = "#4575b4", mid="white", high = "#d73027", midpoint=0, limits=c(-1,1)) +
  labs(x = "Parameter", y = "Partial Rank Correlation Coefficient") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 25, color = "black"),
    axis.ticks = element_line(size = 0.3),
    axis.title.x = element_text(size = 30, margin = margin(t = 10)),
    axis.title.y = element_text(size = 30),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    text = element_text(family = "Arial"),
    axis.text.y = element_text(margin = margin(r = 10)),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(color = "white"),
  ) +
  scale_y_continuous(limits=c(min(prcc_df$PRCC) - 0.1, max(prcc_df$PRCC) + 0.1), expand = c(0,0))

ggsave("../Figures/supfig3.png", plot=plot_prcc, device="png", width=5000, height=2300, units="px")
