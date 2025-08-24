# ------------------------------------------------------------------------------
# Script Name: TIP_dynamics_multiphase_killing.R
# Description: This script reproduces Supplementary Figure 4 from the manuscript
#              "Immune Responses May Make HIV-1 Therapeutic Interfering Particles Less Effective"
#              It simulates TIP dynamics in our immune response model adjusted to include both early
#              (eclipse phase) and late killing (see appendix B.3)
#
# Authors: Griffin Kutler Dodd, Rob J. de Boer
#
# Date Created: 2025-08-17
# Last Modified: 2025-08-24
#
# Requirements:
#   - Packages: ggplot2, patchwork, cowplot
#   - The script grind.R should be in the same directory
#
# Usage:
#   Run the entire script to generate Supplementary Figure 4 from the paper.
#   Outputs will be saved in the '../Figures/' directory.
#
#
# ------------------------------------------------------------------------------

source("grind.R")
library(ggplot2)
library(cowplot)

#Define the model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dTt = lambda - d*Tt - b*V*Tt - b_prime*V_M*Tt - b_prime*V_T*Tt
    dI_1 = b*V*Tt + b_prime*V_M*Tt - delta_1*I_1 - gamma*I_1 - k1*I_1*E - b_prime*V_T*I_1
    dI_2 = gamma*I_1 - delta_2*I_2 - k2*I_2*E
    dI_T1 = b_prime*V_T*Tt - d*I_T1 - gamma*I_T1 - k1*I_T1*E - b*V*I_T1 - b_prime*V_M*I_T1
    dI_T2 = gamma*I_T1 - d*I_T2 - b*V*I_T2 - b_prime*V_M*I_T2 - k2*I_T2*E
    dI_D1 = b*V*(I_T1+I_T2) + b_prime*V_M*(I_T1+I_T2) + b_prime*V_T*I_1 - delta_1*I_D1 - gamma*I_D1 - k1*I_D1*E
    dI_D2 = gamma*I_D1 - delta_2_prime*I_D2 - k2*I_D2*E
    dV = n*delta_2*I_2  - c*V
    dV_M = (f^2)*n_prime*delta_2_prime*I_D2 - c*V_M
    dV_T = ((1-f)^2)*n_prime*delta_2_prime*I_D2 - c*V_T
    dV_H = 2*f*(1-f)*n_prime*delta_2_prime*I_D2 - c*V_H
    dE = (p*(epsilon_1*(I_1+I_T1+I_D1)+epsilon_2*(I_2+I_T2+I_D2))*E)/(h+(epsilon_1*(I_1+I_T1+I_D1)+epsilon_2*(I_2+I_T2+I_D2)+E)) - d_E*E
    
    return(list(c(dTt, dI_1, dI_2, dI_T1, dI_T2, dI_D1, dI_D2, dV, dV_M, dV_T, dV_H, dE)))  
  }) 
}

#Create directory to store figures
dir.create("../Figures", showWarnings = FALSE)

###############################
#Generate Supplementary Figure 4
###############################

#Note that this is a time-intensive simulation; reducing num_points_to_sample will speed it up significantly
num_points_to_sample <- 62500
n_prime_vals <- seq(1000, 10000, length.out=as.integer(sqrt(num_points_to_sample)))
f_vals <- seq(0.01, 0.5, length.out=as.integer(sqrt(num_points_to_sample)))

#Simulate the ODEs and find the SS in the absence of TIPs
get_ss_viral_load <- function(k1, e1, k2, e2, s_init) {
  p <<- c(lambda=2000, d=0.2, b=1e-5, delta_1=0.5, delta_2=1, delta_2_prime=1, gamma=1, n=1000, c=20, f=0, delta_prime=1, n_prime=0, b_prime=0, k1=k1, k2=k2, p=1.1, h=1000, d_E=0.1, epsilon_1=e1, epsilon_2=e2)
  s <<- s_init
  s <<- run(500, timeplot=FALSE)
  newton(silent=TRUE)$state['V']
}

#Get SSs in the absence of TIPs for all four combinations of k1 and k2
ss_VL_without_TIPs_mild_onlyk1 <- get_ss_viral_load(k1=1e-4, e1=1, k2=0, e2=0, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_mild_onlyk1 <- newton(silent=TRUE)$state
ss_VL_without_TIPs_moderate_onlyk1 <- get_ss_viral_load(k1=1e-3, e1=1, k2=0, e2=0, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_moderate_onlyk1 <- newton(silent=TRUE)$state
ss_VL_without_TIPs_mild_onlyk2 <- get_ss_viral_load(k1=0, e1=0, k2=1.3e-4, e2=1, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_mild_onlyk2 <- newton(silent=TRUE)$state
ss_VL_without_TIPs_moderate_onlyk2 <- get_ss_viral_load(k1=0, e1=0, k2=5e-3, e2=1, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_moderate_onlyk2 <- newton(silent=TRUE)$state
ss_VL_without_TIPs_mild_k1andk2 <- get_ss_viral_load(k1=2.4e-5, e1=1, k2=2.4e-5, e2=1, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_mild_k1andk2 <- newton(silent=TRUE)$state
ss_VL_without_TIPs_moderate_k1andk2 <- get_ss_viral_load(k1=4e-4, e1=1, k2=4e-4, e2=1, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_moderate_k1andk2 <- newton(silent=TRUE)$state

#Function to generate each panel
simulate_panel <- function(k1, e1, k2, e2, s_baseline, baseline_VL) {
  ss_vals <- outer(n_prime_vals, f_vals, Vectorize(function(n_p, f) {
    p <<- c(lambda=2000, d=0.2, b=1e-5, delta_1=0.5, delta_2=1, delta_2_prime=1, gamma=1, n=1000,
            c=20, f=f, delta_prime=1, n_prime=n_p, b_prime=1e-5, k1=k1, k2=k2, p=1.1, h=1000, d_E=0.1, epsilon_1=e1, epsilon_2=e2)
    s <<- s_baseline
    s['V_T'] <<- 5
    #Note that SS finding can be troublesome due to strong oscillations, so we instead take the average over a large number of timesteps after a burn-in period
    s <<- run(1000, timeplot=FALSE)
    sim <- run(2500, table=TRUE, timeplot=FALSE)
    (mean(sim$V) + mean(sim$V_M))  # Infectious VL
  }))
  
  ss_vals <- pmin(((ss_vals - baseline_VL) / baseline_VL) * 100, 100)
  df <- expand.grid(x=n_prime_vals, y=f_vals)
  df$z <- as.vector(ss_vals)  
  
  ggplot(df, aes(x=x, y=y, fill=z)) + 
    geom_raster() + 
    scale_fill_gradientn(colors=c("blue", "lightblue", "white", "red", "darkred"), limits=c(-100, 100), na.value = "transparent") +
    scale_x_continuous(expand = c(0, 0), breaks = pretty(n_prime_vals, 5)) +
    scale_y_continuous(expand = c(0, 0), breaks = pretty(f_vals, 5)) +
    labs(x=NULL, y=NULL) +
    theme_minimal() +
    guides(fill="none")
}

panel_settings <- list(
  list(k1=1e-4, e1=1, k2=0, e2=0, s=ss_mild_onlyk1, vl=ss_VL_without_TIPs_mild_onlyk1),
  list(k1=0, e1=0, k2=1.3e-4, e2=1, s=ss_mild_onlyk2, vl=ss_VL_without_TIPs_mild_onlyk2),
  list(k1=2.4e-5, e1=1, k2=2.4e-5, e2=1, s=ss_mild_k1andk2, vl=ss_VL_without_TIPs_mild_k1andk2),
  list(k1=1e-3, e1=1, k2=0, e2=0, s=ss_moderate_onlyk1, vl=ss_VL_without_TIPs_moderate_onlyk1),
  list(k1=0, e1=0, k2=5e-3, e2=1, s=ss_moderate_onlyk2, vl=ss_VL_without_TIPs_moderate_onlyk2),
  list(k1=4e-4, e1=1, k2=4e-4, e2=1, s=ss_moderate_k1andk2, vl=ss_VL_without_TIPs_moderate_k1andk2)
)

plots <- lapply(panel_settings, function(cfg) {
  simulate_panel(cfg$k1, cfg$e1, cfg$k2, cfg$e2, cfg$s, cfg$vl)
})

#Plot figure and add axis labels
(fig_all <- (plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]]) + 
    plot_layout(nrow=2) & 
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 25, color = "black"),
          axis.ticks = element_line(size = 0.3),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
          text = element_text(family = "Arial"),
          plot.margin = margin(25,35,25,35)))

fig_all_labeled <- ggdraw(fig_all) +
  draw_label(bquote(paste("Virus particles produced by a doubly infected cell (",italic("n'"),")")),
             x = 0.5, y = 0.01, vjust = 0, size = 33) +
  draw_label(bquote(paste("Fraction of wild-type RNA in a doubly infected cell (",italic(f),")")),
             x = 0.02, y = 0.5, angle = 90, vjust = 1, size = 33) +
  draw_label("Early killing", x=0.185, y=0.98, size=30) +
  draw_label("Late killing", x=0.512, y=0.98, size=30) +
  draw_label("Constant killing", x=0.835, y=0.98, size=30) +
  draw_label("Mild immune response", x=0.98, y=0.75, angle=-90, size=30) +
  draw_label("Moderate immune response", x=0.98, y=0.26, angle=-90, size=30)
  
ggsave("../Figures/supfig4.png", fig_all_labeled, device="png", width=7500, height=5000, units="px")
