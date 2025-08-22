# ------------------------------------------------------------------------------
# Script Name: TIP_dynamics_VL_eclipse_phase_immune_response.R
# Description: This script reproduces Figure 5 and Supplementary Figure 1 from the manuscript
#              "Immune Responses May Make HIV-1 Therapeutic Interfering Particles Less Effective"
#              It simulates TIP dynamics in our basic model with effector-cell killing only during the eclipse phase
#.             (Eq. 17-28 for Fig. 5, Eq. B.25-B.41 for Fig. S1)
#
# Authors: Griffin Kutler Dodd, Rob J. de Boer
#
# Date Created: 2025-08-16
# Last Modified: 2025-08-21
#
# Requirements:
#   - Packages: ggplot2, patchwork, cowplot
#   - The script grind.R should be in the same directory
#
# Usage:
#   Run the entire script to generate Fiugre 5 and Supplementary Figure 1 from the paper.
#   Outputs will be saved in the '../Figures/' directory.
#
#
# ------------------------------------------------------------------------------

source("grind.R")
library(ggplot2)
library(patchwork)
library(cowplot)

#Create directory to store figures
dir.create("../Figures", showWarnings = FALSE)

###############################
#Generate Figure 5
###############################

#Define the model (uniform killing rates)
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dTt = lambda - d*Tt - b*V*Tt - b_prime*V_M*Tt - b_prime*V_T*Tt
    dI_1 = b*V*Tt + b_prime*V_M*Tt - delta_1*I_1 - gamma*I_1 - k*I_1*E - b_prime*V_T*I_1
    dI_2 = gamma*I_1 - delta_2*I_2
    dI_T1 = b_prime*V_T*Tt - d*I_T1 - gamma*I_T1 - k*I_T1*E - b*V*I_T1 - b_prime*V_M*I_T1
    dI_T2 = gamma*I_T1 - d*I_T2 - b*V*I_T2 - b_prime*V_M*I_T2
    dI_D1 = b*V*(I_T1+I_T2) + b_prime*V_M*(I_T1+I_T2) + b_prime*V_T*I_1 - delta_1*I_D1 - gamma*I_D1 - k*I_D1*E
    dI_D2 = gamma*I_D1 - delta_2_prime*I_D2
    dV = n*delta_2*I_2  - c*V
    dV_M = (f^2)*n_prime*delta_2_prime*I_D2 - c*V_M
    dV_T = ((1-f)^2)*n_prime*delta_2_prime*I_D2 - c*V_T
    dV_H = 2*f*(1-f)*n_prime*delta_2_prime*I_D2 - c*V_H
    dE = (p*(I_1+I_T1+I_D1)*E)/(h+I_1+I_T1+I_D1+E) - d_E*E
    
    return(list(c(dTt, dI_1, dI_2, dI_T1, dI_T2, dI_D1, dI_D2, dV, dV_M, dV_T, dV_H, dE)))  
  }) 
}

num_points_to_sample <- 62500
n_prime_vals <- seq(1000, 10000, length.out=as.integer(sqrt(num_points_to_sample)))
f_vals <- seq(0.01, 0.5, length.out=as.integer(sqrt(num_points_to_sample)))

#Simulate the ODEs and find the SS in the absence of TIPs
get_ss_viral_load <- function(k, s_init) {
  p <<- c(lambda=2000, d=0.2, b=1e-5, delta_1=0.5, delta_2=1, delta_2_prime=1, gamma=1, n=1000, c=20, f=0, delta_prime=1, n_prime=0, b_prime=0, k=k, p=1.1, h=1000, d_E=0.1)
  s <<- s_init
  s <<- run(500, timeplot=FALSE)
  newton(silent=TRUE)$state['V']
}

#Get SSs in the absence of TIPs for both values of k
ss_VL_without_TIPs_k_low <- get_ss_viral_load(1e-4, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_k_low <- newton(silent=TRUE)$state
ss_VL_without_TIPs_k_high <- get_ss_viral_load(1e-3, c(Tt=10000, I_1=0, I_2=0, I_T1=0, I_T2=0, I_D1=0, I_D2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_k_high <- newton(silent=TRUE)$state

#Function to generate each panel
simulate_panel <- function(k, s_baseline, baseline_VL) {
  ss_vals <- outer(n_prime_vals, f_vals, Vectorize(function(n_p, f) {
    p <<- c(lambda=2000, d=0.2, b=1e-5, delta_1=0.5, delta_2=1, delta_2_prime=1, gamma=1, n=1000,
            c=20, f=f, delta_prime=1, n_prime=n_p, b_prime=1e-5, k=k, p=1.1, h=1000, d_E=0.1)
    s <<- s_baseline
    s['V_T'] <<- 5
    s <<- run(500, timeplot=FALSE)
    ss <- newton(silent=TRUE)
    (ss$state['V'] + ss$state['V_M'])  # Infectious VL
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
  list(k=1e-4, s=ss_k_low, vl=ss_VL_without_TIPs_k_low),
  list(k=1e-3, s=ss_k_high, vl=ss_VL_without_TIPs_k_high)
)

plots <- lapply(panel_settings, function(cfg) {
  simulate_panel(cfg$k, cfg$s, cfg$vl)
})

#Assemble and save final figure
(fig_all <- (plots[[1]] + plots[[2]]) + 
    plot_layout(nrow=1) & 
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 25, color = "black"),
          axis.ticks = element_line(size = 0.3),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
          text = element_text(family = "Arial"),
          plot.margin = margin(25,25,25,35)))

fig_all_labeled <- ggdraw(fig_all) +
  draw_plot_label(
    label = c("a", "b"),
    x = c(0.005, 0.51),
    y = c(0.98, 0.98),
    hjust = 0, vjust = 1,
    size = 35, fontface = "bold"
  ) +
  draw_label(bquote(paste("Virus particles produced by a doubly infected cell (",italic("n'"),")")),
             x = 0.5, y = 0.02, vjust = 0, size = 25) +
  draw_label(bquote(paste("Fraction of wild-type RNA in a doubly infected cell (",italic(f),")")),
             x = 0.04, y = 0.5, angle = 90, vjust = 1, size = 22)

ggsave("../Figures/fig5.png", plot=fig_all_labeled, device="png", width=5000, height=2300, units="px")

###############################
#Generate Supplementary Figure 1
###############################

#Model with explicit tracking of most recent infections and virus-specific killing rates (k, k_M, k_T)
#e_TIP is a boolean that indicates whether TIP-infected cells are killed (for the saturation function of effector cell proliferation)
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dTt = lambda - d*Tt - b*V*Tt - b_prime*V_M*Tt - b_prime*V_T*Tt
    dI_1 = b*V*Tt - delta_1*I_1 - gamma*I_1 - k*I_1*E - b_prime*V_T*I_1
    dI_2 = gamma*I_1 - delta_2*I_2
    dI_M1 = b_prime*V_M*Tt - delta_1*I_M1 - gamma*I_M1 - k_M*I_M1*E - b_prime*V_T*I_M1
    dI_M2 = gamma*I_M1 - delta_2*I_M2
    dI_T1 = b_prime*V_T*Tt - d*I_T1 - gamma*I_T1 - k_T*I_T1*E - b*V*I_T1 - b_prime*V_M*I_T1
    dI_T2 = gamma*I_T1 - d*I_T2 - b*V*I_T2 - b_prime*V_M*I_T2
    dI_TV1 = b*V*(I_T1+I_T2) - delta_1*I_TV1 - gamma*I_TV1 - k*I_TV1*E
    dI_TV2 = gamma*I_TV1 - delta_2_prime*I_TV2
    dI_VT1 = b_prime*V_T*I_1 - delta_1*I_VT1 - gamma*I_VT1 - k_T*I_VT1*E
    dI_VT2 = gamma*I_VT1 - delta_2_prime*I_VT2
    dI_TM1 = b_prime*V_M*(I_T1+I_T2) - delta_1*I_TM1 - gamma*I_TM1 - k_M*I_TM1*E
    dI_TM2 = gamma*I_TM1 - delta_2_prime*I_TM2
    dI_MT1 = b_prime*V_T*I_M1 - delta_1*I_MT1 - gamma*I_MT1 - k_T*I_MT1*E
    dI_MT2 = gamma*I_MT1 - delta_2_prime*I_MT2
    dV = n*delta_2*I_2  - c*V
    dV_M = (f^2)*n_prime*delta_2_prime*(I_TV2+I_VT2+I_TM2+I_MT2) - c*V_M
    dV_T = ((1-f)^2)*n_prime*delta_2_prime*(I_TV2+I_VT2+I_TM2+I_MT2) - c*V_T
    dV_H = 2*f*(1-f)*n_prime*delta_2_prime*(I_TV2+I_VT2+I_TM2+I_MT2) - c*V_H
    dE = (p*(I_1+I_M1+I_TV1+I_TM1+e_TIP*(I_T1+I_VT1+I_MT1))*E)/(h+I_1+I_M1+I_TV1+I_TM1+e_TIP*(I_T1+I_VT1+I_MT1)+E) - d_E*E
    
    return(list(c(dTt, dI_1, dI_2, dI_M1, dI_M2, dI_T1, dI_T2, dI_TV1, dI_TV2, dI_VT1, dI_VT2, dI_TM1, dI_TM2, dI_MT1, dI_MT2, dV, dV_M, dV_T, dV_H, dE)))  
  }) 
}

#Note that this is a time-intensive simulation; reducing num_points_to_sample will speed it up significantly
num_points_to_sample <- 62500
n_prime_vals <- seq(1000, 10000, length.out=as.integer(sqrt(num_points_to_sample)))
f_vals <- seq(0.01, 0.5, length.out=as.integer(sqrt(num_points_to_sample)))

get_ss_viral_load <- function(k, s_init) {
  p <<- c(lambda=2000, d=0.2, b=1e-5, delta_1=0.5, delta_2=1, delta_2_prime=1, gamma=1, n=1000, c=20, f=0, delta_prime=1, n_prime=0, b_prime=0, k=k, k_M=0, k_T=0, p=1.1, h=1000, d_E=0.1, e_TIP=0)
  s <<- s_init
  s <<- run(500, timeplot=FALSE)
  newton(silent=TRUE)$state['V']
}

ss_VL_without_TIPs_k_low <- get_ss_viral_load(1e-4, c(Tt=10000, I_1=0, I_2=0, I_M1=0, I_M2=0, I_T1=0, I_T2=0, I_TV1=0, I_TV2=0, I_VT1=0, I_VT2=0, I_TM1=0, I_TM2=0, I_MT1=0, I_MT2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_k_low <- newton(silent=TRUE)$state
ss_VL_without_TIPs_k_high <- get_ss_viral_load(1e-3, c(Tt=10000, I_1=0, I_2=0, I_M1=0, I_M2=0, I_T1=0, I_T2=0, I_TV1=0, I_TV2=0, I_VT1=0, I_VT2=0, I_TM1=0, I_TM2=0, I_MT1=0, I_MT2=0, V=5, V_M=0, V_T=0, V_H=0, E=1))
ss_k_high <- newton(silent=TRUE)$state

#Function to generate each panel
simulate_panel <- function(k, k_M, k_T, e_TIP, s_baseline, baseline_VL) {
  ss_vals <- outer(n_prime_vals, f_vals, Vectorize(function(n_p, f) {
    p <<- c(lambda=2000, d=0.2, b=1e-5, delta_1=0.5, delta_2=1, delta_2_prime=1, gamma=1, n=1000,
            c=20, f=f, delta_prime=1, n_prime=n_p, b_prime=1e-5, k=k, k_M=k_M, k_T=k_T, p=1.1, h=1000, d_E=0.1, e_TIP=e_TIP)
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
  list(k=1e-4, k_M=1e-3, k_T=1e-3, e_TIP=1, s=ss_k_low, vl=ss_VL_without_TIPs_k_low),
  list(k=1e-3, k_M=1e-3, k_T=0, e_TIP=0, s=ss_k_high, vl=ss_VL_without_TIPs_k_high)
)

plots <- lapply(panel_settings, function(cfg) {
  simulate_panel(cfg$k, cfg$k_M, cfg$k_T, cfg$e_TIP, cfg$s, cfg$vl)
})

#Assemble and save final figure
(fig_all <- (plots[[1]] + plots[[2]]) + 
    plot_layout(nrow=1) & 
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 25, color = "black"),
          axis.ticks = element_line(size = 0.3),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
          text = element_text(family = "Arial"),
          plot.margin = margin(25,25,25,35)))

fig_all_labeled <- ggdraw(fig_all) +
  draw_plot_label(
    label = c("a", "b"),
    x = c(0.005, 0.51),
    y = c(0.98, 0.98),
    hjust = 0, vjust = 1,
    size = 35, fontface = "bold"
  ) +
  draw_label(bquote(paste("Virus particles produced by a doubly infected cell (",italic("n'"),")")),
             x = 0.5, y = 0.02, vjust = 0, size = 25) +
  draw_label(bquote(paste("Fraction of wild-type RNA in a doubly infected cell (",italic(f),")")),
             x = 0.04, y = 0.5, angle = 90, vjust = 1, size = 22)

ggsave("../Figures/supfig1.png", plot=fig_all_labeled, device="png", width=5000, height=2300, units="px")