# ------------------------------------------------------------------------------
# Script Name: TIP_dynamics_VL_dependent_production_rate.R
# Description: This script reproduces Supplementary Figure 2 from the manuscript
#              "Immune Responses May Make HIV-1 Therapeutic Interfering Particles Less Effective"
#              It simulates TIP dynamics in our basic model without an immune response and with a viral load-dependent
#.             target cell production rate (Eq. 7-14 in the text with the source term from the figure caption)
#
# Authors: Griffin Kutler Dodd, Rob J. de Boer
#
# Date Created: 2025-08-16
# Last Modified: 2025-08-16
#
# Requirements:
#   - Packages: ggplot2, patchwork
#   - The script grind.R should be in the same directory
#
# Usage:
#   Run the entire script to generate Supplementary Figures 2 from the paper.
#   Outputs will be saved in the '../Figures/' directory.
#
#
# ------------------------------------------------------------------------------

source("grind.R")
library(ggplot2)
library(patchwork)

#Define the model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dTt = lambda*(1+(m*V+m*V_M)/(h+V+V_M)) - d*Tt - b*V*Tt - b_prime*V_M*Tt - b_prime*V_T*Tt
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

#Create directory to store figures
dir.create("../Figures", showWarnings = FALSE)

###############################
#Generate Supplementary Figure 2
###############################

num_points_to_sample <- 62500
n_prime_vals <- seq(1000, 10000, length.out=as.integer(sqrt(num_points_to_sample)))
f_vals <- seq(0.01, 0.5, length.out=as.integer(sqrt(num_points_to_sample)))

#Simulate the ODEs and find the SS in the absence of TIPs
get_ss_viral_load <- function(h, s_init) {
  p <<- c(lambda=2000, d=0.2, b=4e-6, delta=1, n=1000, c=20, f=0, delta_prime=1, n_prime=0, b_prime=0, m=9, h=h)
  s <<- s_init
  s <<- run(500, timeplot=FALSE)
  newton(silent=TRUE)$state['V']
}

#Get SSs in the absence of TIPs for both values of h
ss_VL_without_TIPs_h1 <- get_ss_viral_load(1e4, c(Tt=10000, I=0, I_T=0, I_D=0, V=5, V_M=0, V_T=0, V_H=0))
ss_h1 <- newton(silent=TRUE)$state

ss_VL_without_TIPs_h2 <- get_ss_viral_load(1e6, c(Tt=10000, I=0, I_T=0, I_D=0, V=5, V_M=0, V_T=0, V_H=0))
ss_h2 <- newton(silent=TRUE)$state

#Function to generate each panel
simulate_panel <- function(h, s_baseline, baseline_VL) {
  ss_vals <- outer(n_prime_vals, f_vals, Vectorize(function(n_p, f) {
    p <<- c(lambda=2000, d=0.2, b=4e-6, delta=1, n=1000, c=20, f=f, 
            delta_prime=1, n_prime=n_p, b_prime=4e-6, m=9, h=h)
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
  list(h=1e4, s=ss_h1, vl=ss_VL_without_TIPs_h1),
  list(h=1e6, s=ss_h2, vl=ss_VL_without_TIPs_h2)
)

plots <- lapply(panel_settings, function(cfg) {
  simulate_panel(cfg$h, cfg$s, cfg$vl)
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

ggsave("../Figures/supfig2.png", plot=fig_all_labeled, device="png", width=5000, height=2300, units="px")