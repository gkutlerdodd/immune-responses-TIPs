# ------------------------------------------------------------------------------
# Script Name: TIP_dynamics_without_immune_response.R
# Description: This script reproduces Figures 2, 3, and 4 from the manuscript
#              "Immune Responses May Make HIV-1 Therapeutic Interfering Particles Less Effective"
#              It simulates TIP dynamics in our basic model without an immune response (Eq. 7-14 in the text)
#
# Authors: Griffin Kutler Dodd, Rob J. de Boer
#
# Date Created: 2025-08-15
# Last Modified: 2025-08-22
#
# Requirements:
#   - Packages: ggplot2, patchwork, cowplot, scales, grid, gridExtra
#   - The script grind.R should be in the same directory
#
# Usage:
#   Run the entire script to generate Figures 2, 3, and 4 from the paper.
#   Outputs will be saved in the '../Figures/' directory.
#
#
# ------------------------------------------------------------------------------

source("grind.R")
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)
library(grid)
library(gridExtra)

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

#Create directory to store figures
dir.create("../Figures", showWarnings = FALSE)

###############################
#Generate Figure 2
###############################

#Parameters
num_points_to_sample <- 62500
n_prime_vals <- seq(1000, 10000, length.out=as.integer(sqrt(num_points_to_sample)))
f_vals <- seq(0.01, 0.5, length.out=as.integer(sqrt(num_points_to_sample)))

#Simulate the ODEs and find the SS in the absence of TIPs
get_ss_viral_load <- function(beta, s_init) {
  p <<- c(lambda=2000, d=0.2, b=beta, delta=1, n=1000, c=20, f=0, delta_prime=1, n_prime=0, b_prime=0)
  s <<- s_init
  s <<- run(500, timeplot=FALSE)
  newton(silent=TRUE)$state['V']
}

#Get SSs in the absence of TIPs
ss_VL_without_TIPs_acute <- get_ss_viral_load(4e-6, c(Tt=10000, I=0, I_T=0, I_D=0, V=5, V_M=0, V_T=0, V_H=0))
ss_acute <- newton(silent=TRUE)$state

ss_VL_without_TIPs_chronic <- get_ss_viral_load(1e-5, c(Tt=10000, I=0, I_T=0, I_D=0, V=5, V_M=0, V_T=0, V_H=0))
ss_chronic <- newton(silent=TRUE)$state

#Function to generate each panel
simulate_panel <- function(beta, beta_prime, s_baseline, baseline_VL, label_text) {
  ss_vals <- outer(n_prime_vals, f_vals, Vectorize(function(n_p, f) {
    p <<- c(lambda=2000, d=0.2, b=beta, delta=1, n=1000, c=20, f=f, 
            delta_prime=1, n_prime=n_p, b_prime=beta_prime)
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
    annotate("text", x = 1200, y = 0.46, label = label_text, hjust=0, vjust=0, size=9) +
    theme_minimal() +
    guides(fill="none")
}

panel_settings <- list(
  list(beta=4e-6, beta_prime=4e-6, s=ss_acute, vl=ss_VL_without_TIPs_acute, label=expression(paste(beta," = ",beta,"' = ",4%*%10^-6))),
  list(beta=1e-5, beta_prime=4e-6, s=ss_chronic, vl=ss_VL_without_TIPs_chronic, label=expression(paste(beta > beta,"'"))),
  list(beta=4e-6, beta_prime=1e-5, s=ss_acute, vl=ss_VL_without_TIPs_acute, label=expression(paste(beta < beta,"'"))),
  list(beta=1e-5, beta_prime=1e-5, s=ss_chronic, vl=ss_VL_without_TIPs_chronic, label=expression(paste(beta," = ",beta,"' = ",10^-5)))
)

plots <- lapply(panel_settings, function(cfg) {
  simulate_panel(cfg$beta, cfg$beta_prime, cfg$s, cfg$vl, cfg$label)
})

#Assemble and save final figure
(fig_all <- (plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]]) + 
    plot_layout(nrow=2) & 
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
    label = c("a", "b", "c", "d"),
    x = c(0.005, 0.51, 0.005, 0.51),
    y = c(0.98, 0.98, 0.50, 0.50),
    hjust = 0, vjust = 1,
    size = 35, fontface = "bold"
  ) +
  draw_label(bquote(paste("Virus particles produced by a doubly infected cell (",italic("n'"),")")),
             x = 0.5, y = 0.01, vjust = 0, size = 30) +
  draw_label(bquote(paste("Fraction of wild-type RNA in a doubly infected cell (",italic(f),")")),
             x = 0.03, y = 0.5, angle = 90, vjust = 1, size = 30)


ggsave("../Figures/fig2.png", fig_all_labeled, device="png", width=5250, height=5000, units="px")


###############################
#Generate Figure 3
###############################

#Functions to generate custom axis labels
add_custom_TC_label <- function(plot) {
  ylabel <- textGrob(
    label = expression(atop("Target Cell Count", scriptstyle("(cells mL"^{-1}*")"))),
    rot = 90,
    gp = gpar(fontsize = 35),
    just = c("center"),
    vp = viewport(y = 0.5, height = 1)
  )
  arranged <- arrangeGrob(plot, left = ylabel)
  return(arranged)
}

add_custom_VL_label <- function(plot) {
  ylabel <- textGrob(
    label = expression(atop("Viral Load", scriptstyle("(virions mL"^{-1}*")"))),
    rot = 90,
    gp = gpar(fontsize = 35),
    just = c("center"),
    vp = viewport(y = 0.55, height = 1)
  )
  arranged <- arrangeGrob(plot, left = ylabel)
  return(arranged)
}

#Simulate the model while varying f and tracking target cell counts and VL
simulate_dynamics <- function(f, init_state) {
  n_prime_vals <- seq(1000, 10000, by=20)
  Tt_count <- I_T_count <- healthy_T_count <- total_T_count <- viral_load <- c()
  
  for (i in n_prime_vals) {
    p <<- c(lambda=2000, d=0.2, b=4*10^-6, delta=1, n=1000, c=20,
           f=f, delta_prime=1, n_prime=i, b_prime=4*10^-6)
    s <<- init_state
    s['V_T'] <<- 10
    s <<- run(1000, timeplot=FALSE)
    ss <<- newton(silent=TRUE)$state
    Tt_count <- c(Tt_count, ss["Tt"])
    I_T_count <- c(I_T_count, ss["I_T"])
    healthy_T_count <- c(healthy_T_count, ss["Tt"]+ss["I_T"])
    total_T_count <- c(total_T_count, ss["Tt"]+ss["I_T"]+ss["I_D"]+ss["I"])
    viral_load <- c(viral_load, ss["V"]+ss["V_M"])
  }
  
  data.frame(n_prime_vals, total_T_count, healthy_T_count, 
             Tt_count, I_T_count, viral_load)
}

ss_VL_without_TIPs_acute <- get_ss_viral_load(4e-6, c(Tt=10000, I=0, I_T=0, I_D=0, V=5, V_M=0, V_T=0, V_H=0))
ss_acute <- newton(silent=TRUE)$state

#Simulate the model for f=0.3 and f=0.1 for a range of n' values
df1 <- simulate_dynamics(0.3, ss_acute)
df2 <- simulate_dynamics(0.1, ss_acute)

#Plotting function for target cell count
xlims <- range(df1$n_prime_vals)  # same for all plots
ylim_cells <- c(0, 10500)
ylim_vl <- c(0, 52500)

plot_cells <- function(df, show_x = TRUE, show_y = TRUE) {
  p <- ggplot(df, aes(x=n_prime_vals)) +
    geom_line(aes(y=total_T_count), color="black", size=0.9) +
    geom_line(aes(y=healthy_T_count), color="orange", size=0.9) +
    geom_line(aes(y=Tt_count), color="blue", size=0.9) +
    geom_line(aes(y=I_T_count), color="green", size=0.9) +
    geom_hline(yintercept = 10000, linetype="dashed") +
    theme_cowplot(12) +
    labs(x=NULL, y=NULL) +
    theme(axis.text=element_text(size=25),
          axis.ticks = element_line(size = 0.7)) +
    coord_cartesian(xlim = xlims, ylim = ylim_cells, expand = FALSE) +
    theme(aspect.ratio = 1)
  
  if(!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if(!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  return(p)
}

#Plotting function for viral load
plot_vl <- function(df, show_x = TRUE, show_y = TRUE) {
  p <- ggplot(df, aes(x=n_prime_vals, y=viral_load)) +
    geom_line(size=0.9) +
    theme_cowplot(12) +
    labs(x=NULL, y=NULL) +
    theme(axis.text=element_text(size=25),
          axis.ticks = element_line(size = 0.7)) +
    coord_cartesian(xlim = xlims, ylim = ylim_vl, expand = FALSE) +
    theme(aspect.ratio = 1)
  
  if(!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if(!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  return(p)
}

plot1_cells <- plot_cells(df1, show_x = TRUE, show_y = TRUE)
plot1_cells <- wrap_elements(add_custom_TC_label(plot1_cells))
plot2_cells <- plot_cells(df2, show_x = TRUE, show_y = FALSE)
plot2_cells <- wrap_elements(plot2_cells)

plot1_vl <- plot_vl(df1, show_x = TRUE, show_y = TRUE)
plot1_vl <- wrap_elements(add_custom_VL_label(plot1_vl))
plot2_vl <- plot_vl(df2, show_x = TRUE, show_y = FALSE)
plot2_vl <- wrap_elements(plot2_vl)

#Assemble plot with axis labels
main_plot <- plot1_cells + plot2_cells + plot1_vl + plot2_vl +
  plot_layout(widths = c(1, 1), heights = c(1, 1), ncol=2) &
  theme(
    plot.margin = margin(20, 20, 20, 20),
  )

shared_x_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = as.expression(bquote("Virus particles produced by a doubly infected cell (" * italic(n)*"')")),
           size = 12) +
  theme_void()

final_plot <- main_plot / shared_x_plot + plot_layout(heights = c(1, 0.035))

final_plot_labeled <- ggdraw(final_plot) +
  draw_plot_label(
    label = c("bold(a[1])", "bold(b[1])", "bold(a[2])", "bold(b[2])"), 
    x = c(0.04, 0.53, 0.04, 0.53),
    y = c(0.99, 1, 0.51, 0.52),
    size = 35, family = "Arial", parse=TRUE
  )

ggsave("../Figures/fig3.png", final_plot_labeled, device="png", width=6000, height=5000, units="px", bg="white")


###############################
#Generate Figure 4
###############################

#Compute R0T
r_0 <- function(params) {
  with(as.list(params), {
    return(((((1-f)^2)*n_prime*b_prime)/(n*b))*((1-((c*d)/(n*b*lambda)))))
  })
}

#Scientific label formatting
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#Compute R0T over a range of beta values
compute_r0_over_b <- function(b_vals, b_prime_val) {
  sapply(b_vals, function(bi) {
    params <- c(lambda=2000, d=0.2, b=bi, delta=1, n=1000, c=20, f=0.3, delta_prime=1, n_prime=5000, b_prime=b_prime_val)
    r_0(params)
  })
}

#Run model while varying n' or beta'
simulate_dynamics <- function(var_vals, var_name, f_val, init_state) {
  ss_vals <- numeric(length(var_vals))
  r_0_vals <- numeric(length(var_vals))
  
  for (i in seq_along(var_vals)) {
    p <<- c(lambda=2000, d=0.2, b=4e-6, delta=1, n=1000, c=20, f=f_val, delta_prime=1, n_prime=5000, b_prime=4e-6)
    p[var_name] <<- var_vals[i]
    s <<- init_state
    s['V_T'] <<- 5
    s <<- run(1000, timeplot=FALSE)
    s <<- newton(silent=TRUE)$state
    
    ss_vals[i] <- s["V"] + s["V_M"]
    r_0_vals[i] <- r_0(p)
  }
  
  data.frame(ss_vals = ss_vals, r_0_vals = r_0_vals)
}

b_vals <- seq(1e-6, 1e-4, by=1e-8)
r_0_1_vals <- compute_r0_over_b(b_vals, 4e-6)
r_0_2_vals <- compute_r0_over_b(b_vals, 7e-6)
r_0_3_vals <- compute_r0_over_b(b_vals, 1e-5)

df <- data.frame(b_vals, r_0_1_vals, r_0_2_vals, r_0_3_vals)

ss_VL_without_TIPs_acute <- get_ss_viral_load(4e-6, c(Tt=10000, I=0, I_T=0, I_D=0, V=5, V_M=0, V_T=0, V_H=0))
ss_acute <- newton(silent=TRUE)$state

#Vary n' and beta' for f = 0.3
num_points_to_sample <- 5000

n_prime_vals <- seq(15000, 500, length.out=num_points_to_sample)
b_prime_vals <- seq(1e-6, 1e-4, length.out=num_points_to_sample)

df_2 <- simulate_dynamics(n_prime_vals, "n_prime", f_val=0.3, init_state=ss_acute)
df_3 <- simulate_dynamics(b_prime_vals, "b_prime", f_val=0.3, init_state=ss_acute)

#Vary n' and beta' for f = 0.1
df_4 <- simulate_dynamics(n_prime_vals, "n_prime", f_val=0.1, init_state=ss_acute)
df_5 <- simulate_dynamics(b_prime_vals, "b_prime", f_val=0.1, init_state=ss_acute)

#Function to make the y-axis labels on subplots b and c
add_custom_VL_label <- function(plot) {
  ylabel <- textGrob(
    label = expression(atop("Viral Load", scriptstyle("(virions mL"^{-1}*")"))),
    rot = 90,
    gp = gpar(fontsize = 25),
    just = c("center"),
    vp = viewport(y = 0.6, height = 1)
  )
  arranged <- arrangeGrob(plot, left = ylabel)
  return(arranged)
}

plot1 <- ggplot(df, aes(x=b_vals)) + 
  geom_line(aes(y=r_0_1_vals), color="green", size=0.8) + 
  geom_line(aes(y=r_0_2_vals), color="blue", size=0.8) + 
  geom_line(aes(y=r_0_3_vals), color="orange", size=0.8) + 
  coord_cartesian(xlim=c(2e-6, 2e-5), ylim=c(1, 3.2), expand=TRUE) + 
  scale_x_continuous(breaks = c(2e-6, 8e-6, 1.4e-5, 2e-5), labels=scientific_10) + 
  labs(x=expression("Wild-type infection rate (" * beta * ")"), y=expression(R[0]^T)) +
  theme_cowplot(12) +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20)) +
  geom_vline(xintercept=4e-6, linetype="dashed")

plot2 <- ggplot(df_2, aes(x=r_0_vals, y=ss_vals)) + 
  geom_line(color="red", size=0.8) + 
  geom_line(data=df_3, aes(x=r_0_vals, y=ss_vals), color="blue", linetype="dashed", size=0.8) +
  coord_cartesian(xlim=c(0.8, 3), ylim=c(30000, 60000), expand=TRUE) +
  labs(x = expression(R[0]^T), y = NULL) +
  theme_cowplot(12) +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20))

plot3 <- ggplot(df_4, aes(x=r_0_vals, y=ss_vals)) + 
  geom_line(color="red", size=0.8) + 
  geom_line(data=df_5, aes(x=r_0_vals, y=ss_vals), color="blue", linetype="dashed", size=0.8) +
  coord_cartesian(xlim=c(0.8, 3), ylim=c(10000, 60000), expand=TRUE) +
  labs(x = expression(R[0]^T), y = NULL) +
  theme_cowplot(12) +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20))

plot1 <- wrap_elements(plot1)
plot2 <- wrap_elements(add_custom_VL_label(plot2))
plot3 <- wrap_elements(add_custom_VL_label(plot3))

#Assemble final plot
final_plot <- plot1 | (plot2 / plot3) + 
  plot_layout(widths = c(1, 1), heights = c(1, 1)) &
  theme(
    plot.tag = element_text(size=35, face="bold", family="Arial"),
    plot.tag.position = c(-0.01, 1),
    plot.margin = margin(20, 20, 20, 20),
    axis.text = element_text(size=22),
    axis.title = element_text(size=25),
    axis.ticks = element_line(linewidth=0.7)
  )

final_plot_labeled <- ggdraw(final_plot) +
  draw_plot_label(
  label = c("a", "b", "c"), 
  x = c(0.01, 0.52, 0.52),  # x positions of each panel tag
  y = c(0.99, 0.99, 0.49),  # y positions
  size = 35, fontface = "bold", family = "Arial"
)

ggsave("../Figures/fig4.png", plot=final_plot_labeled, device="png", width=5000, height=3000, units="px")