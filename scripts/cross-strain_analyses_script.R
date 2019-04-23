# Project: Microalgae Cross-strain Growth Medium Reuse Experiment
# Author: Sarah Loftus, sarah.e.loftus@gmail.com

# This code organizes, analyses, and visualizes data from mircoalgae cross-strain reused medium experiments. 

# Input data and metadata are also available on Figshare: https://doi.org/10.6084/m9.figshare.7831913.v1

#######################################################################################

### Install packages ###

# Data management
library(tidyverse) 
library(lubridate)
library(broom)

# Plotting
library(gtable)
library(grid)
library(ggpubr)

### Read Data (see metadata for CSV file column descriptions)
growth_df <- read.csv("input_data/Data1_Growth.csv")        # Algae growth and experimental data for each culture replicate
daily_df <- read.csv("input_data/Data2_Daily.csv")          # Daily sampling data 

#### Parse dates and times, and calculate Days elapsed ####

daily_df %>%
  mutate(Date = mdy(Date),   
         Tstart = hm(Tstart),
         Tend = hm(Tend),
         Tmid = Tstart + as.period(0.5*as.duration(Tend - Tstart))) %>%   # Tmid = midpoint between start and end of sampling time
  do(mutate(., DaysElapsed = ((Date + Tmid) - (Date[Round == 0 & Day == 0] + Tmid[Round == 0 & Day == 0]) )/ 86400) ) %>%   # DaysElapsed =  days elapsed since start of each experiment
  select(Round, Day, Chl_medium, OD750_medium, DaysElapsed) -> daily_df2 
  # Note: warning messages do not pose problem to further analysis

#### Define and calculate additional variables ####

growth_df %>%    
  inner_join(daily_df2, by = c("Round", "Day") ) %>% # Add daily measurement data to growth data frame
  mutate(biomass_OD = OD750 -  OD750_medium,         # biomass_OD = Blank-corrected OD750 of the culture, in arbitrary units
         biomass_chl = Chl - Chl_medium              # biomass_chl = Blank-corrected chlorophyll-a concentration of the culture, in relative fluorescence units
  ) -> growth_df2    

# # Check regression of OD vs algae cell concentrations to compare with other experiments
# OD_vs_cell_Navicula <- lm(biomass_OD ~ AlgaeConc, data = growth_df2[growth_df2$Algae == "Navicula", ])  
# # R^2 = 0.8081, biomass OD = 0.30008*AlgaeConc + 0.01173
# 
# OD_vs_cell_C323 <- lm(biomass_OD ~ AlgaeConc, data = growth_df2[growth_df2$Algae == 'C323', ])  
# # R^2 = 0.5192, biomass OD = 0.03825*AlgaeConc + 0.01979
# 
# OD_vs_cell_D046 <- lm(biomass_OD ~ AlgaeConc, data = growth_df2[growth_df2$Algae == 'D046', ])  
# # R^2 = 0.7837, biomass OD = 0.04281*AlgaeConc + 0.016368
# 
# # Plots
# plot(biomass_OD ~ AlgaeConc, data = growth_df2[growth_df2$Algae == "Navicula", ], 
#      xlab = "Algae Concentration (million cells/mL)", ylab = "OD 750", main = "Navicula sp.", ylim = c(0, 0.13), xlim = c(0,1))
# abline(lm(growth_df2$biomass_OD[growth_df2$Algae == "Navicula"] ~ growth_df2$AlgaeConc[growth_df2$Algae == "Navicula"]))
# 
# plot(biomass_OD ~ AlgaeConc, data = growth_df2[growth_df2$Algae == "D046", ], 
#      xlab = "Algae Concentration (million cells/mL)", ylab = "OD 750", main = "Chlorella sp. D046", xlim = c(0,8))
# abline(lm(growth_df2$biomass_OD[growth_df2$Algae == "D046"] ~ growth_df2$AlgaeConc[growth_df2$Algae == "D046"]))
# 
# plot(biomass_OD ~ AlgaeConc, data = growth_df2[growth_df2$Algae == "C323", ], 
#      xlab = "Algae Concentration (million cells/mL)", ylab = "OD 750", main = "Staurosira sp. C323", ylim = c(0, 0.2))
# abline(lm(growth_df2$biomass_OD[growth_df2$Algae == "C323"] ~ growth_df2$AlgaeConc[growth_df2$Algae == "C323"]))


# Calculate daily averages and standard deviations of all variables in growth data frame. For use in plotting. ###

growth_df2 %>%
  select(-Replicate, -Chl, -OD750, -Chl_medium, -OD750_medium) %>% 
  group_by(Algae, Round, Treatment, Day) %>% 
  summarize_all(funs(mean, sd), na.rm = TRUE) %>%
  replace(., is.na(.), NA)-> growth_df_avgs  

#### Calculate variables for algae growth statistics ####

# Define time periods where cultures are in exponential phase for calculating specific growth rates (*** this is based on OD growth phase)
mu_period_02 <- 0:2   # Days 0-2; for: Navicula Round 3 all treatments/replicates
mu_period_03 <- 0:3   # Days 0-3; for: C323 Round 1 Fresh C, C323 Round 3 Fresh all replciates
mu_period_04 <- 0:4   # Days 0-4; for: C323 Round 1 Recycled B, D046 Round 3 all treatments/replicates
mu_period_13 <- 1:3   # Days 1-3; for: C323 Round 1 Fresh A
mu_period_14 <- 1:4   # Days 1-4; for: C323 Round 1 Fresh B & Recycled A & Recycled C, C323 Round 2 Fresh all replicates
mu_period_25 <- 2:5   # Days 2-5; for: C323 Round 0 all treatments/replicates
mu_period_35 <- 3:5   # Days 3-5; for: C323 Round 2 Recycled all replicates
mu_period_38 <- 3:8   # Days 3-8; for: C323 Round 3 Recycled all replciates

# Specific growth rate 
growth_df2 %>%
  filter( (Algae == "Navicula" & Day %in% mu_period_02) |     # Based on the stated definitions above, select only the days used to calculate specific growth rate
          ( (Algae == "C323" & ( (Round == 1 & Treatment == "F" & Replicate == "C") | (Round == 3 & Treatment == "F") ) ) & Day %in% mu_period_03) |
          ( ((Algae == "C323" & Round == 1 & Treatment == "R" & Replicate == "E") | (Algae == "D046" & Round == 3)) & Day %in% mu_period_04) |
          (Algae == "C323" & Round == 1 & Treatment == "F" & Replicate == "A" & Day %in% mu_period_13) |
          ( (Algae == "C323" & ( ( Round == 1 & ( (Treatment == "F" & Replicate == "B") | (Treatment == "R" & (Replicate == "D" | Replicate == "F")))) | (Round == 2 & Treatment == "F"))) & Day %in% mu_period_14) | 
          (Algae == "C323" & Round == 0 & Day %in% mu_period_25) |
          (Algae == "C323" & Round == 2 & Treatment == "R" & Day %in% mu_period_35) |
          (Algae == "C323" & Round == 3 & Treatment == "R" & Day %in% mu_period_38) ) %>%   
  group_by(Algae, Round, Treatment, Replicate) %>%
  do(tidy(lm(log(biomass_OD) ~ DaysElapsed, data = .))) %>%   # linear regression of natural log algae cell concentration versus time
  filter(term == "DaysElapsed") %>%                          # where "DaysElapsed" is the row term in the linear regression output to indicate the slope, i.e., the specific growth rate
  rename(mu = estimate) %>%                                  # rename the slope of the linear regression ("estimate") as mu, mu = the specific growth rate (units = 1/day)
  select(Algae, Round, Treatment, Replicate, mu) -> mu_data

mu_data %>%
  group_by(Algae, Round, Treatment) %>% 
  select(-Replicate) %>% 
  summarize_all(funs(mean, sd), na.rm = TRUE) -> mu_avgs  

# Maximum OD over the growth period for each replicate, then take avergae and st dev

growth_df2 %>% 
  select(Round, Day, Algae, Treatment, Replicate, biomass_OD) %>% 
  group_by(Round, Algae, Treatment, Replicate) %>% 
  summarize(max_OD = max(biomass_OD)) -> max_OD_data 

max_OD_data %>% 
  group_by(Algae, Round, Treatment) %>% 
  select(-Replicate) %>% 
  summarize_all(funs(mean, sd), na.rm = TRUE) -> max_OD_avgs 

# Change in DOC in Round 3

growth_df2 %>% 
  filter(Round == 3) %>% 
  select(Day, Algae, Treatment, Replicate, DOC, biomass_OD) %>% 
  group_by(Algae, Treatment, Replicate) %>% 
  summarize(DOC_change = DOC[Day == 8] - DOC[Day == 0], # Change in DOC over the growth period
            DOC_per_OD = DOC_change/(biomass_OD[Day == 8] - biomass_OD[Day == 0])) -> DOC_change # Change in DOC per change in OD over the growth period
  
DOC_change %>%   
  select(-Replicate) %>% 
  summarize_all(funs(mean, sd), na.rm = TRUE) -> avg_DOC_change



#### Statistics ####

# Test for equal variances of specific growth rate (mu) and maximum OD, then perform t-tests

# mu - variance test
C323_mu_var <- var.test(mu_data$mu[mu_data$Algae == 'C323' & mu_data$Round == 3 & mu_data$Treatment == "F"], 
                       mu_data$mu[mu_data$Algae == 'C323' & mu_data$Round == 3 & mu_data$Treatment == "R"],
                       alternative = "two.sided")

D046_mu_var <- var.test(mu_data$mu[mu_data$Algae == 'D046' & mu_data$Treatment == "F"], 
                       mu_data$mu[mu_data$Algae == 'D046' & mu_data$Treatment == "R"],
                       alternative = "two.sided")

Navi_mu_test <- var.test(mu_data$mu[mu_data$Algae == 'Navicula' & mu_data$Treatment == "F"], 
                       mu_data$mu[mu_data$Algae == 'Navicula' & mu_data$Treatment == "R"],
                       alternative = "two.sided")

# mu - T test
C323_mu_test <- t.test(mu_data$mu[mu_data$Algae == 'C323' & mu_data$Round == 3 & mu_data$Treatment == "F"], 
                       mu_data$mu[mu_data$Algae == 'C323' & mu_data$Round == 3 & mu_data$Treatment == "R"],
                       alternative = "two.sided", paired = FALSE, var.equal = TRUE)

D046_mu_test <- t.test(mu_data$mu[mu_data$Algae == 'D046' & mu_data$Treatment == "F"], 
                        mu_data$mu[mu_data$Algae == 'D046' & mu_data$Treatment == "R"],
                        alternative = "two.sided", paired = FALSE, var.equal = TRUE)

Navi_mu_test <- t.test(mu_data$mu[mu_data$Algae == 'Navicula' & mu_data$Treatment == "F"], 
                        mu_data$mu[mu_data$Algae == 'Navicula' & mu_data$Treatment == "R"],
                        alternative = "two.sided", paired = FALSE, var.equal = TRUE)

# maximum OD - var test
C323_max_OD_var <- var.test(max_OD_data$max_OD[max_OD_data$Algae == 'C323' & max_OD_data$Round == 3 & max_OD_data$Treatment == "F"], 
                           max_OD_data$max_OD[max_OD_data$Algae == 'C323' & max_OD_data$Round == 3 & max_OD_data$Treatment == "R"],
                           alternative = "two.sided")

D046_max_OD_var <- var.test(max_OD_data$max_OD[max_OD_data$Algae == 'D046' & max_OD_data$Treatment == "F"], 
                           max_OD_data$max_OD[max_OD_data$Algae == 'D046' & max_OD_data$Treatment == "R"],
                           alternative = "two.sided")

Navi_max_OD_var <- var.test(max_OD_data$max_OD[max_OD_data$Algae == 'Navicula' & max_OD_data$Treatment == "F"], 
                           max_OD_data$max_OD[max_OD_data$Algae == 'Navicula' & max_OD_data$Treatment == "R"],
                           alternative = "two.sided")

# maximum OD - T test
C323_max_OD_test <- t.test(max_OD_data$max_OD[max_OD_data$Algae == 'C323' & max_OD_data$Round == 3 & max_OD_data$Treatment == "F"], 
                       max_OD_data$max_OD[max_OD_data$Algae == 'C323' & max_OD_data$Round == 3 & max_OD_data$Treatment == "R"],
                       alternative = "two.sided", paired = FALSE, var.equal = TRUE)

D046_max_OD_test <- t.test(max_OD_data$max_OD[max_OD_data$Algae == 'D046' & max_OD_data$Treatment == "F"], 
                       max_OD_data$max_OD[max_OD_data$Algae == 'D046' & max_OD_data$Treatment == "R"],
                       alternative = "two.sided", paired = FALSE, var.equal = FALSE)

Navi_max_OD_test <- t.test(max_OD_data$max_OD[max_OD_data$Algae == 'Navicula' & max_OD_data$Treatment == "F"], 
                       max_OD_data$max_OD[max_OD_data$Algae == 'Navicula' & max_OD_data$Treatment == "R"],
                       alternative = "two.sided", paired = FALSE, var.equal = TRUE)

#### Figures ####

#### Figure 1: Daily OD ####

  # C323 Round 0-2 ####  
  C323_OD_plot_R02 <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & growth_df_avgs$Round != 3, ], 
                             aes(x = as.numeric(DaysElapsed_mean), y = biomass_OD_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=biomass_OD_mean - biomass_OD_sd, ymax=biomass_OD_mean + biomass_OD_sd), width= 0.4) +
    labs(x = "Days") +
    labs(y = expression("OD"["750"])) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 0.15, 0.05), expand = expand_scale(mult = c(0.05,0.1))) +  
    scale_x_continuous(limits = c(-0.1,20), breaks = seq(0, 21, 2)) +
    annotate("text", x = c(3, 9, 16), y = 0.16, label = c("Fresh medium", "1st reuse", "2nd reuse"), color = "gray40", size = 5) +
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 9),
           legend.position = c(0.1,0.72),
           axis.title.y = element_text(margin = margin(r =10), size = 14),  
           axis.title.x = element_text(margin = margin(r =35), size = 14),   
           panel.background = element_rect(fill = "white"),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           axis.ticks = element_blank(),
           axis.text = element_text(size = 12)) 
  
  # Save plot as pdf file
  pdf("Figures/Fig1_OD_Round02.pdf", 
      width = 7.5, height = 3)
  grid.draw(C323_OD_plot_R02)
  dev.off()
  
  # D046 Round 3 ####    
  D046_OD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046", ], 
                         aes(x = Day, y = biomass_OD_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=biomass_OD_mean - biomass_OD_sd, ymax=biomass_OD_mean + biomass_OD_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 0.20, 0.05), expand = expand_scale(mult = c(0.05,0.1))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs( y = " ") +
    annotate("text", x = 1.9, y = 0.22, label = expression(paste(bold("A"), italic("  Chlorella"), " sp. D046")), size = 4) + 
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 9),
           axis.title.y = element_text(margin = margin(r = 15), size = 14),
           axis.title.x = element_blank(),
           legend.position = c(0.17,0.6),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 12))   
  
  # Navicula Round 3 #### 
  Navicula_OD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula", ], 
                             aes(x = Day, y = biomass_OD_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=biomass_OD_mean - biomass_OD_sd, ymax=biomass_OD_mean + biomass_OD_sd), width= 0.2) +
    labs(y = expression("OD"["750"])) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 0.1, 0.03), expand = expand_scale(mult = c(0.05, 0.05))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    annotate("text", x = 1.7, y = 0.095, label = expression(paste( bold("B"), italic("  Navicula"), " sp. SFP")), size = 4) + 
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 9),
           legend.position = c(0.17,0.6),
           axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
           axis.title.x = element_blank(),   
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA),            
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 12)) 
  
  # C323 Round 3 ####  
  C323_OD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & growth_df_avgs$Round == 3, ], 
                         aes(x = Day, y = biomass_OD_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=biomass_OD_mean - biomass_OD_sd, ymax=biomass_OD_mean + biomass_OD_sd), width= 0.2) +
    labs(x = "Days") +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 0.1, 0.03), expand = expand_scale(mult = c(0.05,0.15))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    annotate("text", x = 2, y = 0.085, label = expression(paste(bold("C"), italic("  Staurosira"), " sp. C323")), size = 4) + 
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 9),
           legend.position = c(0.17,0.6),
           axis.title.y = element_blank(),  # moves axis title away from axis label
           axis.title.x = element_text(margin = margin(r = 35), size = 14),   
           panel.background = element_rect(fill = "white"),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           axis.ticks = element_blank(),
           axis.text = element_text(size = 12))   
  
  # Combine Round 3 Daily OD plots in a grid & Save as pdf ####
  
  # Convert ggplots to grobs
  D046_OD_grob <- ggplotGrob(D046_OD_plot)
  C323_OD_grob <- ggplotGrob(C323_OD_plot)
  Navi_OD_grob <- ggplotGrob(Navicula_OD_plot)
  
  # Combine in grid
  OD_gridplot <- rbind(D046_OD_grob, Navi_OD_grob, C323_OD_grob, size = "first")
  
  grid.newpage()
  grid.draw(OD_gridplot)
  
  # Save plot as pdf file
  pdf("Figures/Fig2_ODGrid.pdf", 
      width = 4, height = 6)
  grid.draw(OD_gridplot)
  dev.off()
  
#### Figure S1: Specific growth rate and maximum OD in phase 2 ####
  
  # reorder the x-axis to D046, Navi, C323
  mu_avgs$Algae <- factor(mu_avgs$Algae, levels = c("D046", "Navicula", "C323"))
  max_OD_avgs$Algae <- factor(max_OD_avgs$Algae, levels = c("D046", "Navicula", "C323"))
  
  specific_growth_plot <- ggplot(data = mu_avgs[mu_avgs$Round == 3, ], 
                         aes(x = Algae, y = mean, fill = Treatment))  +
    geom_col(position = position_dodge(0.9), color = "gray15") +
    geom_errorbar(aes(ymin=mean - sd, ymax= mean + sd), width= 0.4, position = position_dodge(0.9)) +
    labs(y = expression(paste("Specific growth rate ", ("d"^-1)))) +
    scale_fill_manual(labels = c("Fresh", "Reused"), values = c('gray53','gray33')) +
    scale_x_discrete(labels=c("D046" = expression(paste(italic("Chlorella"), " sp.")), "Navicula" = expression(paste(italic("Navicula"), " sp.")), "C323" = expression(paste(italic("Staurosira"), " sp.")))) +
    annotate("text", x = 0.6, y = 1.2, label = expression(paste(bold("A"))), size = 5) + 
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05)), breaks = seq(0,1.2, by = 0.2)) +
    theme( legend.key = element_rect(fill = NA),  
           legend.title = element_blank(),
           legend.text = element_text(size = 11, margin = margin(r = 8)),
           legend.position = "top",
           legend.key.size = unit(5, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 10),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           panel.background = element_rect(fill = NA),
           axis.line = element_line(color = "black", size = 0.4),
           axis.ticks.x = element_blank(),
           axis.title.y = element_text(margin = margin(r = 6),size = 12),
           axis.text.y = element_text(size = 12),
           axis.title.x = element_blank(),
           axis.text.x = element_blank()) 
  
  max_OD_plot <- ggplot(data = max_OD_avgs[max_OD_avgs$Round == 3, ], 
                                 aes(x = Algae, y = mean, fill = Treatment))  +
    geom_col(position = position_dodge(0.9), color = "gray15") +
    geom_errorbar(aes(ymin=mean - sd, ymax= mean + sd), width= 0.4, position = position_dodge(0.9)) +
    labs(y = expression("Maximum OD"["750"])) +
    scale_fill_manual(labels = c("Fresh", "Reused"), values = c('gray53','gray33')) +
    scale_x_discrete(labels=c("D046" = expression(paste(italic("Chlorella"), " sp.")), "Navicula" = expression(paste(italic("Navicula"), " sp.")), "C323" = expression(paste(italic("Staurosira"), " sp.")))) +
    annotate("text", x = 0.6, y = 0.22, label = expression(paste(bold("B"))), size = 5) + 
    scale_y_continuous(expand = expand_scale(mult = c(0.02, 0.05)), breaks = seq(0,0.2, by = 0.05)) +
    theme( legend.position = "none",                              
           panel.background = element_rect(fill = NA),
           axis.line = element_line(color = "black", size = 0.4),
           axis.ticks.x = element_blank(),
           axis.title = element_text(margin = margin(r = 6),size = 12),
           axis.text = element_text(size = 11),
           axis.title.x = element_blank(),
           axis.title.y = element_text(margin = margin(r = 6),size = 12))   
  
  # Combine plots
  mu_grob <- ggplotGrob(specific_growth_plot)
  max_OD_grob <- ggplotGrob(max_OD_plot)
  
  mu_maxOD_grid <- rbind(mu_grob, max_OD_grob, size = "first")
  
  # Save plot as pdf file
  pdf("Figures/FigS1_mu_maxOD_Grid.pdf", 
      width = 4.5, height = 6.5)
  grid.draw(mu_maxOD_grid)
  dev.off()
  
#### Nutrients versus Time ####    
#### PO4
  # D046 Round 3 ####    
  D046_PO4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !(is.na(growth_df_avgs$PO4_mean)), ], 
                         aes(x = Day, y = PO4_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=PO4_mean - PO4_sd, ymax=PO4_mean + PO4_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 175, 50), expand = expand_scale(mult = c(0.05,0.05))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs( y = "Phosphate (µM)", title = expression(paste(italic("  Chlorella"), " sp. D046"))) +
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 9),
           axis.title.y = element_text(margin = margin(r = 15), size = 11),
           axis.title.x = element_blank(),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 11)) 

  # Navicula Round 3 ####    
  Navi_PO4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !(is.na(growth_df_avgs$PO4_mean)), ], 
                          aes(x = Day, y = PO4_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=PO4_mean - PO4_sd, ymax=PO4_mean + PO4_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 175, 50), expand = expand_scale(mult = c(0.05,0.05))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs(title = expression(paste(italic("  Navicula"), " sp."))) +
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 9),
           axis.title.y = element_blank(),
           axis.title.x = element_blank(),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 11))   
  
  # C323 Round 3 ####    
  C323_PO4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !(is.na(growth_df_avgs$PO4_mean)) & growth_df_avgs$Round == 3, ], 
                          aes(x = Day, y = PO4_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=PO4_mean - PO4_sd, ymax=PO4_mean + PO4_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 175, 50), expand = expand_scale(mult = c(0.05,0.01)), limits = c(0, 200)) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs(title = expression(paste(italic("Staurosira"), " sp. C323"))) +
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 9),
           axis.title.y = element_blank(),
           axis.title.x = element_blank(),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 11)) 

# NH4
  # D046 Round 3 ####    
  D046_NH4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !(is.na(growth_df_avgs$NH4_mean)), ], 
                          aes(x = Day, y = NH4_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=NH4_mean - NH4_sd, ymax=NH4_mean + NH4_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 450, 100), expand = expand_scale(mult = c(0.05,0.05))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs( y = "Ammonium (µM)") +
    theme( legend.position="none",
           axis.title.y = element_text(margin = margin(r = 15), size = 11),
           axis.title.x = element_blank(),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 11)) 
  
  # Navicula Round 3 ####    
  Navi_NH4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !(is.na(growth_df_avgs$NH4_mean)), ], 
                          aes(x = Day, y = NH4_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=NH4_mean - NH4_sd, ymax=NH4_mean + NH4_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 450, 100), expand = expand_scale(mult = c(0.05,0.05))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    theme( legend.position="none",
           axis.title.y = element_blank(),
           axis.title.x = element_blank(),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 11))   
  
  # C323 Round 3 ####    
  C323_NH4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !(is.na(growth_df_avgs$NH4_mean)) & growth_df_avgs$Round == 3, ], 
                          aes(x = Day, y = NH4_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=NH4_mean - NH4_sd, ymax=NH4_mean + NH4_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 400, 100), expand = expand_scale(mult = c(0.05,0.01)), limits = c(0, 450)) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    theme( legend.position="none",
           axis.title.y = element_blank(),
           axis.title.x = element_blank(),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 11)) 
 
# Silicate
  # D046 Round 3 ####    
  D046_Si_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !(is.na(growth_df_avgs$Si_mean)), ], 
                          aes(x = Day, y = Si_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=Si_mean - Si_sd, ymax=Si_mean + Si_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 1400, 200), expand = expand_scale(mult = c(0.05,0.05)), limits = c(0, 900)) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs( y = "Silicate (µM)", x = "Days") +
    theme( legend.position="none",
           axis.title.y = element_text(margin = margin(r = 15), size = 11),
           axis.title.x = element_text(margin = margin(r = 15), size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 11),
           axis.text.y = element_text(size = 11)) 
  
  # Navicula Round 3 ####    
  Navi_Si_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !(is.na(growth_df_avgs$Si_mean)), ], 
                          aes(x = Day, y = Si_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=Si_mean - Si_sd, ymax=Si_mean + Si_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 1400, 200), expand = expand_scale(mult = c(0.05,0.05))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs(x = "Days") +
    theme( legend.position="none",
           axis.title.y = element_blank(),
           axis.title.x = element_text(margin = margin(r = 15), size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 11),
           axis.text.y = element_text(size = 11))   
  
  # C323 Round 3 ####    
  C323_Si_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !(is.na(growth_df_avgs$Si_mean)) & growth_df_avgs$Round == 3, ], 
                          aes(x = Day, y = Si_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=Si_mean - Si_sd, ymax=Si_mean + Si_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 1400, 200), expand = expand_scale(mult = c(0.05,0.01)), limits = c(0, 900)) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 1)) +
    labs(x = "Days") +
    theme( legend.position="none",
           axis.title.y = element_blank(),
           axis.title.x = element_text(margin = margin(r = 15), size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 11),
           axis.text.y = element_text(size = 11)) 

  # Combine 9 plots in a 3x3 grid
  # Convert plots to grobs
  D046_PO4_grob <- ggplotGrob(D046_PO4_plot)
  Navi_PO4_grob <- ggplotGrob(Navi_PO4_plot)
  C323_PO4_grob <- ggplotGrob(C323_PO4_plot)
  D046_NH4_grob <- ggplotGrob(D046_NH4_plot)
  Navi_NH4_grob <- ggplotGrob(Navi_NH4_plot)
  C323_NH4_grob <- ggplotGrob(C323_NH4_plot)
  D046_Si_grob <- ggplotGrob(D046_Si_plot)
  Navi_Si_grob <- ggplotGrob(Navi_Si_plot)
  C323_Si_grob <- ggplotGrob(C323_Si_plot)  
  
  # Make plot columns
  nuts_col1 <- rbind(D046_PO4_grob, D046_NH4_grob, D046_Si_grob, size = "first")
  nuts_col2 <- rbind(Navi_PO4_grob, Navi_NH4_grob, Navi_Si_grob, size = "first")
  nuts_col3 <- rbind(C323_PO4_grob, C323_NH4_grob, C323_Si_grob, size = "first")
  # Set widths
  nuts_col1$widths <- unit.pmax(D046_PO4_grob$widths, D046_NH4_grob$widths, D046_Si_grob$widths)
  nuts_col2$widths <- unit.pmax(Navi_PO4_grob$widths, Navi_NH4_grob$widths, Navi_Si_grob$widths)
  nuts_col3$widths <- unit.pmax(C323_PO4_grob$widths, C323_NH4_grob$widths, C323_Si_grob$widths)
  # Combine plot columns
  nutrients_plot <- cbind(nuts_col1, nuts_col2, nuts_col3, size = "first")
  
  grid.newpage()
  grid.draw(nutrients_plot)
  
  # Save plot as pdf file
  pdf("Figures/Nutrients.pdf", 
      width = 7, height = 5)
  grid.draw(nutrients_plot)
  dev.off()
  
#### intial nutrients ####
  
  growth_df_avgs$Algae <- factor(growth_df_avgs$Algae, levels = c("D046", "Navicula", "C323"))

  initial_NH4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Round == 3 & growth_df_avgs$Day == 0, ], 
                                 aes(x = Algae, y = NH4_mean, fill = Treatment))  +
    geom_col(position = position_dodge(0.9), color = "gray15") +
    geom_errorbar(aes(ymin=NH4_mean - NH4_sd, ymax= NH4_mean + NH4_sd), width= 0.4, position = position_dodge(0.9)) +
    labs(y = "Ammonium (µM)") +
    scale_fill_manual(labels = c("Fresh", "Reused"), values = c('gray53','gray33')) +
    scale_x_discrete(labels=c("D046" = expression(paste(italic("Chlorella"), " sp.")), "Navicula" = expression(paste(italic("Navicula"), " sp.")), "C323" = expression(paste(italic("Staurosira"), " sp.")))) +
    annotate("text", x = 0.6, y = 450, label = expression(paste(bold("A"))), size = 5) + 
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05)), breaks = seq(0,450, by = 100)) +
    theme( legend.key = element_rect(fill = NA),  
           legend.title = element_blank(),
           legend.text = element_text(size = 11, margin = margin(r = 8)),
           legend.position = "top",
           legend.key.size = unit(5, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 10),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           panel.background = element_rect(fill = NA),
           axis.line = element_line(color = "black", size = 0.4),
           axis.ticks.x = element_blank(),
           axis.title.y = element_text(margin = margin(r = 6),size = 12),
           axis.text.y = element_text(size = 12),
           axis.title.x = element_blank(),
           axis.text.x = element_blank()) 
  
  initial_PO4_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Round == 3 & growth_df_avgs$Day == 0, ], 
                        aes(x = Algae, y = PO4_mean, fill = Treatment))  +
    geom_col(position = position_dodge(0.9), color = "gray15") +
    geom_errorbar(aes(ymin=PO4_mean - PO4_sd, ymax= PO4_mean + PO4_sd), width= 0.4, position = position_dodge(0.9)) +
    labs(y = "Phosphate (µM)") +
    scale_fill_manual(labels = c("Fresh", "Reused"), values = c('gray53','gray33')) +
    scale_x_discrete(labels=c("D046" = expression(paste(italic("Chlorella"), " sp.")), "Navicula" = expression(paste(italic("Navicula"), " sp.")), "C323" = expression(paste(italic("Staurosira"), " sp.")))) +
    annotate("text", x = 0.6, y = 200, label = expression(paste(bold("B"))), size = 5) + 
    scale_y_continuous(expand = expand_scale(mult = c(0.02, 0.05)), breaks = seq(0,200, by = 50)) +
    theme( legend.position = "none",                              
           panel.background = element_rect(fill = NA),
           axis.line = element_line(color = "black", size = 0.4),
           axis.ticks.x = element_blank(),
           axis.title.y = element_text(margin = margin(r = 6),size = 12),
           axis.text.y = element_text(size = 12),
           axis.title.x = element_blank(),
           axis.text.x = element_blank()) 
  
  initial_Si_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Round == 3 & growth_df_avgs$Day == 0, ], 
                             aes(x = Algae, y = Si_mean, fill = Treatment))  +
    geom_col(position = position_dodge(0.9), color = "gray15") +
    geom_errorbar(aes(ymin=Si_mean - Si_sd, ymax= Si_mean + Si_sd), width= 0.4, position = position_dodge(0.9)) +
    labs(y = "Silicate (µM)") +
    scale_fill_manual(labels = c("Fresh", "Reused"), values = c('gray53','gray33')) +
    scale_x_discrete(labels=c("D046" = expression(paste(italic("Chlorella"), " sp.")), "Navicula" = expression(paste(italic("Navicula"), " sp.")), "C323" = expression(paste(italic("Staurosira"), " sp.")))) +
    annotate("text", x = 0.6, y = 900, label = expression(paste(bold("C"))), size = 5) + 
    scale_y_continuous(expand = expand_scale(mult = c(0.02, 0.05)), breaks = seq(0,1500, by = 200)) +
    theme( legend.position = "none",                              
           panel.background = element_rect(fill = NA),
           axis.line = element_line(color = "black", size = 0.4),
           axis.ticks.x = element_blank(),
           axis.title = element_text(margin = margin(r = 6),size = 12),
           axis.text = element_text(size = 12),
           axis.title.x = element_blank(),
           axis.title.y = element_text(margin = margin(r = 6),size = 12))
  
  # Combine plots
  init_NH4_grob <- ggplotGrob(initial_NH4_plot)
  init_PO4_grob <- ggplotGrob(initial_PO4_plot)
  init_Si_grob <- ggplotGrob(initial_Si_plot)
  
  init_nuts_grid <- rbind(init_NH4_grob, init_PO4_grob, init_Si_grob, size = "first")
  
  # Save plot as pdf file
  pdf("Figures/FigS2_init_nuts.pdf", 
      width = 4, height = 7.5)
  grid.draw(init_nuts_grid)
  dev.off()

  
  #### DOC ####
  
  D046_DOC_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !(is.na(growth_df_avgs$DOC_mean)), ], 
                         aes(x = Day, y = DOC_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=DOC_mean - DOC_sd, ymax=DOC_mean + DOC_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 500, 50), expand = expand_scale(mult = c(0.05,0.1))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 2)) +
    annotate("text", x = 0.2, y = 220, label = expression(paste(bold("A"))), size = 4) + 
    labs( y = "Biologically-derived DOC (µM)", x = "Days", title = expression(paste(italic("  Chlorella"), " sp. D046"))) +
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 10),
           axis.title = element_text(margin = margin(r = 15), size = 11),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 12),
           axis.text.y = element_text(size = 12))   
  
  Navicula_DOC_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !(is.na(growth_df_avgs$DOC_mean)), ], 
                             aes(x = Day, y = DOC_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=DOC_mean - DOC_sd, ymax=DOC_mean + DOC_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 500, 50), expand = expand_scale(mult = c(0.05,0.1))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 2)) +
    annotate("text", x = 0.2, y = 260, label = expression(paste(bold("B"))), size = 4) + 
    labs(x = "Days", title = expression(paste(italic("  Navicula"), " sp. SFP"))) +
    theme( legend.key = element_rect(fill = NA),  
           legend.title = element_blank(),
           legend.text = element_text(size = 10),
           axis.title.x = element_text(margin = margin(r = 15), size = 11),
           axis.title.y = element_blank(),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 12),
           axis.text.y = element_text(size = 12))
  
  C323_DOC_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !(is.na(growth_df_avgs$DOC_mean)) & growth_df_avgs$Round == 3, ], 
                         aes(x = Day, y = DOC_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=DOC_mean - DOC_sd, ymax=DOC_mean + DOC_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 500, 50), expand = expand_scale(mult = c(0.05,0.1))) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 2)) +
    annotate("text", x = 0.2, y = 175, label = expression(paste(bold("C"))), size = 4) + 
    labs(x = "Days", title = expression(paste(italic("Staurosira"), " sp. C323"))) +
    theme( legend.key = element_rect(fill = NA),  
           legend.title = element_blank(),
           legend.text = element_text(size = 10),
           axis.title.x = element_text(margin = margin(r = 15), size = 11),
           axis.title.y = element_blank(),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 12),
           axis.text.y = element_text(size = 12))
  
  # Combine Round 3 Daily OD plots in a grid & Save as pdf ####
  
  # Convert ggplots to grobs
  D046_DOC_grob <- ggplotGrob(D046_DOC_plot)
  C323_DOC_grob <- ggplotGrob(C323_DOC_plot)
  Navi_DOC_grob <- ggplotGrob(Navicula_DOC_plot)
  
  # Combine in grid
  DOC_gridplot <- cbind(D046_DOC_grob, Navi_DOC_grob, C323_DOC_grob, size = "first")
  
  grid.newpage()
  grid.draw(DOC_gridplot)
  
  # Save plot as pdf file
  pdf("Figures/FigS4_DOCvsTime.pdf", 
      width = 7.5, height = 3)
  grid.draw(DOC_gridplot)
  dev.off()
  
#### Fv/Fm  ####
  
  D046_FvFm_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !(is.na(growth_df_avgs$FvFm_mean)), ], 
                          aes(x = Day, y = FvFm_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=FvFm_mean - FvFm_sd, ymax=FvFm_mean + FvFm_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expand_scale(mult = c(0.05,0.1)), lim = c(0,0.7)) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 2)) +
    annotate("text", x = 0.2, y = 0.7, label = expression(paste(bold("A"))), size = 4) + 
    labs( y = expression(paste("F"["v"]*"/F"["m"])), x = "Days", title = expression(paste(italic("  Chlorella"), " sp. D046"))) +
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.title = element_blank(),
           legend.text = element_text(size = 10),
           axis.title = element_text(margin = margin(r = 15), size = 11),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 12),
           axis.text.y = element_text(size = 12))   
  
  Navicula_FvFm_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !(is.na(growth_df_avgs$FvFm_mean)), ], 
                              aes(x = Day, y = FvFm_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=FvFm_mean - FvFm_sd, ymax=FvFm_mean + FvFm_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expand_scale(mult = c(0.05,0.1)), lim = c(0, 0.7)) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 2)) +
    annotate("text", x = 0.2, y = 0.7, label = expression(paste(bold("B"))), size = 4) + 
    labs(x = "Days", title = expression(paste(italic("  Navicula"), " sp. SFP"))) +
    theme( legend.key = element_rect(fill = NA),  
           legend.title = element_blank(),
           legend.text = element_text(size = 10),
           axis.title.x = element_text(margin = margin(r = 15), size = 11),
           axis.title.y = element_blank(),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 12),
           axis.text.y = element_blank())
  
  C323_FvFm_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !(is.na(growth_df_avgs$FvFm_mean)) & growth_df_avgs$Round == 3, ], 
                          aes(x = Day, y = FvFm_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=FvFm_mean - FvFm_sd, ymax=FvFm_mean + FvFm_sd), width= 0.2) +
    scale_color_manual(labels = c("Fresh", "Reused"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Reused"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expand_scale(mult = c(0.05,0.1)), lim = c(0, 0.7)) +  
    scale_x_continuous(limits = c(-0.1,8.1), breaks = seq(0, 8, 2)) +
    annotate("text", x = 0.2, y = 0.7, label = expression(paste(bold("C"))), size = 4) + 
    labs(x = "Days", title = expression(paste(italic("Staurosira"), " sp. C323"))) +
    theme( legend.key = element_rect(fill = NA),  
           legend.title = element_blank(),
           legend.text = element_text(size = 10),
           axis.title.x = element_text(margin = margin(r = 15), size = 11),
           axis.title.y = element_blank(),
           legend.position = "top",
           legend.key.size = unit(4, "mm"),
           legend.justification = "center",
           legend.margin=margin(t = 0, b = 0, l = 4),
           legend.box.margin=margin(t = 0, b = 0),
           legend.spacing.x = unit(1, "mm"),
           plot.title = element_text(hjust = 0.5, size = 11),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA), 
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank(),
           axis.text.x = element_text(size = 12),
           axis.text.y = element_blank())
  
  # Combine Round 3 Daily plots in a grid & Save as pdf ####
  
  # Convert ggplots to grobs
  D046_FvFm_grob <- ggplotGrob(D046_FvFm_plot)
  C323_FvFm_grob <- ggplotGrob(C323_FvFm_plot)
  Navi_FvFm_grob <- ggplotGrob(Navicula_FvFm_plot)
  
  # Combine in grid
  FvFm_gridplot <- cbind(D046_FvFm_grob, Navi_FvFm_grob, C323_FvFm_grob, size = "first")
  
  grid.newpage()
  grid.draw(FvFm_gridplot)
  
  # Save plot as pdf file
  pdf("Figures/FigS3_FvFmvsTime.pdf", 
      width = 7.5, height = 3)
  grid.draw(FvFm_gridplot)
  dev.off()
  
#### mu #### with algae on x-axis and r and f as colors of the bars
  