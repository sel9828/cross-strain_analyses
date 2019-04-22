# Code for creating Bray-Curtis NMDS plots and calculating community-wide changes (OTU richness, ANOSIM) from bacteria community OTU data.
# Author: Sarah Loftus, sarah.e.loftus@gmail.com

# Input data and metadata are also available on Figshare: https://doi.org/10.6084/m9.figshare.7831913.v1

# Load packages
library(vegan)
library(tidyverse)
library(RColorBrewer)
library(devtools)
#install_github("johannesbjork/LaCroixColoR")
library(LaCroixColoR)
library(gtable)
library(grid)


# Abolsute abundance OTU table
otutable <- read.delim("input_data/Data3_Absolute_OTUs.txt", skip = 1, header = T)  
otutable_rownames <- data.frame(otutable[,-1], row.names=otutable[,1])
otutable_rownames_f <- t(otutable_rownames[, 1:(ncol(otutable_rownames)-1)]) # transpose & remove taxonomy data (last column)

# Relative abundance OTU table
otutable_relative <- read.delim("input_data/Data4_Relative_OTUs.txt", skip = 1, header = T)  
otutable_rownames_relative <- data.frame(otutable_relative[,-1], row.names=otutable_relative[,1])
otutable_rownames_f_relative <- t(otutable_rownames_relative[, 1:(ncol(otutable_rownames_relative)-1)]) # transpose & remove taxonomy data (last column)

    # remove inoculum samples
    otutable_relative_no_inoc <- otutable_rownames_f_relative[-c(1, 8, 15), ]

# Relative abundance OTU table recalculated without the family Spirochaetacaceae (OTUs 16 and 27) 
no_spiro_relative <- read.delim("input_data/Data5_Relative_OTUs_noSpiro.txt", skip = 1, header = T)  
no_spiro_rownames_relative <- data.frame(no_spiro_relative[,-1], row.names=no_spiro_relative[,1])
no_spiro_rownames_f_relative <- t(no_spiro_rownames_relative[, 1:(ncol(no_spiro_rownames_relative)-1)]) # transpose & remove taxonomy data (last column)

    # remove inoculum samples
    no_spiro_relative_no_inoc <- no_spiro_rownames_f_relative[-c(1, 8, 15), ]


# Create a vector to use for differentiating points by color + shape on the figure
ID_all <- c("C323 Inoculum", rep("C323 Fresh", 3), rep("C323 Recycled", 3), 
            "D046 Inoculum", rep("D046 Fresh", 3), rep("D046 Recycled", 3), 
            "Navicula Inoculum", rep("Navicula Fresh", 3), rep("Navicula Recycled", 3))

ID_reps <- c("C323 Inoculum", "C323 Fresh A", "C323 Fresh B","C323 Fresh C","C323 Recycled A", "C323 Recycled B","C323 Recycled C",
             "D046 Inoculum", "D046 Fresh A", "D046 Fresh B", "D046 Fresh C","D046 Recycled A","D046 Recycled B","D046 Recycled C",
             "Navicula Inoculum", "Navicula Fresh A", "Navicula Fresh B","Navicula Fresh C","Navicula Recycled A", "Navicula Recycled B", "Navicula Recycled C")

ID_labels <- c("C323 Inoculum", "C323 Fresh","C323 Recycled",
               "D046 Inoculum", "D046 Fresh", "D046 Recycled",
               "Navicula Inoculum", "Navicula Fresh", "Navicula Recycled")

# Vectors of colors & shapes to use for plots
color_all <- c("gray79", rep('gray63', 3), rep('gray43', 3), 
                "darkolivegreen1",rep('green3', 3), rep('green4', 3),
                "lightblue3", rep('dodgerblue2', 3), rep('dodgerblue4', 3))

color <- c('gray63', "gray79", 'gray43', # Note: the order is different here to match the ggplot default of alphabetically ordering the IDs, so the Inoc and Fresh colors are switched
           'green3', "darkolivegreen1", 'green4',
           'dodgerblue2', "lightblue3", 'dodgerblue4') 

color_assign <- c("C323 Inoculum" = "gray79", "C323 Fresh" = 'gray63',"C323 Recycled" = 'gray43', 
                  "D046 Inoculum" = "greenyellow", "D046 Fresh" = 'green3', "D046 Recycled" = 'green4', 
                  "Navicula Inoculum" = "lightblue3", "Navicula Fresh" = 'dodgerblue2', "Navicula Recycled" = 'dodgerblue4') 

color_no_inoc <- c('gray63', 'gray43', 
                   'green3', 'green4', 
                   'dodgerblue2', 'dodgerblue4')

shape_all <- c(rep(16, 7), rep(17, 7), rep(15, 7))

shape <- c(rep(16, 3), rep(17, 3), rep(15, 3))

shape_no_inoc <- c(rep(16, 2), rep(17, 2), rep(15, 2))


# Legend text to use for plots
legend <- c(expression(paste(italic("Staurosira"), " sp. Inoculum")), expression(paste(italic("Staurosira"), " sp. Fresh")), expression(paste(italic("Staurosira"), " sp. Reused")), 
            expression(paste(italic("Chlorella"), " sp. Inoculum")), expression(paste(italic("Chlorella"), " sp. Fresh")), expression(paste(italic("Chlorella"), " sp. Reused")),
            expression(paste(italic("Navicula"), " sp. Inoculum")), expression(paste(italic("Navicula"), " sp. Fresh")), expression(paste(italic("Navicula"), " sp. Reused")))

### NMDS ###

# vegan package to calculate nMDS values using Relative Abundance data

bray_MDS <- metaMDS(comm = as.matrix(otutable_rownames_f_relative), distance = "bray", k = 2)

bray_MDS_noSpiro <- metaMDS(comm = as.matrix(no_spiro_rownames_f_relative), distance = "bray", k = 2)
      

# NMDS plots
      
# #ggplot (whole community, no ellipses)
# ggplot_NMDS = data.frame(NMDS1 = bray_MDS$points[,1], NMDS2 = bray_MDS$points[,2])      
# 
# ggplot(data = ggplot_NMDS, aes(x = NMDS1, y = NMDS2, color = ID_all, shape = ID_all)) +
#   geom_point(size = 4) +
#   scale_color_manual(breaks = ID_labels, values = color) +
#   scale_shape_manual(breaks = ID_labels, values = shape) +
#   theme( legend.background = element_rect(fill = "white"),  # removes color from behind legned points/lines
#          legend.key = element_rect(fill = "white"),
#          legend.title = element_blank(),
#          legend.text = element_text(size = 10),
#          legend.text.align = 0,
#          panel.background = element_rect(fill = NA),
#          panel.grid.major = element_line(colour = "gray87", size = 0.2),
#          panel.grid.minor = element_line(colour = "gray90", size = 0.2),
#          panel.border = element_rect(color = "gray60", fil = NA), 
#          axis.ticks = element_blank(),
#          axis.text = element_text(size = 12))

# points + ellipses - whole community
pdf("figures/Fig4_NMDS.pdf", width = 6, height = 4.5)

par(xpd = TRUE, mar = par()$mar + c(0,0,0,7)) # add room for legend
ellipse_plot <- ordiplot(bray_MDS, type = "none", las = 1, tck = 0.025)  # save the plot data in a plot object
points(ellipse_plot, "sites", pch= shape_all, col= color_all, cex = 1.1)  # add points, color & shape vectors that have a specific shape/color for each individual point
# add ellipses
ordiellipse(bray_MDS, groups = ID_all, display = "sites", draw = "lines",
            col = color, kind = "sd", conf = 0.95, lwd = 1.75)
legend(x = 1, y = 0.8, legend = legend, col = color_assign, pch = shape, cex = 0.8, pt.cex = 1.1, bty = "n") 
text(x= -1, y = 1, cex = 0.8, labels = paste("Stress =", format(bray_MDS$stress, digits=2)))

dev.off()

# points + ellipses - no Spirochaetes
pdf("figures/FigS6_NMDS_noSpiro.pdf", width = 6, height = 4.5)

par(xpd = TRUE, mar = par()$mar + c(0,0,0,7)) # add room for legend
ellipse_plot <- ordiplot(bray_MDS_noSpiro, type = "none", las = 1, tck = 0.025)  # save the plot data in a plot object
points(ellipse_plot, "sites", pch= shape_all, col= color_all, cex = 1.1)  # add points, color & shape vectors that have a specific shape/color for each individual point
# add ellipses
ordiellipse(bray_MDS_noSpiro, groups = ID_all, display = "sites", draw = "lines",
            col = color, kind = "sd", conf = 0.95, lwd = 1.75)
legend(x = 1, y = 0.8, legend = legend, col = color_assign, pch = shape, cex = 0.8, pt.cex = 1.1, bty = "n") 
text(x= -1, y = 1, cex = 0.8, labels = paste("Stress =", format(bray_MDS_noSpiro$stress, digits=2)))

dev.off()


### OTU richness ###

# species richness
spec_rich <- as.data.frame(specnumber(otutable_rownames_f))
spec_rich$sample <- ID_reps
names(spec_rich)[1] <- "spec_rich"

# reorder the levels
spec_rich$sample <- factor(spec_rich$sample, levels = ID_reps)

# plot bar graph
spec_rich_plot <- ggplot(data = spec_rich, aes(x = sample, y = spec_rich, fill = ID_all)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(breaks = ID_labels, values = color_assign, labels = legend) +
  geom_text(aes(label = spec_rich), position = position_dodge(width = 1), vjust = -0.5, size = 3.5) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05)), breaks = seq(0, 30, by = 5)) +
  labs(x = "Sample", y = "Number of OTUs") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.4), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text.align = 0)

pdf("figures/FigS5_OTU_richness.pdf", 
    width = 7, height = 4.5)
grid.draw(spec_rich_plot)
dev.off()
  

### Community level analyses ###

# Test Fresh vs Reused final communities
# Create data frame for treatment
treatment <- as.factor(c(rep("Fresh", 3), rep("Reused", 3)))

# Create permutation matrix to use for sample combinations
  # Number of possible combinations of the 6 samples is (6!/((6-3)!3!) )/2
  # Use the combn function to get the combinations, but append the final 10 columns to the first ten so that it is a combination of 6 instead of 3 (2 groups of 3)

combos <- combn(6,3)
perm_matrix <- t(rbind(combos[, 1:10], combos[, 20:11]))[2:10, ] # append the combos matrix & transpose so each row is a combination of 1-6; remove the first row because this is the given data set


# Anosim
C323.anosim <- anosim(otutable_rownames_f_relative[2:7, ], grouping = treatment, distance = "bray", permutations = perm_matrix)

D046.anosim <- anosim(otutable_rownames_f_relative[9:14, ], grouping = treatment, distance = "bray", permutations = perm_matrix)

Navi.anosim <- anosim(otutable_rownames_f_relative[16:21, ], grouping = treatment, distance = "bray", permutations = perm_matrix)


## Also test with taking out the Spirochetes since these were contaminants. 

# Anosim
C323.anosim_no_spiro <- anosim(no_spiro_rownames_f_relative[2:7, ], grouping = treatment, distance = "bray", permutations = perm_matrix)

D046.anosim_no_spiro <- anosim(no_spiro_rownames_f_relative[9:14, ], grouping = treatment, distance = "bray", permutations = perm_matrix)

Navi.anosim_no_spiro <- anosim(no_spiro_rownames_f_relative[16:21, ], grouping = treatment, distance = "bray", permutations = perm_matrix)


### Community composition plot at family level ###

# Read in fmaily classification table
otu_family_legend <- read.delim("input_data/RDP_OTU_family_key.txt", header = T)  

# Make absolute OTU table
otutable_absolute <- otutable[, 1:(ncol(otutable)-1)]

# Merge family classification and absolute otu table
otutable_absolute %>% 
  rename(OTU_ID = X.OTU.ID) %>% 
  inner_join(otu_family_legend, by = "OTU_ID") %>% 
  select(-taxonomy) %>% 
  gather(rownames(otutable_rownames_f), key = "sample", value = "abundance") -> otu_table_family # need to transform to tidy dataset first

# Sum Families within samples
otu_table_family %>% 
  select(-OTU_ID) %>% 
  group_by(sample, Family_Legend) %>% 
  summarize(fam_abundance = sum(abundance)) -> summed_families_by_sample

# Find the top families in the whole dataset. Choose to keep original name for taxa that make up at least 1% of the whole dataset, 
# rename the others as "other" OR group into phylum level Other group if they will collectively make up more than 1% as a collective other group
summed_families_by_sample %>% 
  select(-sample) %>% 
  group_by(Family_Legend) %>% 
  summarize(tot_abundance = sum(fam_abundance)) %>% 
  arrange(desc(tot_abundance))  %>% 
  mutate(rel_tot_abundance = 100*tot_abundance/sum(tot_abundance)) %>% # relative abundance in the total dataset
  mutate(legend_label = ifelse(rel_tot_abundance > 1, as.character(Family_Legend),  # Create new column for legend name
                               ifelse(Family_Legend == "HTCC2188" | Family_Legend == "Order_HTCC2188" | Family_Legend == "Order_Chromatiales" | Family_Legend == "Piscirickettsiaceae" | Family_Legend == "Saccharospirillaceae" , "Other_Gammaproteobacteria",
                                      ifelse(Family_Legend == "Pelagibacteraceae" | Family_Legend == "Hyphomicrobiaceae", "Other_Alphaproteobacteria", "Other")))) -> families_total_label 

    # Re-calculate abundances of new legend groupings to check if all > 1%
    families_total_label %>% 
      select(-tot_abundance, -Family_Legend) %>% 
      group_by(legend_label) %>% 
      summarize(rel_legend_abundance = sum(rel_tot_abundance)) -> legend_total_label
      # All groups are above 1% except the "Other" & "Other Gammaproteobacteria" category, which is OK because it's a catch-all for the low abundance taxa

# Merge the samples table with the new legend names. Sum the other categories
summed_families_by_sample %>% 
  inner_join(families_total_label, by = "Family_Legend") %>% 
  select(-tot_abundance, -rel_tot_abundance, -Family_Legend) %>% 
  group_by(sample, legend_label) %>% 
  summarize(figure_abundance = sum(fam_abundance)) -> summed_families_by_sample_legend

# Need to manually change names of samples and the order/names of the taxa

# reorder samples
samp_order <-  c("SL8", "SL9", "SL10", "SL11", "SL12", "SL13", "SL14", #D046
                  "SL15", "SL16", "SL17", "SL18", "SL19", "SL20", "SL21", # Navi
                  "SL1", "SL2", "SL3", "SL4", "SL5", "SL6", "SL7") #C323
summed_families_by_sample_legend$sample <- factor(summed_families_by_sample_legend$sample, levels = samp_order)

# legend names in order
legend_order <- c ( "Spirochaetaceae",
                   "Plantomycetaceae",
                   "Cryomorphaceae", "Balneolaceae", # Bacteroidetes
                   "Phyllobacteriaceae", "Hyphomonadaceae", "Rhodobacteraceae", "Rhodospirillaceae", "Erythrobacteraceae", "Other_Alphaproteobacteria", # Alphaproteobacteria
                   "Alteromonadaceae", "Other_Gammaproteobacteria", # Gammaproteobacteria
                   "Other") # Other (taxa less than 1% of total dataset and not belonging to a group below)

summed_families_by_sample_legend$legend_label <- factor(summed_families_by_sample_legend$legend_label, levels = legend_order)

# Updated legend names
legend_name <- c ( "Spirochaetaceae",
                   "Plantomycetaceae",
                   "Cryomorphaceae", "[Balneolaceae]", # Bacteroidetes
                   "Phyllobacteriaceae", "Hyphomonadaceae", "Rhodobacteraceae", "Rhodospirillaceae", "Erythrobacteraceae", "Other Alphaproteobacteria", # Alphaproteobacteria
                   "Alteromonadaceae", "Other Gammaproteobacteria", # Gammaproteobacteria
                   "Other Bacteria") # Other (taxa less than 1% of total dataset and not belonging to a group above) 
   

# Make color palette
colorCount <-  length(unique(summed_families_by_sample_legend$legend_label))

# Plot each algae separately then combine. Legend retains all families even if they're not in the sample.

# D046 (first plot)
D046_family_plot <- ggplot(data = filter(summed_families_by_sample_legend, sample %in% samp_order[1:7]),
                           aes(x = sample, y = figure_abundance, fill = legend_label, label = figure_abundance))+
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Relative Abundance") +
  scale_fill_manual(name = "Family",
                    labels = legend_name,
                    values = lacroix_palette(type = "paired", n = colorCount)) + 
  scale_y_continuous(labels = scales::percent, expand = expand_scale(mult = 0)) +
  theme(panel.background = element_blank(),
        axis.line = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# Navi (second plot)
Navi_family_plot <- ggplot(data = filter(summed_families_by_sample_legend, sample %in% samp_order[8:14]),
                           aes(x = sample, y = figure_abundance, fill = legend_label, label = figure_abundance))+
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Relative Abundance") +
  scale_fill_manual(name = "Family",
                    labels = legend_name,
                    values = lacroix_palette(type = "paired", n = colorCount)) + 
  scale_y_continuous(labels = scales::percent, expand = expand_scale(mult = 0)) +
  theme(panel.background = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

# C323 (third plot)
  C323_family_plot <- ggplot(data = filter(summed_families_by_sample_legend, sample %in% samp_order[15:21]),
                           aes(x = sample, y = figure_abundance, fill = legend_label, label = figure_abundance))+
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Relative Abundance") +
  scale_fill_manual(name = "Family",
                    labels = legend_name,
                    values = lacroix_palette(type = "paired", n = colorCount)) + 
  scale_y_continuous(labels = scales::percent, expand = expand_scale(mult = 0)) +
  theme(panel.background = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_blank(),
        legend.spacing.x = unit(1, "mm"),
        legend.key.size = unit(5, "mm"),
        axis.title = element_blank(),
        legend.title= element_text(size = 11),
        legend.title.align = 0.25)

# Combine plots
  # Convert ggplots to grobs
  D046_family_grob <- ggplotGrob(D046_family_plot)
  C323_family_grob <- ggplotGrob(C323_family_plot)
  Navi_family_grob <- ggplotGrob(Navi_family_plot)
  
  # Combine in grid
  family_gridplot <- cbind(D046_family_grob, Navi_family_grob, C323_family_grob, size = "first")
  
  grid.newpage()
  grid.draw(family_gridplot)
  
  # Save plot as pdf file. Edit x-axis using visual software. 
  pdf("figures/Fig3_family_composition.pdf", 
      width = 7.5, height = 4.5)
  grid.draw(family_gridplot)
  dev.off()  
  
  


