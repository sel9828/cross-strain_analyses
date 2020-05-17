# cross-strain_analyses

Data organizing, statistical analyses, and data visualization for experimental data and bacteria community data from cross-strain microalgae experiments in reused growth medium.

This repository contains 4 R scripts within an R project:
1) cross-strain_analyses_script.R carries out data wrangling and statistical analyses for experimental data related to microalgae growth, and produces figures from microalgae growth data. 
2) bacteria_community_analyses.R analyzes community-level and OTU-level data, using the vegan package, generated from sequencing bacteria communities associated with microalgae. It also generates NMDS plots of communities and a bacteria composition plot at the family taxonomic level.
3) deseq2_analyses.R calculates the differential abundance of bacteria OTUs between different treatments (fresh versus recycled growth medium) using the DESeq2 package.
4) deseq2_picrust_output.R is similar to deseq2_analyses.R, but calclulates the differential abundance of predicted metabolic pathways between treatments (fresh versus recycled growth medium). 

All data and metadata are provided in the input_data folder and are also citable from Figshare (https://doi.org/10.6084/m9.figshare.7831913.v3).

The R project also contains output folders "figures" and "DESeq2_results."
