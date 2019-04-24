# cross-strain_analyses

Data organizing, statistical analyses, and data visualization for experimental data and bacteria community data from cross-strain microalgae experiments in reused growth medium.

This repository contains 3 R scripts within an R project:
1) cross-strain_analyses_script.R carries out data wrangling and statistical analyses for experimental data related to microalgae growth, and produces figures from microalgae growth data. 
2) bacteria_community_analyses.R analyzes community-level and OTU-level data, using the vegan package, generated from sequencing bacteria communities associated with microalgae. It also generates NMDS plots of communities and a bacteria composition plot at the family taxonomic level.
3) deseq2_analyses.R calculates the differential abundance of bacteria OTUs between different treatments (fresh versus reused growth medium) using the DESeq2 package.

All data and metadata is provided in the input_data folder and is also citable on Figshare (https://doi.org/10.6084/m9.figshare.7831913.v1).

The R project also contains output folders "figures" and for "DESeq2_results."
