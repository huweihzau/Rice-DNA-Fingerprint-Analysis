# Genetic Diversity and Phenotypic Analysis of Shanlan Rice

This repository contains the R scripts and datasets used for the genetic and phenotypic analysis of Shanlan Rice resources. The project utilizes 38 InDel markers to assess genetic diversity, population structure, and phenotypic correlations among different rice accessions.

## 📂 Repository Structure

```text
.
├── data/
│   ├── genotype_of_38_InDel_markers.csv   # Binary genotype matrix for 38 InDel markers
│   └── phenotype.csv                      # Phenotypic traits data for correlation analysis
├── output/                                # Directory for storing generated plots and results
├── scripts/
│   ├── 0_install_dependencies.R           # Setup script to install required R packages
│   ├── 1_Correlation_coefficient.R        # Calculates and plots phenotypic correlations
│   ├── 2_DNA_fingerprint.R                # Visualizes DNA fingerprints of the population
│   ├── 3_Heatmap_of_SMC.R                 # Calculates Genetic Similarity (SMC) and plots heatmap
│   ├── 4_Principal_component_analysis.R   # Performs PCA based on genetic markers
│   ├── 5_Select_similar_varieties.R       # Filters and selects varieties based on genetic similarity
│   └── 6_Phylogenetic_tree.R              # Constructs UPGMA phylogenetic trees (Fan/Encircle style)
├── .gitignore                             # Specifies files to be ignored by Git
├── LICENSE                                # License information
├── README.md                              # Project documentation
└── requirements.txt                       # List of dependencies

🛠️ Prerequisites & Installation
The analysis is performed using R. Please ensure you have R (and RStudio recommended) installed.

Quick Setup
You can install all necessary dependencies by running the provided setup script:
R

source("scripts/0_install_dependencies.R")
Alternatively, install manually:
R

install.packages(c("ggplot2", "dplyr", "corrplot", "pheatmap", "ape", "RColorBrewer", "factoextra"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggtree")

🚀 Usage
The scripts are numbered to suggest a logical workflow, though they can be run independently.

1. Phenotypic Analysis
Script: scripts/1_Correlation_coefficient.R

Input: data/phenotype.csv

Output: Correlation matrix plots showing relationships between agronomic traits.

2. Genotypic Visualization
Script: scripts/2_DNA_fingerprint.R

Input: data/genotype_of_38_InDel markers.csv

Output: Visual representation of InDel marker distribution across varieties.

3. Genetic Similarity (SMC)
Script: scripts/3_Heatmap_of_SMC.R

Description: Calculates the Simple Matching Coefficient (SMC) between accessions and generates a heatmap to visualize genetic relatedness.

4. Population Structure (PCA)
Script: scripts/4_Principal_component_analysis.R

Output: PCA plots clustering the rice varieties based on their genotypic profiles.

5. Variety Selection
Script: scripts/5_Select_similar_varieties.R

Description: Identifies and filters rice varieties that share high genetic similarity.

6. Phylogenetic Analysis
Script: scripts/6_Phylogenetic_tree.R

Output: Generates UPGMA phylogenetic trees with subpopulation highlighting.
