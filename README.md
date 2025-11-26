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
