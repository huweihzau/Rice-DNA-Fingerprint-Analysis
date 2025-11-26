# Script: 0_install_dependencies.R
# Description: Install all R packages required for this project.

# List of CRAN packages
cran_packages <- c(
  "ggplot2",
  "dplyr",
  "corrplot",
  "pheatmap",
  "ape",
  "RColorBrewer",
  "factoextra",
  "reshape2"
)

# Install CRAN packages if not present
new_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# List of Bioconductor packages
bio_packages <- c("ggtree")

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in bio_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

message("All dependencies installed successfully.")