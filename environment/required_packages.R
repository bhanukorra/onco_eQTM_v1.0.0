# Install all R packages required by Onco-eQTM pipeline

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_packages <- c("data.table", "MatrixEQTL", "dplyr", "RSQLite", "DBI", "matrixStats")
bioc_packages <- c("IlluminaHumanMethylation450kanno.ilmn12.hg19", "minfi")

install.packages(cran_packages[!cran_packages %in% installed.packages()[,"Package"]])
BiocManager::install(bioc_packages[!bioc_packages %in% installed.packages()[,"Package"]])

cat("All packages installed.\n")
