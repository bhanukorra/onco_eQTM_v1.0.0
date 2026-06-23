# One-click installer for all R packages required by Onco-eQTM pipeline

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_packages <- c(
  "data.table",
  "MatrixEQTL",
  "dplyr",
  "RSQLite",
  "DBI",
  "matrixStats"
)

bioc_packages <- c(
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "minfi"
)

new_cran <- cran_packages[!cran_packages %in% installed.packages()[, "Package"]]
if (length(new_cran)) install.packages(new_cran)

new_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[, "Package"]]
if (length(new_bioc)) BiocManager::install(new_bioc)

cat("All packages installed successfully.\n")
