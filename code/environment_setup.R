library(renv)
renv::init(project = "~/Documents/multires_bhicect/HiC_enrichment/HiC_zscore_explore/")

renv::install("tidyverse")
renv::install("furrr")
renv::install("Matrix")
renv::snapshot()
