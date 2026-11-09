# Run once before the analysis; installation is separate from the research script.
packages <- c("vars", "tsDyn", "MTS", "TSA", "mvtnorm", "matrixcalc", "fastmatrix", "MASS", "dplyr", "portes", "urca", "rugarch", "aod", "readxl", "ggplot2", "fBasics", "fUnitRoots", "timeSeries", "svars")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
