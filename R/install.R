source("https://bioconductor.org/biocLite.R")
biocLite(c("xcms", "MSnbase", "doParallel", "msdata", "magrittr",
           "devtools"))
## Need xcms version > 3.3.1
if (packageVersion("xcms") < "3.3.1")
  devtools::install_github("sneumann/xcms", ref = "master")
