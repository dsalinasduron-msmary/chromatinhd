library(scCustomize)

args <- commandArgs(trailingOnly=TRUE)
integrated.multiome.rds <- readRDS(args[1])
WriteH5MU(integrated.multiome.rds,args[2])
