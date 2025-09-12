setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
#source("util/save_model_outputs.R")
require(brms)
load("output/brms_coverid_models_small.rda")
#load("output/diversity_change.rda")

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.3
cex_level = 0.8
xsq = seq(0, 1, length = 1e3)

d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
d = d[d$trt=="Control",]



# count for each plot and surveyor pair:
# total number of cases with zero given nonzero
# frequency of shared congener
# frequency of shared other
# frequency of no shared

# work out probability
# pr. no shared = ?