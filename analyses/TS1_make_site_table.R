# The following code creates a table with all sites used in the analysis (Table S1)

#setwd("~/Dropbox/Projects/117_ObservationError/src")
d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
rm(list = ls())
source("util/diversity_change_functions.R")

site_table = unique(d[,c("SITE", "SITE_CODE", "block", "plot")])
site_table$n_resurveys = NA
site_table$n_surveyors = NA

for(i in 1:nrow(site_table)) {
  ps = which(d$SITE_CODE==site_table$SITE_CODE[i] &
               d$plot==site_table$plot[i] & 
               d$block==site_table$block[i])
  
  tmp = unique(d[ps,c("surveyor_1", "surveyor_2", "surveyor_3", "surveyor_4")])
  site_table$n_resurveys[i] = sum(!is.na(tmp))
  site_table$n_surveyors[i] = length(unique(tmp[!is.na(tmp)]))
}

write.csv(site_table, "output/site_table.csv", row.names=FALSE)

table(site_table$SITE)
length(unique(site_table$SITE)) # number of sites (25)
nrow(site_table) # number of plots (91)

tmp = c(paste(d$SITE, d$plot, d$block, d$surveyor_1),
        paste(d$SITE, d$plot, d$block, d$surveyor_2),
        paste(d$SITE, d$plot, d$block, d$surveyor_3),
        paste(d$SITE, d$plot, d$block, d$surveyor_4))
tmp = tmp[-grep("NA", tmp)]
length(unique(tmp)) # number of surveys (143)

# format to number of plots and resurveys per site
site_table_summary = unique(d[,c("SITE", "SITE_CODE")])
site_table_summary$n_plots = NA
site_table_summary$n_samesurvey = NA
site_table_summary$n_differentsurvey = NA

for(i in 1:nrow(site_table_summary)) {
  ps = which(d$SITE_CODE==site_table_summary$SITE_CODE[i])
  
  tmp = unique(d[ps,c("block", "plot","surveyor_1", "surveyor_2", "surveyor_3", "surveyor_4")])
  site_table_summary$n_plots[i] = nrow(tmp)
  
  tot_surveys = apply(tmp[,-(1:2)],1,function(x) sum(!is.na(x)))
  unique_surveys = apply(tmp[,-(1:2)],1,function(x) length(unique(x[(!is.na(x))])))
  
  site_table_summary$n_samesurvey[i] = sum(tot_surveys)-sum(unique_surveys)
  site_table_summary$n_differentsurvey[i] = sum(tot_surveys)-site_table_summary$n_samesurvey[i]
}
site_table_summary = site_table_summary[order(site_table_summary$SITE),]
write.csv(site_table_summary, "output/site_table_summary.csv", row.names=FALSE)
