d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
times = unique(d[,c("date_1", "date_2", "date_3", "date_4")])

d_hours = rep(NA, nrow(times))
for(i in 1:ncol(times)) {
  for(j in 1:ncol(times)) {
    if(i!=j) {
      tmp = abs(difftime(as.Date(times[,i],"%d/%m/%Y"),
                         as.Date(times[,j],"%d/%m/%Y"),
                         units = "days"))
      tmp = abs(as.numeric(tmp))
      ps = !is.na(tmp)
      if(sum(ps)>0) {
        d_hours[ps] = pmax(d_hours[ps], tmp[ps], na.rm=TRUE)
      }
    }
  }
}

d_hours = d_hours[!is.na(d_hours)]
hist(d_hours)
mean(d_hours)
sd(d_hours)
quantile(d_hours, c(0.025, 0.5, 0.975))
