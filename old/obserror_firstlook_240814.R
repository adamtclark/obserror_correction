setwd("~/Dropbox/Projects/117_ObservationError/src/")
rm(list = ls())

d_resurvey = read.csv("data/NutNet_ReSurvey_processed_240814.csv")

plot_sites_fun = function(agdat, nameslist = c("div_1", "div_2", "div_3"), alphalevel = 0.2) {
  stlst = sort(unique(agdat$SITE_CODE))
  require(lmodel2)
  
  colnames(agdat)[which(colnames(agdat) %in% nameslist)] = c("div_1", "div_2", "div_3")
  
  slope_lst = matrix(nrow = length(stlst), ncol = 3, data = NA)
  
  suppressWarnings({
    for(i in 1:length(stlst)) {
      ps = which(agdat$SITE_CODE == stlst[i])
      xsq = seq(min(agdat$div_1[ps], na.rm=TRUE), max(agdat$div_1[ps], na.rm=TRUE), length=100)
      collst = viridis(3)
      collst2 = adjustcolor(viridis(3), alpha.f = alphalevel)
      
      if(sum(!is.na(agdat$div_2[ps]+agdat$div_1[ps]))>2) {
        mod1 = suppressMessages(lmodel2(div_2~div_1, agdat[ps,], "relative", "relative"))
        pd1 = mod1$regression.results[4,2]+mod1$regression.results[4,3]*xsq
        #predict(mod1, newdata = data.frame(div_1=xsq))
        lines(xsq, pd1, col = collst2[1])
        slope_lst[i,1] = mod1$regression.results[4,3]
      }
      
      if(sum(!is.na(agdat$div_3[ps]+agdat$div_1[ps]))>2) {
        mod2 = suppressMessages(lmodel2(div_3~div_1, agdat[ps,], "relative", "relative"))
        pd2 = mod2$regression.results[4,2]+mod2$regression.results[4,3]*xsq
        #predict(mod1, newdata = data.frame(div_1=xsq))
        lines(xsq, pd2, col = collst2[2])
        slope_lst[i,2] = mod2$regression.results[4,3]
      }
      if(sum(!is.na(agdat$div_4[ps]+agdat$div_1[ps]))>2) {
        mod3 = suppressMessages(lmodel2(div_4~div_1, agdat[ps,], "relative", "relative"))
        pd3 = mod3$regression.results[4,2]+mod3$regression.results[4,3]*xsq
        #predict(mod1, newdata = data.frame(div_1=xsq))
        lines(xsq, pd3, col = collst2[3])
        slope_lst[i,3] = mod3$regression.results[4,3]
      }
    }})
  return(slope_lst)
}

pdf("figures/species_cover_estimates.pdf", width = 4, height = 10)
  par(mar=c(4,4,2,2), mfrow = c(3,1))
  
  # plot
  require(viridis)
  collst = viridis(3)
  plot(d_resurvey$cover_1, d_resurvey$cover_2, xlab = "cover, first survey", ylab = "cover, repeat surveys", col = collst[1])
  points(d_resurvey$cover_1, d_resurvey$cover_3, col = collst[2])
  points(d_resurvey$cover_1, d_resurvey$cover_4, col = collst[3])
  abline(a=0, b = 1, lty = 2)
  
  slopes = plot_sites_fun(d_resurvey, c("cover_1", "cover_2", "cover_3"))
  legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
  #hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
  
  
  # estimate sd
  #sd_est=2.7
  #A1 = rnorm(1e5,0,sd_est); A2 = rnorm(1e5,0,sd_est)
  #sqrt(mean((A1-A2)^2)/2)
  #sqrt(mean(((A1-A2)^2)/2))
  #mean(((A1-A2)^2)/2)
  #sd_est^2
  
  mean_cover = rowMeans(d_resurvey[,grep("cover", colnames(d_resurvey))], na.rm=TRUE)
  var_cover = apply(d_resurvey[,grep("cover", colnames(d_resurvey))],1,function(x) var(x, na.rm=TRUE))
  
  mean_cover = mean_cover[var_cover>0]
  var_cover = var_cover[var_cover>0]
  
  plot(mean_cover, var_cover, log = "xy", xlab = "mean cover est.", ylab = "cover est. variance")
  mnsq = exp(seq(log(min(mean_cover, na.rm=TRUE)), log(max(mean_cover, na.rm=TRUE)), length=1e3))
  lines(mnsq, mnsq^2, lty = 2, lwd = 2, col = "red")
  
  lines(mnsq, (mnsq*0.1)^2, lty = 3, lwd = 2, col = "red")
  
  mod = loess(log(var_cover)~log(mean_cover))
  pred = predict(mod, newdata = data.frame(mean_cover = mnsq))
  
  lines(mnsq, exp(pred), col = "blue", lwd = 2)
  legend("topleft", c("sd(error) = mean est.", "sd(error) = 10% of mean est.", "loess fit"), lty = c(2,3,1), lwd = 2, col = c("red", "red", "blue"), bty = "n")
  
  
  
  
  
  # missing species
  prob_zero = apply(d_resurvey[,grep("cover", colnames(d_resurvey))], 1, function(x) mean(x[!is.na(x)]==0))
  n_obs = apply(d_resurvey[,grep("cover", colnames(d_resurvey))], 1, function(x) sum(!is.na(x)))
  nonzero_mean = apply(d_resurvey[,grep("cover", colnames(d_resurvey))], 1, function(x) mean(x[!is.na(x) & x>0]))
  zero_obs = n_obs*prob_zero
  
  plot(nonzero_mean, jitter(prob_zero), xlab = "mean non-zero cover", ylab = "probability of zero cover")
  
  mod = glm(cbind(zero_obs, n_obs)~nonzero_mean, family = "binomial")
  summary(mod)
  
  mnsq = (seq((min(nonzero_mean)), max((nonzero_mean)), length=1e3))
  pred = predict(mod, newdata = data.frame(nonzero_mean = mnsq), type = "response")
  
  lines(mnsq, pred, col = "blue", lwd = 2)
  legend("topright", c("logit regression"), lty = c(1), lwd = 2, col = c("blue"), bty = "n")
dev.off()



###################
### number of species
hillfun = function(x, q = 0) {
  x[x==0] = NA
  if(is.null(dim(x))) {
    x = matrix(x)
  }
  if(nrow(x)>1) {
    p = x/rep(colSums(x, na.rm=TRUE), each = nrow(x))
    if(q == 0) {
      D = colSums(x>0, na.rm=TRUE)
    } else if(q == 1) {
      D = exp(-colSums(p*log(p), na.rm = TRUE))
    } else {
      D = colSums(p^q,na.rm=TRUE)^(1/(1-q))
    }
    D[!is.finite(D)] = NA
    D[apply(x, 2, function(y) all(is.na(y)))] = NA
    return(D)
  } else {
    return(NA)
  }
}

## richness
d_resurvey_ag = with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
                                div_2 = cover_2,
                                div_3 = cover_3,
                                div_4 = cover_4),
                          by = list(SITE_CODE = SITE_CODE, trt = trt,
                                    first_nutrient_year = first_nutrient_year,
                                    block = block, plot = plot),
                          FUN = function(x) hillfun(x, q = 0)))
agfun = function(q) {
  with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
                     div_2 = cover_2,
                     div_3 = cover_3,
                     div_4 = cover_4),
                by = list(SITE_CODE = SITE_CODE, trt = trt,
                          first_nutrient_year = first_nutrient_year,
                          block = block, plot = plot),
                FUN = function(x) hillfun(x, q = q)))
}

d_resurvey_ag_RI = agfun(0)

collst = viridis(3)

pdf("figures/diversity_estimates.pdf", width = 10, height = 4)
  par(mfcol = c(2,3), mar = c(4,4,2,2))
  plot(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_2, xlab = "richness, first survey", ylab = "richness, repeat surveys", col = collst[1])
  points(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_3, col = collst[2])
  points(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_4, col = collst[3])
  abline(a=0, b = 1, lty = 2)
  slopes = plot_sites_fun(d_resurvey_ag_RI)
  legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
  hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
  abline(v = mean(slopes, na.rm=TRUE), lty = 2)
  
  
  ## shannon
  d_resurvey_ag_SH = agfun(1)
  
  collst = viridis(3)
  plot(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_2, xlab = "Shannon, first survey", ylab = "Shannon, repeat surveys", col = collst[1])
  points(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_3, col = collst[2])
  points(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_4, col = collst[3])
  abline(a=0, b = 1, lty = 2)
  slopes = plot_sites_fun(d_resurvey_ag_SH)
  legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
  hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
  abline(v = mean(slopes, na.rm=TRUE), lty = 2)
  
  ## Simpson
  d_resurvey_ag_SI = agfun(2)
  
  collst = viridis(3)
  plot(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_2, xlab = "Shannon, first survey", ylab = "Shannon, repeat surveys", col = collst[1])
  points(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_3, col = collst[2])
  points(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_4, col = collst[3])
  abline(a=0, b = 1, lty = 2)
  slopes = plot_sites_fun(d_resurvey_ag_SI)
  legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
  hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
  abline(v = mean(slopes, na.rm=TRUE), lty = 2)
  
  # note - higher error, due to impact of cover?
dev.off()

# total cover
pdf("figures/total_cover_estimates.pdf", width = 5, height = 4)
  par(mar=c(4,4,2,2), mfrow = c(1,1))
  
  cover_resurvey_ag = with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
                             div_2 = cover_2,
                             div_3 = cover_3,
                             div_4 = cover_4),
                        by = list(SITE_CODE = SITE_CODE, trt = trt,
                                  first_nutrient_year = first_nutrient_year,
                                  block = block, plot = plot),
                        FUN = function(x) sum(x, na.rm=TRUE)))
  cover_resurvey_ag[,6:9][cover_resurvey_ag[,6:9]==0] = NA
  
  collst = viridis(3)
  plot(cover_resurvey_ag$div_1, cover_resurvey_ag$div_2, xlab = "total cover, first survey", ylab = "total cover, repeat surveys", col = collst[1])
  points(cover_resurvey_ag$div_1, cover_resurvey_ag$div_3, col = collst[2])
  points(cover_resurvey_ag$div_1, cover_resurvey_ag$div_4, col = collst[3])
  abline(a=0, b = 1, lty = 2)
  legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
dev.off()


## add in lines per site
# average site-level R2, average cross-site R2

# error same vs. different people

# error vs. time between treatments

# error vs. nutrient or fencing treatment, or for high vs. low-diversity plots?

# confusion matrix?


