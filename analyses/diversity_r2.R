# TODO: permutation test for r2, and slope, and prediction interval?


setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
require(lmodel2)
#require(brms)

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.3
cex_level = 0.8
xsq = seq(0, 1, length = 1e3)

d_resurvey = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
d_resurvey = d_resurvey[d_resurvey$trt=="Control",]

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
agfun = function(q, ps_agg = NULL, survey_type = "same") {
  if(is.null(ps_agg)) {
    ps_agg = 1:nrow(d_resurvey)
  }
  tmp = d_resurvey[ps_agg,]
  if(survey_type=="same") {
    covdat = tmp[,grep("cover_", colnames(tmp))]
    covdat = covdat/100
    
    surveyor = tmp[,grep("surveyor_", colnames(tmp))]
    ps_cv = surveyor!=surveyor[,1]
    ps_cv[is.na(ps_cv)] = TRUE
    
    covdat_same = covdat
    covdat_same[ps_cv] = NA
    
    # remove cases with all zeros
    covdat_same[which(rowSums(covdat_same, na.rm=TRUE)==0),] = NA
    
    # add back to data
    tmp[,grep("cover_", colnames(tmp))] = covdat_same
  } else if(survey_type=="different") {
    covdat = tmp[,grep("cover_", colnames(tmp))]
    covdat = covdat/100
    
    surveyor = tmp[,grep("surveyor_", colnames(tmp))]
    ps_cv = surveyor==surveyor[,1]
    ps_cv[,-1][is.na(ps_cv[,-1])] = TRUE
    ps_cv[,1][!is.na(ps_cv[,1])] = FALSE
    
    covdat_diff = covdat
    covdat_diff[ps_cv] = NA
    
    # remove cases with all zeros
    covdat_diff[which(rowSums(covdat_diff, na.rm=TRUE)==0),] = NA
    
    # add back to data
    tmp[,grep("cover_", colnames(tmp))] = covdat_diff
  }
  
  with(tmp, aggregate(list(div_1 = cover_1,
             div_2 = cover_2,
             div_3 = cover_3,
             div_4 = cover_4),
        by = list(SITE_CODE = SITE_CODE, trt = trt,
                  first_nutrient_year = first_nutrient_year,
                  block = block, plot = plot),
        FUN = function(x) hillfun(x, q = q)))
}

plot_sites_fun = function(agdat, ps = NULL, nameslist = c("div_1", "div_2", "div_3", "div_4"), coluse = "black") {
  require(lmodel2)
  
  colnames(agdat)[which(colnames(agdat) %in% nameslist)] = nameslist
  
  if(is.null(ps)) {
    ps = 1:nrow(agdat)
  }
  xsq = seq(min(agdat$div_1[ps], na.rm=TRUE), max(agdat$div_1[ps], na.rm=TRUE), length=100)
  
  obs_1 = rep(agdat$div_1[ps], 3)
  obs_2 = c(agdat$div_2[ps], agdat$div_3[ps], agdat$div_4[ps])
  moddat = data.frame(obs_2 = obs_2, obs_1 = obs_1)
  
  mod = suppressMessages(lmodel2(obs_2~obs_1, moddat, "relative", "relative"))
  pd1 = mod$regression.results[4,2]+mod$regression.results[4,3]*xsq
  
  slope = unlist(c(mod$regression.results[4,3], mod$confidence.intervals[4,4:5]))
  
  # add in CI
  x_mu = mean(obs_1[!is.na(obs_1) & !is.na(obs_2)])
  y_mu = mean(obs_2[!is.na(obs_1) & !is.na(obs_2)])
  
  b0_low = y_mu-slope[2]*x_mu
  b0_hi = y_mu-slope[3]*x_mu
  
  pd2 = b0_low+slope[2]*xsq
  pd3 = b0_hi+slope[3]*xsq
  
  polygon(c(xsq, rev(xsq)), c(pd2, rev(pd3)), col = adjustcolor(coluse, alpha.f = 0.5), border = NA)
  lines(xsq, pd1, col = coluse, lwd = 1.4)
  
  r2 = mod$rsquare
  
  return(list(slope = slope, r2 = r2))
}

d_resurvey_ag_RI_same = agfun(0, survey_type = "same")
d_resurvey_ag_RI_diff = agfun(0, survey_type = "different")

set.seed(1234)
alpha_level = 0.8
cex_level = 0.8

#pdf("figures/diversity_estimates.pdf", width = 10, height = 4)
#par(mfcol = c(2,3), mar = c(4,4,2,2))
plot(jitter(d_resurvey_ag_RI_same$div_1),
     jitter(d_resurvey_ag_RI_same$div_2),
     xlab = "richness, first survey", ylab = "richness, repeat surveys",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_RI_same$div_1), jitter(d_resurvey_ag_RI_same$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_RI_same$div_1), jitter(d_resurvey_ag_RI_same$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_RI_diff$div_1), jitter(d_resurvey_ag_RI_diff$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff$div_1), jitter(d_resurvey_ag_RI_diff$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff$div_1), jitter(d_resurvey_ag_RI_diff$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

slopes_same = plot_sites_fun(d_resurvey_ag_RI_same, coluse = collst[3])
slopes_diff = plot_sites_fun(d_resurvey_ag_RI_diff, coluse = collst[4])

hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
abline(v = mean(slopes, na.rm=TRUE), lty = 2)

# add in R2

# graminoid vs. non-graminoid





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
