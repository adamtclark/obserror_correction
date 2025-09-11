setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
#load("output/diversity_e2.rda")

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

plot_sites_fun = function(agdat, ps = NULL, nameslist = c("div_1", "div_2", "div_3", "div_4"), niter = 1e3) {
  #require(lmodel2)
  
  colnames(agdat)[which(colnames(agdat) %in% nameslist)] = nameslist
  
  if(is.null(ps)) {
    ps = 1:nrow(agdat)
  }
  xsq = seq(min(agdat$div_1[ps], na.rm=TRUE), max(agdat$div_1[ps], na.rm=TRUE), length=100)
  
  obs_1 = rep(agdat$div_1[ps], 3)
  obs_2 = c(agdat$div_2[ps], agdat$div_3[ps], agdat$div_4[ps])
  moddat = data.frame(obs_2 = obs_2, obs_1 = obs_1)
  
  predmat = matrix(nrow = length(xsq), ncol = niter)
  slopemat = numeric(niter)
  r2lst = numeric(niter)
  for(i in 1:niter) {
    ps_iter = sample(1:nrow(moddat), replace = TRUE)
    mod = nondirfit(moddat[ps_iter,])
    predmat[,i] = (-mod$pars[1]-xsq*mod$pars[2])/mod$pars[3]
    #obs2 = (-intercept - obs1*p1)/p2
    r2lst[i] = mod$rsq$r_square_simple
    slopemat[i] = -mod$pars[2]/mod$pars[3]
    
    if(i/20 == floor(i/20)) {
      cat("\r", round(i/niter, 2))
    }
    if(i == niter) {
      cat("\r", "\n")
    }
  }
  return(list(xsq = xsq, predmat = predmat, slopemat=slopemat, r2lst=r2lst))
}

source("util/nondirfit.R")

set.seed(1234)
alpha_level = 0.8
cex_level = 0.8

# get values for plotting
# richness, all species
d_resurvey_ag_RI_same = agfun(0, survey_type = "same")
d_resurvey_ag_RI_diff = agfun(0, survey_type = "different")

mod_same_rich = plot_sites_fun(d_resurvey_ag_RI_same)
pred_same_rich = t(apply(mod_same_rich$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_rich = quantile(mod_same_rich$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_rich = quantile(mod_same_rich$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_rich = mod_same_rich$xsq

mod_diff_rich = plot_sites_fun(d_resurvey_ag_RI_diff)
pred_diff_rich = t(apply(mod_diff_rich$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_rich = quantile(mod_diff_rich$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_rich = quantile(mod_diff_rich$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_rich = mod_diff_rich$xsq

# richness, graminoids
d_resurvey_ag_RI_same_g = agfun(0, survey_type = "same", ps_agg = which(d_resurvey$group=="Graminoid"))
d_resurvey_ag_RI_diff_g = agfun(0, survey_type = "different", ps_agg = which(d_resurvey$group=="Graminoid"))

mod_same_rich_g = plot_sites_fun(d_resurvey_ag_RI_same_g)
pred_same_rich_g = t(apply(mod_same_rich_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_rich_g = quantile(mod_same_rich_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_rich_g = quantile(mod_same_rich_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_rich_g = mod_same_rich_g$xsq

mod_diff_rich_g = plot_sites_fun(d_resurvey_ag_RI_diff_g)
pred_diff_rich_g = t(apply(mod_diff_rich_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_rich_g = quantile(mod_diff_rich_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_rich_g = quantile(mod_diff_rich_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_rich_g = mod_diff_rich_g$xsq

# richness, non-graminoids
d_resurvey_ag_RI_same_ng = agfun(0, survey_type = "same", ps_agg = which(d_resurvey$group=="Non-Graminoid"))
d_resurvey_ag_RI_diff_ng = agfun(0, survey_type = "different", ps_agg = which(d_resurvey$group=="Non-Graminoid"))

mod_same_rich_ng = plot_sites_fun(d_resurvey_ag_RI_same_ng)
pred_same_rich_ng = t(apply(mod_same_rich_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_rich_ng = quantile(mod_same_rich_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_rich_ng = quantile(mod_same_rich_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_rich_ng = mod_same_rich_ng$xsq

mod_diff_rich_ng = plot_sites_fun(d_resurvey_ag_RI_diff_ng)
pred_diff_rich_ng = t(apply(mod_diff_rich_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_rich_ng = quantile(mod_diff_rich_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_rich_ng = quantile(mod_diff_rich_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_rich_ng = mod_diff_rich_ng$xsq

# shannon, all species
d_resurvey_ag_SH_same = agfun(1, survey_type = "same")
d_resurvey_ag_SH_diff = agfun(1, survey_type = "different")

mod_same_shannon = plot_sites_fun(d_resurvey_ag_SH_same)
pred_same_shannon = t(apply(mod_same_shannon$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_shannon = quantile(mod_same_shannon$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_shannon = quantile(mod_same_shannon$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_shannon = mod_same_shannon$xsq

mod_diff_shannon = plot_sites_fun(d_resurvey_ag_SH_diff)
pred_diff_shannon = t(apply(mod_diff_shannon$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_shannon = quantile(mod_diff_shannon$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_shannon = quantile(mod_diff_shannon$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_shannon = mod_diff_shannon$xsq

# shannon, graminoids
d_resurvey_ag_SH_same_g = agfun(1, survey_type = "same", ps_agg = which(d_resurvey$group=="Graminoid"))
d_resurvey_ag_SH_diff_g = agfun(1, survey_type = "different", ps_agg = which(d_resurvey$group=="Graminoid"))

mod_same_shannon_g = plot_sites_fun(d_resurvey_ag_SH_same_g)
pred_same_shannon_g = t(apply(mod_same_shannon_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_shannon_g = quantile(mod_same_shannon_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_shannon_g = quantile(mod_same_shannon_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_shannon_g = mod_same_shannon_g$xsq

mod_diff_shannon_g = plot_sites_fun(d_resurvey_ag_SH_diff_g)
pred_diff_shannon_g = t(apply(mod_diff_shannon_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_shannon_g = quantile(mod_diff_shannon_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_shannon_g = quantile(mod_diff_shannon_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_shannon_g = mod_diff_shannon_g$xsq

# shannon, non-graminoids
d_resurvey_ag_SH_same_ng = agfun(1, survey_type = "same", ps_agg = which(d_resurvey$group=="Non-Graminoid"))
d_resurvey_ag_SH_diff_ng = agfun(1, survey_type = "different", ps_agg = which(d_resurvey$group=="Non-Graminoid"))

mod_same_shannon_ng = plot_sites_fun(d_resurvey_ag_SH_same_ng)
pred_same_shannon_ng = t(apply(mod_same_shannon_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_shannon_ng = quantile(mod_same_shannon_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_shannon_ng = quantile(mod_same_shannon_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_shannon_ng = mod_same_shannon_ng$xsq

mod_diff_shannon_ng = plot_sites_fun(d_resurvey_ag_SH_diff_ng)
pred_diff_shannon_ng = t(apply(mod_diff_shannon_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_shannon_ng = quantile(mod_diff_shannon_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_shannon_ng = quantile(mod_diff_shannon_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_shannon_ng = mod_diff_shannon_ng$xsq

# simpson, all species
d_resurvey_ag_SI_same = agfun(2, survey_type = "same")
d_resurvey_ag_SI_diff = agfun(2, survey_type = "different")

mod_same_simpson = plot_sites_fun(d_resurvey_ag_SI_same)
pred_same_simpson = t(apply(mod_same_simpson$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_simpson = quantile(mod_same_simpson$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_simpson = quantile(mod_same_simpson$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_simpson = mod_same_simpson$xsq

mod_diff_simpson = plot_sites_fun(d_resurvey_ag_SI_diff)
pred_diff_simpson = t(apply(mod_diff_simpson$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_simpson = quantile(mod_diff_simpson$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_simpson = quantile(mod_diff_simpson$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_simpson = mod_diff_simpson$xsq

# simpson, graminoids
d_resurvey_ag_SI_same_g = agfun(2, survey_type = "same", ps_agg = which(d_resurvey$group=="Graminoid"))
d_resurvey_ag_SI_diff_g = agfun(2, survey_type = "different", ps_agg = which(d_resurvey$group=="Graminoid"))

mod_same_simpson_g = plot_sites_fun(d_resurvey_ag_SI_same_g)
pred_same_simpson_g = t(apply(mod_same_simpson_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_simpson_g = quantile(mod_same_simpson_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_simpson_g = quantile(mod_same_simpson_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_simpson_g = mod_same_simpson_g$xsq

mod_diff_simpson_g = plot_sites_fun(d_resurvey_ag_SI_diff_g)
pred_diff_simpson_g = t(apply(mod_diff_simpson_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_simpson_g = quantile(mod_diff_simpson_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_simpson_g = quantile(mod_diff_simpson_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_simpson_g = mod_diff_simpson_g$xsq

# simpson, non-graminoids
d_resurvey_ag_SI_same_ng = agfun(2, survey_type = "same", ps_agg = which(d_resurvey$group=="Non-Graminoid"))
d_resurvey_ag_SI_diff_ng = agfun(2, survey_type = "different", ps_agg = which(d_resurvey$group=="Non-Graminoid"))

mod_same_simpson_ng = plot_sites_fun(d_resurvey_ag_SI_same_ng)
pred_same_simpson_ng = t(apply(mod_same_simpson_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_simpson_ng = quantile(mod_same_simpson_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_simpson_ng = quantile(mod_same_simpson_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_simpson_ng = mod_same_simpson_ng$xsq

mod_diff_simpson_ng = plot_sites_fun(d_resurvey_ag_SI_diff_ng)
pred_diff_simpson_ng = t(apply(mod_diff_simpson_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_simpson_ng = quantile(mod_diff_simpson_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_simpson_ng = quantile(mod_diff_simpson_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_simpson_ng = mod_diff_simpson_ng$xsq





## make plot
pdf("figures/diversity_e2.pdf", width = 12, height = 12)
#png("figures/diversity_e2.png", width = 12, height = 12, units = "in", res = 200)

par(mfrow = c(3,3), mar = c(3,2,3,2), oma = c(2,2,0,0))
# richness, all species
pltlim = c(1, max(unlist(d_resurvey_ag_RI_same[,grep(pattern = "div_", x = colnames(d_resurvey_ag_RI_same), fixed = "TRUE")],
                    d_resurvey_ag_RI_diff[,grep(pattern = "div_", x = colnames(d_resurvey_ag_RI_diff), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_RI_same$div_1),
     jitter(d_resurvey_ag_RI_same$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Richness, repeat surveys", 2, outer = FALSE, line = 2.5)

points(jitter(d_resurvey_ag_RI_same$div_1), jitter(d_resurvey_ag_RI_same$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_RI_same$div_1), jitter(d_resurvey_ag_RI_same$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_RI_diff$div_1), jitter(d_resurvey_ag_RI_diff$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff$div_1), jitter(d_resurvey_ag_RI_diff$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff$div_1), jitter(d_resurvey_ag_RI_diff$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_rich, rev(xsq_same_rich)), c(pred_same_rich[,2], rev(pred_same_rich[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_rich, pred_same_rich[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_rich, rev(xsq_diff_rich)), c(pred_diff_rich[,2], rev(pred_diff_rich[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_rich, pred_diff_rich[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_rich, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_rich, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_rich[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_rich[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("A. All Species", adj = 0, line = 0.75)

# richness, graminoids
pltlim = c(1, max(unlist(d_resurvey_ag_RI_same_g[,grep(pattern = "div_", x = colnames(d_resurvey_ag_RI_same_g), fixed = "TRUE")],
                         d_resurvey_ag_RI_diff_g[,grep(pattern = "div_", x = colnames(d_resurvey_ag_RI_diff_g), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_RI_same_g$div_1),
     jitter(d_resurvey_ag_RI_same_g$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Richness, first survey", 1, line = 2.5)
points(jitter(d_resurvey_ag_RI_same_g$div_1), jitter(d_resurvey_ag_RI_same_g$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_RI_same_g$div_1), jitter(d_resurvey_ag_RI_same_g$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_RI_diff_g$div_1), jitter(d_resurvey_ag_RI_diff_g$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff_g$div_1), jitter(d_resurvey_ag_RI_diff_g$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff_g$div_1), jitter(d_resurvey_ag_RI_diff_g$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_rich_g, rev(xsq_same_rich_g)), c(pred_same_rich_g[,2], rev(pred_same_rich_g[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_rich_g, pred_same_rich_g[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_rich_g, rev(xsq_diff_rich_g)), c(pred_diff_rich_g[,2], rev(pred_diff_rich_g[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_rich_g, pred_diff_rich_g[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_rich_g, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_rich_g, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_rich_g[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_rich_g[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("B. Graminoids", adj = 0, line = 0.75)

# richness, non-graminoids
pltlim = c(1, max(unlist(d_resurvey_ag_RI_same_ng[,grep(pattern = "div_", x = colnames(d_resurvey_ag_RI_same_ng), fixed = "TRUE")],
                         d_resurvey_ag_RI_diff_ng[,grep(pattern = "div_", x = colnames(d_resurvey_ag_RI_diff_ng), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_RI_same_ng$div_1),
     jitter(d_resurvey_ag_RI_same_ng$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_RI_same_ng$div_1), jitter(d_resurvey_ag_RI_same_ng$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_RI_same_ng$div_1), jitter(d_resurvey_ag_RI_same_ng$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_RI_diff_ng$div_1), jitter(d_resurvey_ag_RI_diff_ng$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff_ng$div_1), jitter(d_resurvey_ag_RI_diff_ng$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_RI_diff_ng$div_1), jitter(d_resurvey_ag_RI_diff_ng$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_rich_ng, rev(xsq_same_rich_ng)), c(pred_same_rich_ng[,2], rev(pred_same_rich_ng[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_rich_ng, pred_same_rich_ng[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_rich_ng, rev(xsq_diff_rich_ng)), c(pred_diff_rich_ng[,2], rev(pred_diff_rich_ng[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_rich_ng, pred_diff_rich_ng[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_rich_ng, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_rich_ng, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_rich_ng[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_rich_ng[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("C. Non-Graminoids", adj = 0, line = 0.75)

# shannon, all species
pltlim = c(1, max(unlist(d_resurvey_ag_SH_same[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SH_same), fixed = "TRUE")],
                         d_resurvey_ag_SH_diff[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SH_diff), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_SH_same$div_1),
     jitter(d_resurvey_ag_SH_same$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Shannon, repeat surveys", 2, outer = FALSE, line = 2.5)

points(jitter(d_resurvey_ag_SH_same$div_1), jitter(d_resurvey_ag_SH_same$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SH_same$div_1), jitter(d_resurvey_ag_SH_same$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_SH_diff$div_1), jitter(d_resurvey_ag_SH_diff$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SH_diff$div_1), jitter(d_resurvey_ag_SH_diff$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SH_diff$div_1), jitter(d_resurvey_ag_SH_diff$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_shannon, rev(xsq_same_shannon)), c(pred_same_shannon[,2], rev(pred_same_shannon[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_shannon, pred_same_shannon[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_shannon, rev(xsq_diff_shannon)), c(pred_diff_shannon[,2], rev(pred_diff_shannon[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_shannon, pred_diff_shannon[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_shannon, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_shannon, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_shannon[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_shannon[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("D. All Species", adj = 0, line = 0.75)

# shannon, graminoids
pltlim = c(1, max(unlist(d_resurvey_ag_SH_same_g[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SH_same_g), fixed = "TRUE")],
                         d_resurvey_ag_SH_diff_g[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SH_diff_g), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_SH_same_g$div_1),
     jitter(d_resurvey_ag_SH_same_g$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Shannon, first survey", 1, line = 2.5)
points(jitter(d_resurvey_ag_SH_same_g$div_1), jitter(d_resurvey_ag_SH_same_g$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SH_same_g$div_1), jitter(d_resurvey_ag_SH_same_g$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_SH_diff_g$div_1), jitter(d_resurvey_ag_SH_diff_g$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SH_diff_g$div_1), jitter(d_resurvey_ag_SH_diff_g$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SH_diff_g$div_1), jitter(d_resurvey_ag_SH_diff_g$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_shannon_g, rev(xsq_same_shannon_g)), c(pred_same_shannon_g[,2], rev(pred_same_shannon_g[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_shannon_g, pred_same_shannon_g[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_shannon_g, rev(xsq_diff_shannon_g)), c(pred_diff_shannon_g[,2], rev(pred_diff_shannon_g[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_shannon_g, pred_diff_shannon_g[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_shannon_g, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_shannon_g, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_shannon_g[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_shannon_g[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("E. Graminoids", adj = 0, line = 0.75)

# shannon, non-graminoids
pltlim = c(1, max(unlist(d_resurvey_ag_SH_same_ng[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SH_same_ng), fixed = "TRUE")],
                         d_resurvey_ag_SH_diff_ng[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SH_diff_ng), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_SH_same_ng$div_1),
     jitter(d_resurvey_ag_SH_same_ng$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SH_same_ng$div_1), jitter(d_resurvey_ag_SH_same_ng$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SH_same_ng$div_1), jitter(d_resurvey_ag_SH_same_ng$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_SH_diff_ng$div_1), jitter(d_resurvey_ag_SH_diff_ng$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SH_diff_ng$div_1), jitter(d_resurvey_ag_SH_diff_ng$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SH_diff_ng$div_1), jitter(d_resurvey_ag_SH_diff_ng$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_shannon_ng, rev(xsq_same_shannon_ng)), c(pred_same_shannon_ng[,2], rev(pred_same_shannon_ng[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_shannon_ng, pred_same_shannon_ng[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_shannon_ng, rev(xsq_diff_shannon_ng)), c(pred_diff_shannon_ng[,2], rev(pred_diff_shannon_ng[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_shannon_ng, pred_diff_shannon_ng[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_shannon_ng, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_shannon_ng, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_shannon_ng[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_shannon_ng[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("F. Non-Graminoids", adj = 0, line = 0.75)


# simpson, all species
pltlim = c(1, max(unlist(d_resurvey_ag_SI_same[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SI_same), fixed = "TRUE")],
                         d_resurvey_ag_SI_diff[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SI_diff), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_SI_same$div_1),
     jitter(d_resurvey_ag_SI_same$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Simpson, repeat surveys", 2, outer = FALSE, line = 2.5)

points(jitter(d_resurvey_ag_SI_same$div_1), jitter(d_resurvey_ag_SI_same$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SI_same$div_1), jitter(d_resurvey_ag_SI_same$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_SI_diff$div_1), jitter(d_resurvey_ag_SI_diff$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SI_diff$div_1), jitter(d_resurvey_ag_SI_diff$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SI_diff$div_1), jitter(d_resurvey_ag_SI_diff$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_simpson, rev(xsq_same_simpson)), c(pred_same_simpson[,2], rev(pred_same_simpson[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_simpson, pred_same_simpson[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_simpson, rev(xsq_diff_simpson)), c(pred_diff_simpson[,2], rev(pred_diff_simpson[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_simpson, pred_diff_simpson[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_simpson, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_simpson, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_simpson[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_simpson[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("G. All Species", adj = 0, line = 0.75)

# simpson, graminoids
pltlim = c(1, max(unlist(d_resurvey_ag_SI_same_g[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SI_same_g), fixed = "TRUE")],
                         d_resurvey_ag_SI_diff_g[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SI_diff_g), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_SI_same_g$div_1),
     jitter(d_resurvey_ag_SI_same_g$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Simpson, first survey", 1, line = 2.5)
points(jitter(d_resurvey_ag_SI_same_g$div_1), jitter(d_resurvey_ag_SI_same_g$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SI_same_g$div_1), jitter(d_resurvey_ag_SI_same_g$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_SI_diff_g$div_1), jitter(d_resurvey_ag_SI_diff_g$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SI_diff_g$div_1), jitter(d_resurvey_ag_SI_diff_g$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SI_diff_g$div_1), jitter(d_resurvey_ag_SI_diff_g$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_simpson_g, rev(xsq_same_simpson_g)), c(pred_same_simpson_g[,2], rev(pred_same_simpson_g[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_simpson_g, pred_same_simpson_g[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_simpson_g, rev(xsq_diff_simpson_g)), c(pred_diff_simpson_g[,2], rev(pred_diff_simpson_g[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_simpson_g, pred_diff_simpson_g[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_simpson_g, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_simpson_g, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_simpson_g[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_simpson_g[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("H. Graminoids", adj = 0, line = 0.75)

# simpson, non-graminoids
pltlim = c(1, max(unlist(d_resurvey_ag_SI_same_ng[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SI_same_ng), fixed = "TRUE")],
                         d_resurvey_ag_SI_diff_ng[,grep(pattern = "div_", x = colnames(d_resurvey_ag_SI_diff_ng), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_SI_same_ng$div_1),
     jitter(d_resurvey_ag_SI_same_ng$div_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SI_same_ng$div_1), jitter(d_resurvey_ag_SI_same_ng$div_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_SI_same_ng$div_1), jitter(d_resurvey_ag_SI_same_ng$div_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_SI_diff_ng$div_1), jitter(d_resurvey_ag_SI_diff_ng$div_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SI_diff_ng$div_1), jitter(d_resurvey_ag_SI_diff_ng$div_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_SI_diff_ng$div_1), jitter(d_resurvey_ag_SI_diff_ng$div_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_simpson_ng, rev(xsq_same_simpson_ng)), c(pred_same_simpson_ng[,2], rev(pred_same_simpson_ng[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_simpson_ng, pred_same_simpson_ng[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_simpson_ng, rev(xsq_diff_simpson_ng)), c(pred_diff_simpson_ng[,2], rev(pred_diff_simpson_ng[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_simpson_ng, pred_diff_simpson_ng[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_simpson_ng, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_simpson_ng, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: E2 = ", round(r2_same_simpson_ng[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: E2 = ", round(r2_diff_simpson_ng[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("I. Non-Graminoids", adj = 0, line = 0.75)

dev.off()


#save(list = ls(), file = "output/diversity_e2.rda")
