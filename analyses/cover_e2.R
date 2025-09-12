setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
#load("output/cover_e2.rda")

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.3
cex_level = 0.8
xsq = seq(0, 1, length = 1e3)

d_resurvey = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
d_resurvey = d_resurvey[d_resurvey$trt=="Control",]

agfun = function(ps_agg = NULL, survey_type = "same") {
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
  
  agfun_sum = function(x) {
    if(all(is.na(x))) {
      NA
    } else {
      sum(x, na.rm=TRUE)
    }
  }
  
  with(tmp, aggregate(list(cover_1 = cover_1,
                           cover_2 = cover_2,
                           cover_3 = cover_3,
                           cover_4 = cover_4),
        by = list(SITE_CODE = SITE_CODE, trt = trt,
                  first_nutrient_year = first_nutrient_year,
                  block = block, plot = plot),
        FUN = agfun_sum))
}

plot_sites_fun = function(agdat, ps = NULL, nameslist = c("cover_1", "cover_2", "cover_3", "cover_4"), niter = 1e3) {
  #require(lmodel2)
  
  colnames(agdat)[which(colnames(agdat) %in% nameslist)] = nameslist
  
  if(is.null(ps)) {
    ps = 1:nrow(agdat)
  }
  xsq = seq(min(agdat$cover_1[ps], na.rm=TRUE), max(agdat$cover_1[ps], na.rm=TRUE), length=100)
  
  obs_1 = rep(agdat$cover_1[ps], 3)
  obs_2 = c(agdat$cover_2[ps], agdat$cover_3[ps], agdat$cover_4[ps])
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
d_resurvey_ag_COV_same = agfun(survey_type = "same")
d_resurvey_ag_COV_diff = agfun(survey_type = "different")

mod_same_cover = plot_sites_fun(d_resurvey_ag_COV_same)
pred_same_cover = t(apply(mod_same_cover$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_cover = quantile(mod_same_cover$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_cover = quantile(mod_same_cover$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_cover = mod_same_cover$xsq

mod_diff_cover = plot_sites_fun(d_resurvey_ag_COV_diff)
pred_diff_cover = t(apply(mod_diff_cover$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_cover = quantile(mod_diff_cover$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_cover = quantile(mod_diff_cover$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_cover = mod_diff_cover$xsq

# richness, graminoids
d_resurvey_ag_COV_same_g = agfun(survey_type = "same", ps_agg = which(d_resurvey$group=="Graminoid"))
d_resurvey_ag_COV_diff_g = agfun(survey_type = "different", ps_agg = which(d_resurvey$group=="Graminoid"))

mod_same_cover_g = plot_sites_fun(d_resurvey_ag_COV_same_g)
pred_same_cover_g = t(apply(mod_same_cover_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_cover_g = quantile(mod_same_cover_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_cover_g = quantile(mod_same_cover_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_cover_g = mod_same_cover_g$xsq

mod_diff_cover_g = plot_sites_fun(d_resurvey_ag_COV_diff_g)
pred_diff_cover_g = t(apply(mod_diff_cover_g$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_cover_g = quantile(mod_diff_cover_g$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_cover_g = quantile(mod_diff_cover_g$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_cover_g = mod_diff_cover_g$xsq

# richness, non-graminoids
d_resurvey_ag_COV_same_ng = agfun(survey_type = "same", ps_agg = which(d_resurvey$group=="Non-Graminoid"))
d_resurvey_ag_COV_diff_ng = agfun(survey_type = "different", ps_agg = which(d_resurvey$group=="Non-Graminoid"))

mod_same_cover_ng = plot_sites_fun(d_resurvey_ag_COV_same_ng)
pred_same_cover_ng = t(apply(mod_same_cover_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_same_cover_ng = quantile(mod_same_cover_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_same_cover_ng = quantile(mod_same_cover_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_same_cover_ng = mod_same_cover_ng$xsq

mod_diff_cover_ng = plot_sites_fun(d_resurvey_ag_COV_diff_ng)
pred_diff_cover_ng = t(apply(mod_diff_cover_ng$predmat, 1, function(x) quantile(x, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)))
slope_diff_cover_ng = quantile(mod_diff_cover_ng$slopemat, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
r2_diff_cover_ng = quantile(mod_diff_cover_ng$r2lst, c(0.025, pnorm(-1:1), 0.975), na.rm=TRUE)
xsq_diff_cover_ng = mod_diff_cover_ng$xsq




## make plot
pdf("figures/cover_e2.pdf", width = 12, height = 4)
#png("figures/cover_e2.png", width = 12, height = 4, units = "in", res = 200)

par(mfrow = c(1,3), mar = c(3,2,3,2), oma = c(2,2,0,0))
# richness, all species
pltlim = c(0, max(unlist(d_resurvey_ag_COV_same[,grep(pattern = "cover_", x = colnames(d_resurvey_ag_COV_same), fixed = "TRUE")],
                    d_resurvey_ag_COV_diff[,grep(pattern = "cover_", x = colnames(d_resurvey_ag_COV_diff), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_COV_same$cover_1),
     jitter(d_resurvey_ag_COV_same$cover_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Total cover, repeat surveys", 2, outer = FALSE, line = 2.5)

points(jitter(d_resurvey_ag_COV_same$cover_1), jitter(d_resurvey_ag_COV_same$cover_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_COV_same$cover_1), jitter(d_resurvey_ag_COV_same$cover_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_COV_diff$cover_1), jitter(d_resurvey_ag_COV_diff$cover_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_COV_diff$cover_1), jitter(d_resurvey_ag_COV_diff$cover_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_COV_diff$cover_1), jitter(d_resurvey_ag_COV_diff$cover_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_cover, rev(xsq_same_cover)), c(pred_same_cover[,2], rev(pred_same_cover[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_cover, pred_same_cover[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_cover, rev(xsq_diff_cover)), c(pred_diff_cover[,2], rev(pred_diff_cover[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_cover, pred_diff_cover[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_cover, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_cover, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: R-sq = ", round(r2_same_cover[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: R-sq = ", round(r2_diff_cover[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("A. All Species", adj = 0, line = 0.75)

# richness, graminoids
pltlim = c(0, max(unlist(d_resurvey_ag_COV_same_g[,grep(pattern = "cover_", x = colnames(d_resurvey_ag_COV_same_g), fixed = "TRUE")],
                         d_resurvey_ag_COV_diff_g[,grep(pattern = "cover_", x = colnames(d_resurvey_ag_COV_diff_g), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_COV_same_g$cover_1),
     jitter(d_resurvey_ag_COV_same_g$cover_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
mtext("Total cover, first survey", 1, line = 2.5)
points(jitter(d_resurvey_ag_COV_same_g$cover_1), jitter(d_resurvey_ag_COV_same_g$cover_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_COV_same_g$cover_1), jitter(d_resurvey_ag_COV_same_g$cover_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_COV_diff_g$cover_1), jitter(d_resurvey_ag_COV_diff_g$cover_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_COV_diff_g$cover_1), jitter(d_resurvey_ag_COV_diff_g$cover_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_COV_diff_g$cover_1), jitter(d_resurvey_ag_COV_diff_g$cover_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_cover_g, rev(xsq_same_cover_g)), c(pred_same_cover_g[,2], rev(pred_same_cover_g[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_cover_g, pred_same_cover_g[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_cover_g, rev(xsq_diff_cover_g)), c(pred_diff_cover_g[,2], rev(pred_diff_cover_g[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_cover_g, pred_diff_cover_g[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_cover_g, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_cover_g, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: R-sq = ", round(r2_same_cover_g[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: R-sq = ", round(r2_diff_cover_g[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("B. Graminoids", adj = 0, line = 0.75)

# richness, non-graminoids
pltlim = c(0, max(unlist(d_resurvey_ag_COV_same_ng[,grep(pattern = "cover_", x = colnames(d_resurvey_ag_COV_same_ng), fixed = "TRUE")],
                         d_resurvey_ag_COV_diff_ng[,grep(pattern = "cover_", x = colnames(d_resurvey_ag_COV_diff_ng), fixed = "TRUE")]), na.rm=TRUE))
plot(jitter(d_resurvey_ag_COV_same_ng$cover_1),
     jitter(d_resurvey_ag_COV_same_ng$cover_2),
     xlim = pltlim, ylim = pltlim,
     xlab = "", ylab = "",
     col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_COV_same_ng$cover_1), jitter(d_resurvey_ag_COV_same_ng$cover_3), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)
points(jitter(d_resurvey_ag_COV_same_ng$cover_1), jitter(d_resurvey_ag_COV_same_ng$cover_4), col = adjustcolor(collst[3], alpha.f = alpha_level), cex = cex_level)

points(jitter(d_resurvey_ag_COV_diff_ng$cover_1), jitter(d_resurvey_ag_COV_diff_ng$cover_2), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_COV_diff_ng$cover_1), jitter(d_resurvey_ag_COV_diff_ng$cover_3), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)
points(jitter(d_resurvey_ag_COV_diff_ng$cover_1), jitter(d_resurvey_ag_COV_diff_ng$cover_4), col = adjustcolor(collst[4], alpha.f = alpha_level), pch = 2, cex = cex_level)

abline(a=0, b = 1, lty = 2, lwd = 1.2)

polygon(c(xsq_same_cover_ng, rev(xsq_same_cover_ng)), c(pred_same_cover_ng[,2], rev(pred_same_cover_ng[,4])),
        col = adjustcolor(collst[3], alpha.f = 0.5), border = NA)
lines(xsq_same_cover_ng, pred_same_cover_ng[,3], col = collst[3], lwd = 1.4)

polygon(c(xsq_diff_cover_ng, rev(xsq_diff_cover_ng)), c(pred_diff_cover_ng[,2], rev(pred_diff_cover_ng[,4])),
        col = adjustcolor(collst[4], alpha.f = 0.5), border = NA)
lines(xsq_diff_cover_ng, pred_diff_cover_ng[,3], col = collst[4], lwd = 1.4)

tmp_same = unname(round(quantile(slope_same_cover_ng, c(0.025, 0.975)),2))
tmp_diff = unname(round(quantile(slope_diff_cover_ng, c(0.025, 0.975)),2))
legend("bottomright",
       c(paste("Same Surveyor: R-sq = ", round(r2_same_cover_ng[3],2), ", Slope [",tmp_same[1], ",", tmp_same[2] ,"]", sep = ""),
         paste("Different Surveyors: R-sq = ", round(r2_diff_cover_ng[3],2), ", Slope [",tmp_diff[1], ",", tmp_diff[2], "]", sep = "")),
       col = collst[3:4], lty = 1, pch = 1:2, cex = 1, bty = "n")
title("C. Non-Graminoids", adj = 0, line = 0.75)

dev.off()


#save(list = ls(), file = "output/cover_e2.rda")
