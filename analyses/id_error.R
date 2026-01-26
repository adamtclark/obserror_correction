# The following code replicates Fig. 2 in the manuscript

#setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
require(brms)
#load("output/id_error.rda")

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.3
cex_level = 0.8
xsq = seq(0, 1, length = 1e3)

d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
d = d[d$trt=="Control",]

# get observation likelihood
covdat = d[,grep("cover_", colnames(d))]
covdat = covdat/100
d$cov_mean = apply(covdat, 1, function(x) mean(x, na.rm=TRUE))
d$nonzero_cov_mean = apply(covdat, 1, function(x) mean(x[!is.na(x) & x>0], na.rm=TRUE))

d$ZeroObs = apply(covdat, 1, function(x) sum(!is.na(x) & x==0))
d$Ntrials = apply(covdat, 1, function(x) sum(!is.na(x)))-1
ps = which(d$Ntrials>0)

# full model
mod_p0 <- brm(bf(ZeroObs | trials(Ntrials) ~ nonzero_cov_mean + (1|SITE_CODE)),
              #   zi ~ group + (1|SITE_CODE)),
              save_pars = save_pars(all = TRUE),
              data = d[ps,], family = binomial(), cores = 2)

mod_p1 <- brm(bf(ZeroObs | trials(Ntrials) ~ nonzero_cov_mean,
                 zi ~ (1|SITE_CODE)),
              save_pars = save_pars(all = TRUE),
              data = d[ps,], family = zero_inflated_binomial(), cores = 2)

mod_p2 <- brm(bf(ZeroObs | trials(Ntrials) ~ nonzero_cov_mean,
              zi ~ (1|SITE_CODE/plot)),
              save_pars = save_pars(all = TRUE),
              data = d[ps,], family = zero_inflated_binomial(), cores = 2)

loo_p0_p1_p2 = loo(mod_p0, mod_p1, mod_p2, moment_match = TRUE)
loo_p0_p1_p2 #p2 much better
mod_p0 = mod_p2
rm(mod_p2, mod_p1)

# graminoids vs. non-graminoids
mod_pg <- brm(bf(ZeroObs | trials(Ntrials) ~ nonzero_cov_mean*group,
                 zi ~ (1|SITE_CODE/plot)),
              save_pars = save_pars(all = TRUE),
              data = d[ps,], family = zero_inflated_binomial(), cores = 2)
loo_p0_pg = loo(mod_p0, mod_pg, moment_match = TRUE)
loo_p0_pg # p0 is better

# same vs. different surveyors
# same surveyor
surveyor = d[,grep("surveyor_", colnames(d))]
ps = surveyor!=surveyor[,1]
ps[is.na(ps)] = TRUE

covdat_same = covdat
covdat_same[ps] = NA

# remove cases with all zeros
covdat_same[which(rowSums(covdat_same, na.rm=TRUE)==0),] = NA

d$ZeroObs_same = apply(covdat_same, 1, function(x) sum(!is.na(x) & x==0))
d$Ntrials_same = apply(covdat_same, 1, function(x) sum(!is.na(x)))-1
d$nonzero_cov_mean_same = apply(covdat_same, 1, function(x) mean(x[!is.na(x) & x>0], na.rm=TRUE))
ps_same = which(d$Ntrials_same>0)

mod_ps <- brm(bf(ZeroObs_same | trials(Ntrials_same) ~ nonzero_cov_mean_same,
                 zi ~ (1|SITE_CODE/plot)),
              save_pars = save_pars(all = TRUE),
              data = d[ps_same,], family = zero_inflated_binomial(), cores = 2)

mod_pgs <- brm(bf(ZeroObs_same | trials(Ntrials_same) ~ nonzero_cov_mean_same*group,
                 zi ~ (1|SITE_CODE/plot)),
              save_pars = save_pars(all = TRUE),
              data = d[ps_same,], family = zero_inflated_binomial(), cores = 2)

loo_ps_pgs = loo(mod_ps, mod_pgs, moment_match = TRUE)
loo_ps_pgs
# pgs marginally better

# different surveyor
ps = surveyor==surveyor[,1]
ps[,-1][is.na(ps[,-1])] = TRUE
ps[,1][!is.na(ps[,1])] = FALSE

covdat_diff = covdat
covdat_diff[ps] = NA

# remove cases with all zeros
covdat_diff[which(rowSums(covdat_diff, na.rm=TRUE)==0),] = NA

d$ZeroObs_diff = apply(covdat_diff, 1, function(x) sum(!is.na(x) & x==0))
d$Ntrials_diff = apply(covdat_diff, 1, function(x) sum(!is.na(x)))-1
d$nonzero_cov_mean_diff = apply(covdat_diff, 1, function(x) mean(x[!is.na(x) & x>0], na.rm=TRUE))
ps_diff = which(d$Ntrials_diff>0)

mod_pd <- brm(bf(ZeroObs_diff | trials(Ntrials_diff) ~ nonzero_cov_mean_diff,
                 zi ~ (1|SITE_CODE/plot)),
              save_pars = save_pars(all = TRUE),
              data = d[ps_diff,], family = zero_inflated_binomial(), cores = 2)

mod_pgd <- brm(bf(ZeroObs_diff | trials(Ntrials_diff) ~ nonzero_cov_mean_diff*group,
                  zi ~ (1|SITE_CODE/plot)),
               save_pars = save_pars(all = TRUE),
               data = d[ps_diff,], family = zero_inflated_binomial(), cores = 2)

loo_pd_pgd = loo(mod_pd, mod_pgd, moment_match = TRUE)
loo_pd_pgd
# pd quite a bit better

# model predictions
ntrial = 1
tmp = posterior_epred(mod_p0, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean = xsq, Ntrials = ntrial))
pred_p0 = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_p0 = pred_p0/ntrial
colnames(pred_p0) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_pg, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean = xsq, Ntrials = ntrial, group = "Graminoid"))
pred_pg = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_pg = pred_pg/ntrial
colnames(pred_pg) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_pg, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean = xsq, Ntrials = ntrial, group = "Non-Graminoid"))
pred_png = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_png = pred_png/ntrial
colnames(pred_png) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_ps, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean_same = xsq, Ntrials_same = ntrial))
pred_ps = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_ps = pred_ps/ntrial
colnames(pred_ps) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_pgs, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean_same = xsq, Ntrials_same = ntrial, group = "Graminoid"))
pred_pgs = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_pgs = pred_pgs/ntrial
colnames(pred_pgs) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_pgs, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean_same = xsq, Ntrials_same = ntrial, group = "Non-Graminoid"))
pred_pngs = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_pngs = pred_pngs/ntrial
colnames(pred_pngs) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_pd, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean_diff = xsq, Ntrials_diff = ntrial))
pred_pd = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_pd = pred_pd/ntrial
colnames(pred_pd) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_pgd, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean_diff = xsq, Ntrials_diff = ntrial, group = "Graminoid"))
pred_pgd = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_pgd = pred_pgd/ntrial
colnames(pred_pgd) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

tmp = posterior_epred(mod_pgd, re_formula = ~ 0, newdata = data.frame(nonzero_cov_mean_diff = xsq, Ntrials_diff = ntrial, group = "Non-Graminoid"))
pred_pngd = t(apply(tmp, 2, function(x) c(mean(x), quantile(x, c(0.025, pnorm(c(-1,1)), 0.975)))))
pred_pngd = pred_pngd/ntrial
colnames(pred_pngd) = c("Estimate", "Q2.5", "QmSD", "QpSD", "Q97.5")

# plotting

# save models
#save(list = ls(), file = "output/id_error.rda")

#####################
# plot total
pdf("figures/id_error.pdf", width = 8, height = 4)
#png("figures/id_error.png", width = 8, height = 4, units = "in", res = 200)

X = cbind(c(1,1), c(1,1), c(1,1), c(2,3),  c(2,3))
layout(X)
par(mar=c(2,2,2,2), oma = c(2,2,0,0))
ps_same = which(d$Ntrials_same>0)
ps_diff = which(d$Ntrials_diff>0)
plot(d$nonzero_cov_mean_same[ps_same], (d$ZeroObs_same/d$Ntrials_same)[ps_same],
     xaxs = "i",
     ylim = c(0, 1),
     xlim = c(0,1),
     cex = cex_level, col = adjustcolor(collst[3], alpha_level),
     xlab = "", ylab = "")
mtext("Mean non-zero species-level cover", 1, line = 2.4)
mtext("Prob. detection error", 2, line = 2.4)
abline(h = c(0, 1), v = c(0,1), lty = 3)

polygon(c(xsq, rev(xsq)), c(pred_ps[,c(3)], rev(pred_ps[,c(4)])), col = adjustcolor(collst[3], 0.5), border = NA)
lines(xsq, pred_ps[,1], lty = c(1), col = collst[3], lwd = 1.8)

points((d$nonzero_cov_mean_diff[ps_diff]), (d$ZeroObs_diff/d$Ntrials_diff)[ps_diff],
       pch = 2,
       cex = cex_level, col = adjustcolor(collst[4], alpha_level))
polygon(c(xsq, rev(xsq)), c(pred_pd[,c(3)], rev(pred_pd[,c(4)])), col = adjustcolor(collst[4], 0.5), border = NA)
lines(xsq, pred_pd[,1], lty = c(1), col = collst[4], lwd = 1.8)
#matlines(xsq, pred_cd[,c(1,4,5)], lty = c(1,2,2), col = collst[4])

legend(0.6, 0.97, c("Same Surveyor", "Different Surveyors"),
       lty = 1, pch = c(1,2), col = collst[3:4], bty = "n", cex = 1.2)
title("A. All Species", adj = 0)

# add in mean lines
ctps = seq(0, 1, length=20)
ctps_plot = ctps[-length(ctps)]
ct_same = cut(d$nonzero_cov_mean_same[ps_same], ctps)
y_same = (d$ZeroObs_same/d$Ntrials_same)[ps_same]
mod_same = lm(y_same~ct_same, data = d[ps_same,])
pred_same = predict(mod_same, newdata = data.frame(ct_same = sort(unique(ct_same))))
tmp = table(ct_same)>0
lines(ctps_plot[tmp], pred_same, lwd = 1.5, lty = 2, col = collst[3])

ct_diff = cut(d$nonzero_cov_mean_diff[ps_diff], ctps)
y_diff = (d$ZeroObs_diff/d$Ntrials_diff)[ps_diff]
mod_diff = lm(y_diff~ct_diff, data = d[ps_diff,])
pred_diff = predict(mod_diff, newdata = data.frame(ct_diff = sort(unique(ct_diff))))
tmp = table(ct_diff)>0
lines(ctps_plot[tmp], pred_diff, lwd = 1.5, lty = 2, col = collst[4])


# plot graminoids
ps_same = which(d$Ntrials_same>0 & d$group == "Graminoid")
ps_diff = which(d$Ntrials_diff>0 & d$group == "Graminoid")

plot(d$nonzero_cov_mean_same[ps_same], (d$ZeroObs_same/d$Ntrials_same)[ps_same],
     xaxs = "i",
     ylim = c(0, 1),
     xlim = c(0,1),
     cex = cex_level, col = adjustcolor(collst[3], alpha_level),
     xlab = "", ylab = "")
abline(h = c(0, 1), v = c(0,1), lty = 3)

polygon(c(xsq, rev(xsq)), c(pred_pgs[,c(3)], rev(pred_pgs[,c(4)])), col = adjustcolor(collst[3], 0.5), border = NA)
lines(xsq, pred_pgs[,1], lty = c(1), col = collst[3], lwd = 1.8)

points((d$nonzero_cov_mean_diff[ps_diff]), (d$ZeroObs_diff/d$Ntrials_diff)[ps_diff],
       pch = 2,
       cex = cex_level, col = adjustcolor(collst[4], alpha_level))
polygon(c(xsq, rev(xsq)), c(pred_pd[,c(3)], rev(pred_pd[,c(4)])), col = adjustcolor(collst[4], 0.5), border = NA)
lines(xsq, pred_pd[,1], lty = c(1), col = collst[4], lwd = 1.8)

# add in mean lines
ctps = seq(0, 1, length=20)
ctps_plot = ctps[-length(ctps)]
ct_same = cut(d$nonzero_cov_mean_same[ps_same], ctps)
y_same = (d$ZeroObs_same/d$Ntrials_same)[ps_same]
mod_same = lm(y_same~ct_same, data = d[ps_same,])
pred_same = predict(mod_same, newdata = data.frame(ct_same = sort(unique(ct_same))))
tmp = table(ct_same)>0
lines(ctps_plot[tmp], pred_same, lwd = 1.5, lty = 2, col = collst[3])

ct_diff = cut(d$nonzero_cov_mean_diff[ps_diff], ctps)
y_diff = (d$ZeroObs_diff/d$Ntrials_diff)[ps_diff]
mod_diff = lm(y_diff~ct_diff, data = d[ps_diff,])
pred_diff = predict(mod_diff, newdata = data.frame(ct_diff = sort(unique(ct_diff))))
tmp = table(ct_diff)>0
lines(ctps_plot[tmp], pred_diff, lwd = 1.5, lty = 2, col = collst[4])

title("B. Graminoids", adj = 0)


# plot non-graminoids
ps_same = which(d$Ntrials_same>0 & d$group == "Non-Graminoid")
ps_diff = which(d$Ntrials_diff>0 & d$group == "Non-Graminoid")

plot(d$nonzero_cov_mean_same[ps_same], (d$ZeroObs_same/d$Ntrials_same)[ps_same],
     xaxs = "i",
     ylim = c(0, 1),
     xlim = c(0,1),
     cex = cex_level, col = adjustcolor(collst[3], alpha_level),
     xlab = "", ylab = "")
abline(h = c(0, 1), v = c(0,1), lty = 3)

polygon(c(xsq, rev(xsq)), c(pred_pngs[,c(3)], rev(pred_pngs[,c(4)])), col = adjustcolor(collst[3], 0.5), border = NA)
lines(xsq, pred_pngs[,1], lty = c(1), col = collst[3], lwd = 1.8)

points((d$nonzero_cov_mean_diff[ps_diff]), (d$ZeroObs_diff/d$Ntrials_diff)[ps_diff],
       pch = 2,
       cex = cex_level, col = adjustcolor(collst[4], alpha_level))
polygon(c(xsq, rev(xsq)), c(pred_pd[,c(3)], rev(pred_pd[,c(4)])), col = adjustcolor(collst[4], 0.5), border = NA)
lines(xsq, pred_pd[,1], lty = c(1), col = collst[4], lwd = 1.8)

# add in mean lines
ctps = seq(0, 1, length=20)
ctps_plot = ctps[-length(ctps)]
ct_same = cut(d$nonzero_cov_mean_same[ps_same], ctps)
y_same = (d$ZeroObs_same/d$Ntrials_same)[ps_same]
mod_same = lm(y_same~ct_same, data = d[ps_same,])
pred_same = predict(mod_same, newdata = data.frame(ct_same = sort(unique(ct_same))))
tmp = table(ct_same)>0
lines(ctps_plot[tmp], pred_same, lwd = 1.5, lty = 2, col = collst[3])

ct_diff = cut(d$nonzero_cov_mean_diff[ps_diff], ctps)
y_diff = (d$ZeroObs_diff/d$Ntrials_diff)[ps_diff]
mod_diff = lm(y_diff~ct_diff, data = d[ps_diff,])
pred_diff = predict(mod_diff, newdata = data.frame(ct_diff = sort(unique(ct_diff))))
tmp = table(ct_diff)>0
lines(ctps_plot[tmp], pred_diff, lwd = 1.5, lty = 2, col = collst[4])

title("C. Non-Graminoids", adj = 0)

mtext("Mean non-zero species-level cover", 1, line = 2.4)
mtext("Prob. detection error", 2, line = 2.4, adj = -1)

dev.off()