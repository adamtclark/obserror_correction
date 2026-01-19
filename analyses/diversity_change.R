setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
#source("util/save_model_outputs.R")
require(brms)
require(nlme)
load("output/brms_coverid_models_small.rda")
#load("output/diversity_change.rda")

set.seed(123432)

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.5
cex_level = 0.8
lwduse = 1.5

d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")
d = d[d$trt=="Control",]

# get cover data
covdat = d[,grep("cover_", colnames(d))]
covdat = covdat/100
d$cover = covdat[,1] # save first survey
d$nonzero_cov_mean = apply(covdat, 1, function(x) mean(x[!is.na(x) & x>0], na.rm=TRUE))

# get self- vs. non-self data
# same surveyor
surveyor = d[,grep("surveyor_", colnames(d))]
ps = surveyor!=surveyor[,1]
ps[is.na(ps)] = TRUE
covdat_same = covdat
covdat_same[ps] = NA
covdat_same[which(rowSums(covdat_same, na.rm=TRUE)==0),] = NA

d$ZeroObs_same = apply(covdat_same, 1, function(x) sum(!is.na(x) & x==0))
d$Ntrials_same = apply(covdat_same, 1, function(x) sum(!is.na(x)))-1
d$nonzero_cov_mean_same = apply(covdat_same, 1, function(x) mean(x[!is.na(x) & x>0], na.rm=TRUE))

# different surveyor
ps = surveyor==surveyor[,1]
ps[,-1][is.na(ps[,-1])] = TRUE
ps[,1][!is.na(ps[,1])] = FALSE
covdat_diff = covdat
covdat_diff[ps] = NA
covdat_diff[which(rowSums(covdat_diff, na.rm=TRUE)==0),] = NA

d$ZeroObs_diff = apply(covdat_diff, 1, function(x) sum(!is.na(x) & x==0))
d$Ntrials_diff = apply(covdat_diff, 1, function(x) sum(!is.na(x)))-1
d$nonzero_cov_mean_diff = apply(covdat_diff, 1, function(x) mean(x[!is.na(x) & x>0], na.rm=TRUE))

divdat = d[,c("SITE_CODE", "plot", "species", "cover", "group")]
divdat = divdat[!is.na(divdat$cover) & divdat$cover!=0,]

# simulate diversity change
simulate_change = function(divdat, deltaS = -1, process_noise = TRUE) {
  divdat_new = divdat
  
  # Add process noise following Taylor Power Law
  a = 0.01 # intercept
  b = 1.2  # scaling exponent
  
  if(process_noise) {
    proc_noise_sd = sqrt(a*divdat_new$cover^b)
    #plot(divdat_new$cover, proc_noise_sd); abline(a=0,b=0.1,lty=2)
    divdat_new$cover = pmin(abs(rnorm(nrow(divdat_new), mean = divdat_new$cover, sd = proc_noise_sd)),1)
  }
  
  if(deltaS==0) {
    return(divdat_new)
  } else {
    index_lst = paste(divdat$SITE_CODE, divdat$plot)
    uindex_lst = sort(unique(index_lst))
    
    for(i in 1:length(uindex_lst)) {
      ps = which(index_lst==uindex_lst[i])
      dtmp = divdat[ps,]
      S = length(ps)
      S_target = S+deltaS
      dS_percent = deltaS/S
      
      if(deltaS<0) { # species loss
        if(S_target > 0) {
          while(S != S_target) {
            #rnd = sample(x = which(dtmp$cover>0), size = 1)
            #dtmp$cover[rnd] = pmax(dtmp$cover[rnd] + dS_percent,0)
            
            rnd = sample(which(dtmp$cover > 0),1)
            dtmp$cover[rnd] = pmin(pmax(rnorm(1, dtmp$cover[rnd], sqrt(a*dtmp$cover[rnd]^b)),0),1)
            
            if(dtmp$cover[rnd] == 0) {
              S = S-1
            }
          }
        } else {
          # skip this entry, since deltaS is too small
          dtmp$cover = 0
        }
        divdat_new[ps,] = dtmp
      } else if(deltaS>0) { # species gain
        site = sort(unique(dtmp$SITE_CODE))
        
        new_species = divdat[divdat$SITE_CODE==site & !divdat$species %in% dtmp$species,]
        new_species = tapply(new_species$cover, new_species$species, mean)
        new_groups = divdat$group[match(names(new_species), divdat$species)]
        
        S_target = S+deltaS
        S_max = S + length(new_species)
        if(S_target <= S_max) {
          rnd = sample(x = length(new_species), size = deltaS)
          uplot = sort(unique(dtmp$plot))
          dnew = data.frame(SITE_CODE=site, plot = uplot,
                                  species = names(new_species[rnd]),
                                  cover = unname(new_species[rnd]),
                                  group = unname(new_groups[rnd]))
          divdat_new = rbind(divdat_new, dnew)
        } else {
          # skip this entry, since deltaS is too large
          dtmp$cover = 0
          divdat_new[ps,] = dtmp
        }
      }
    }
    # remove zeros
    divdat_new = divdat_new[divdat_new$cover>0,]
    
    #re-sort
    divdat_new = divdat_new[order(divdat_new$SITE_CODE, divdat_new$plot, divdat_new$species),]
    
    return(divdat_new)
    # check
    #checknum = 65
    #divdat[index_lst==uindex_lst[checknum],]
    #divdat_new[paste(divdat_new$SITE_CODE, divdat_new$plot)==uindex_lst[checknum],]
  }
}

# calculate diversity indices
hillfun = function(x, q = 0) {
  if(sum(x, na.rm=T)>0) {
    x[x==0] = NA
    if(is.null(dim(x))) {
      x = matrix(x)
    }
    if(nrow(x)>=1) {
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
  } else {
    return(0)
  }
}

# add noise to observations
simulate_noise = function(divdat) {
  divdat_out = divdat[,c("SITE_CODE", "plot", "species", "cover", "group")]
  
  ## make presence/absence error
  # make estimates for a single resurvey
  divdat$Ntrials_same = 1
  divdat$Ntrials_diff = 1
  
  # get pq for self-resurvey
  divdat$nonzero_cov_mean_same =  divdat$cover
  regression_same = colMeans(posterior_epred(mod_pgs, newdata = divdat, re_formula = NA))
  divdat$pzero_same = regression_same
  divdat$pq_hat_same = (regression_same)/(2-regression_same)
  
  # get q for non-self survey
  divdat$nonzero_cov_mean_diff =  divdat$cover
  regression_diff =colMeans(posterior_epred(mod_pd, newdata = divdat, re_formula = NA))
  divdat$pzero_diff = regression_diff
  divdat$p_hat_diff = with(divdat, (sqrt(2) *sqrt(-pq_hat_same* regression_diff^2 + 3 *pq_hat_same* regression_diff - 2 *pq_hat_same + regression_diff^2 - 3 *regression_diff + 2) + regression_diff - 2)/(regression_diff - 2))
  divdat$q_hat_diff = with(divdat, (regression_diff*(1-p_hat_diff^2+2*p_hat_diff)-(-2*p_hat_diff^2+4*p_hat_diff))/(-2+regression_diff*2)/p_hat_diff)
  
  # simulate noise
  # same surveyor, only missing errors
  error_same = rbinom(nrow(divdat), 1, prob = divdat$pq_hat_same)
  missed_species_same = error_same

  # different surveyor, missing errors AND ID errors
  error_diff = rbinom(nrow(divdat), 1, prob = divdat$p_hat_diff)
  tmp = rbinom(sum(error_diff), 1, prob = divdat$q_hat_diff[error_diff==1])
  missed_species_diff = id_error_diff = rep(0, nrow(divdat))
  missed_species_diff[error_diff==1][tmp==1] = 1
  id_error_diff[error_diff==1][tmp==0] = 1
  
  # mask missing species
  divdat_out$cover_resurvey_diff = divdat_out$cover_resurvey_same = divdat_out$cover
  divdat_out$cover_resurvey_same[missed_species_same==1] = 0
  divdat_out$cover_resurvey_diff[missed_species_diff==1] = 0
  
  # add in new rows for ID errors
  if(sum(id_error_diff)>0) {
    tmp = divdat_out[id_error_diff==1,]
    tmp$cover = tmp$cover_resurvey_same = 0
    tmp$species=paste(tmp$species, "_ID_error", sep = "")
    divdat_out$cover_resurvey_diff[id_error_diff==1] = 0
    divdat_out = rbind(divdat_out, tmp)
  }
  
  # re-sort dataframe
  divdat_out = divdat_out[order(divdat_out$SITE_CODE, divdat_out$plot, divdat_out$species),]
  
  ## add in cover error
  moddat=divdat_out
  moddat$cov_mean_same =  moddat$cover_resurvey_same
  moddat = moddat[moddat$cov_mean_same>0,]
  covermodel_same = colMeans(posterior_epred(mod_cgs, newdata = moddat, re_formula = NA))
  ps = divdat_out$cover_resurvey_same>0
  # recall model predicts CV, so updated value is (rnorm(CV)+1)*mu
  # setting minimum observation to 0.01%
  divdat_out$cover_resurvey_same[ps] =
    pmax(1e-04, (rnorm(length(covermodel_same), 0, covermodel_same)+1)*divdat_out$cover_resurvey_same[ps])
  
  moddat=divdat_out
  moddat$cov_mean_diff =  moddat$cover_resurvey_diff
  moddat = moddat[moddat$cov_mean_diff>0,]
  covermodel_diff = colMeans(posterior_epred(mod_cgd, newdata = moddat, re_formula = NA))
  ps = divdat_out$cover_resurvey_diff>0
  divdat_out$cover_resurvey_diff[ps] =
    pmax(1e-04, (rnorm(length(covermodel_diff), 0, covermodel_diff)+1)*divdat_out$cover_resurvey_diff[ps])
  
  # return output
  # new columns include estimates for
  # resurveys from same vs. different surveyor
  return(divdat_out)
}


get_pq_est =  function(divdat, type = "same", ntrials = 1, maskzeros = TRUE) {
  if(sum(colnames(divdat)=="cover")>0) {
    oldcover = divdat$cover
  } else {
    oldcover = NULL
  }
  
  if(type == "same") {
    divdat$cover = divdat$cover_resurvey_same
  } else if(type == "diff") {
    divdat$cover = divdat$cover_resurvey_diff
  } else {
    return("Error: type must be 'same' or 'diff'.")
  }
  
  # get pq for self-resurvey
  divdat$nonzero_cov_mean_same =  divdat$cover
  divdat$Ntrials_same = ntrials
  regression_same = colMeans(posterior_epred(mod_pgs, newdata = divdat, re_formula = NA))
  divdat$pzero_same = regression_same
  divdat$pq_hat_same = (regression_same)/(2-regression_same)
  
  # get q for non-self survey
  divdat$nonzero_cov_mean_diff =  divdat$cover
  divdat$Ntrials_diff = ntrials
  regression_diff =colMeans(posterior_epred(mod_pd, newdata = divdat, re_formula = NA))
  divdat$pzero_diff = regression_diff
  divdat$p_hat_diff = with(divdat, (sqrt(2) *sqrt(-pq_hat_same* regression_diff^2 + 3 *pq_hat_same* regression_diff - 2 *pq_hat_same + regression_diff^2 - 3 *regression_diff + 2) + regression_diff - 2)/(regression_diff - 2))
  divdat$q_hat_diff = with(divdat, (regression_diff*(1-p_hat_diff^2+2*p_hat_diff)-(-2*p_hat_diff^2+4*p_hat_diff))/(-2+regression_diff*2)/p_hat_diff)

  divdat$nonzero_cov_mean_same = divdat$Ntrials_same = NULL
  divdat$nonzero_cov_mean_diff = divdat$Ntrials_diff = NULL
  divdat$pzero_same = divdat$pzero_diff = NULL
  
  if(type == "same") {
    divdat$pi = divdat$pq_hat_same
    divdat$qi = 1
  } else if(type == "diff") {
    divdat$pi = divdat$p_hat_diff
    divdat$qi = divdat$q_hat_diff
  }
  
  if(maskzeros) {
    divdat$pi[divdat$cover==0] = NA
    divdat$qi[divdat$cover==0] = NA
  }
  
  divdat$p_hat_diff=divdat$q_hat_diff = NULL
  divdat$pq_hat_same = NULL
  
  divdat$cover = oldcover

  return(divdat)
}

# calculate alpha diversity
#divdat_out=simulate_noise(divdat)

get_alpha = function(divdat_out, q = 0) {
  alpha_out = unique(divdat_out[,c("SITE_CODE", "plot")])
  alpha_out$alpha_different = alpha_out$alpha_same = alpha_out$alpha_true = NA
  
  indexa = paste(alpha_out$SITE_CODE, alpha_out$plot)
  indexd = paste(divdat_out$SITE_CODE, divdat_out$plot)
  uindex = sort(unique(indexa))

  for(i in 1:length(uindex)) {
    ps = which(indexd==uindex[i])
    ps_a = which(indexa==uindex[i])
    alpha_out$alpha_true[ps_a] = hillfun(divdat_out$cover[ps], q = q)
    alpha_out$alpha_same[ps_a] = hillfun(divdat_out$cover_resurvey_same[ps], q = q)
    alpha_out$alpha_different[ps_a] = hillfun(divdat_out$cover_resurvey_diff[ps], q = q)
  }
  
  return(alpha_out)
}

# calculate turnover between two samples
#divdat_out_1=simulate_noise(divdat)
#divdat_out_2=simulate_noise(divdat)

get_beta = function(divdat_out_1, divdat_out_2, q = 0) {
  beta_out = unique(rbind(divdat_out_1[,c("SITE_CODE", "plot")],
                           divdat_out_2[,c("SITE_CODE", "plot")]))
  beta_out$beta_different = beta_out$beta_same = beta_out$beta_true = NA
  beta_out$gamma_different = beta_out$gamma_same = beta_out$gamma_true = NA
  beta_out$alpha_1_different = beta_out$alpha_1_same = beta_out$alpha_1_true = NA
  beta_out$alpha_2_different = beta_out$alpha_2_same = beta_out$alpha_2_true = NA
  
  # create combined dataset
  divdat_out_agg = with(rbind(divdat_out_1, divdat_out_2),
       aggregate(x = list(cover=cover, cover_resurvey_same=cover_resurvey_same, cover_resurvey_diff=cover_resurvey_diff),
                 by = list(SITE_CODE=SITE_CODE, plot=plot, species=species, group=group),
                 FUN = function(x) sum(x, na.rm=TRUE)))
  
  indexb = paste(beta_out$SITE_CODE, beta_out$plot)
  indexd_1 = paste(divdat_out_1$SITE_CODE, divdat_out_1$plot)
  indexd_2 = paste(divdat_out_2$SITE_CODE, divdat_out_2$plot)
  indexd_agg = paste(divdat_out_agg$SITE_CODE, divdat_out_agg$plot)
  uindex = sort(unique(indexb))
  
  no_comparison_possible = rep(0, length(uindex))
  
  for(i in 1:length(uindex)) {
    ps_1 = which(indexd_1==uindex[i])
    ps_2 = which(indexd_2==uindex[i])
    ps_b = which(indexb==uindex[i])
    ps_agg = which(indexd_agg==uindex[i])
    
    if(all(length(ps_1)>0,
           length(ps_2)>0,
           length(ps_b)>0,
           length(ps_agg)>0)) {
      alpha_true_1 = hillfun(divdat_out_1$cover[ps_1], q = q)
      alpha_same_1 = hillfun(divdat_out_1$cover_resurvey_same[ps_1], q = q)
      alpha_different_1 = hillfun(divdat_out_1$cover_resurvey_diff[ps_1], q = q)
      
      alpha_true_2 = hillfun(divdat_out_2$cover[ps_2], q = q)
      alpha_same_2 = hillfun(divdat_out_2$cover_resurvey_same[ps_2], q = q)
      alpha_different_2 = hillfun(divdat_out_2$cover_resurvey_diff[ps_2], q = q)
      
      gamma_true = hillfun(divdat_out_agg$cover[ps_agg], q = q)
      gamma_same = hillfun(divdat_out_agg$cover_resurvey_same[ps_agg], q = q)
      gamma_different = hillfun(divdat_out_agg$cover_resurvey_diff[ps_agg], q = q)
      
      beta_out$beta_true[ps_b] = gamma_true/((alpha_true_1+alpha_true_2)/2)
      beta_out$beta_same[ps_b] = gamma_same/((alpha_same_1+alpha_same_2)/2)
      beta_out$beta_different[ps_b] = gamma_different/((alpha_different_1+alpha_different_2)/2)
      
      beta_out$gamma_true[ps_b] = gamma_true
      beta_out$gamma_same[ps_b] = gamma_same
      beta_out$gamma_different[ps_b] = gamma_different
      
      beta_out$alpha_1_true[ps_b] = alpha_true_1
      beta_out$alpha_1_same[ps_b] = alpha_same_1
      beta_out$alpha_1_different[ps_b] = alpha_different_1
      
      beta_out$alpha_2_true[ps_b] = alpha_true_2
      beta_out$alpha_2_same[ps_b] = alpha_same_2
      beta_out$alpha_2_different[ps_b] = alpha_different_2
    } else {
      no_comparison_possible[ps_b] = 1
    }
  }
  beta_out = beta_out[no_comparison_possible==0,]
  
  return(beta_out)
}

# loop: test observed change across deltaS values
deltaS_lvls = -10:10

simout_rich = NULL
simout_shannon = NULL
simout_simpson = NULL
for(i in 1:length(deltaS_lvls)) {
  divdat_new = simulate_change(divdat, deltaS = deltaS_lvls[i])
  
  divdat_0 = simulate_noise(divdat)
  divdat_1 = simulate_noise(divdat_new)
  
  beta_out_rich = get_beta(divdat_out_1 = divdat_0, divdat_out_2 = divdat_1, q = 0)
  simout_rich = rbind(simout_rich,
                 data.frame(beta_out_rich, dS = deltaS_lvls[i]))
  
  beta_out_shannon = get_beta(divdat_out_1 = divdat_0, divdat_out_2 = divdat_1, q = 1)
  simout_shannon = rbind(simout_shannon,
                      data.frame(beta_out_shannon, dS = deltaS_lvls[i]))
  
  beta_out_simpson = get_beta(divdat_out_1 = divdat_0, divdat_out_2 = divdat_1, q = 2)
  simout_simpson = rbind(simout_simpson,
                      data.frame(beta_out_simpson, dS = deltaS_lvls[i]))
  
  cat("\r", round(i/length(deltaS_lvls),2))
}

## analyse and plot results
plot_delta_function = function(simout, dobs, colu = 1, ylim=c(-1,1), doplot=TRUE) {
  simout$derror = dobs
  mod = lme(derror~-1+as.factor(dS),
            data = simout,
            random = ~1|SITE_CODE/plot)
  cfs = summary(mod)$tTable
  dS_lst = sort(unique(simout$dS))
  
  if(doplot) {
    plot(0, 0,
         type = "n",,
         xlim = range(dS_lst),
         ylim = ylim,
         ylab = "", xlab = "")
    abline(h=0, lty=2)
    abline(v=0,lty=3)
  }
  polygon(c(dS_lst, rev(dS_lst)),
          c(cfs[,1]+cfs[,2]*qnorm(0.025),
            rev(cfs[,1]+cfs[,2]*qnorm(0.975))),
          col = adjustcolor(colu,alpha.f = alpha_level),
          border = NA)
  lines(dS_lst, cfs[,1], lwd = lwduse,
        col = adjustcolor(colu,alpha.f = alpha_level))
}

annotate_plot = function(xlm = c(-10, 10), ylm = c(-0.14, 0.14)) {
  text(xlm[1], ylm[2], "underestimate species loss", cex = 0.8, pos = 4)
  text(xlm[2], ylm[1], "underestimate species gain", cex = 0.8, pos = 2)
  text(xlm[1], ylm[1], "overestimate species loss", cex = 0.8, pos = 4)
  text(xlm[2], ylm[2], "overestimate species gain", cex = 0.8, pos = 2)
}

annotate_plot2 = function(xlm = c(-10, 10), ylm = c(-0.04, 0.14)) {
  text(xlm[2], ylm[2], "overestimate turnover", cex = 0.8, pos = 2)
  text(xlm[2], ylm[1], "underestimate turnover", cex = 0.8, pos = 2)
}

pdf("figures/diversity_change.pdf", width = 7, height = 6)
#png("figures/diversity_change.png", width = 7, height = 6, units = "in", res = 200)

par(mfrow = c(3,2), mar=c(2,2,1.5,1), oma = c(2,2,0,0))
lm_alpha = c(-0.15, 0.15)
lm_beta = c(-0.05, 0.15)

#TOTAL:
#rich
dtrue = (simout_rich$alpha_2_true-simout_rich$alpha_1_true)
dobs_same = ((simout_rich$alpha_2_same-simout_rich$alpha_1_same)-dtrue)/simout_rich$alpha_1_same
plot_delta_function(simout = simout_rich,dobs = dobs_same, ylim = lm_alpha, colu = collst[3])
annotate_plot()

dobs_diff = ((simout_rich$alpha_2_diff-simout_rich$alpha_1_diff)-dtrue)/simout_rich$alpha_1_diff
plot_delta_function(simout = simout_rich,dobs = dobs_diff,ylim=lm_alpha, doplot = FALSE, colu = collst[4])
title("A. Richness, alpha", adj = 0)

dtrue = (simout_rich$beta_true)
dobs_same = (simout_rich$beta_same)-dtrue
plot_delta_function(simout = simout_rich,dobs = dobs_same, ylim = lm_beta, colu = collst[3])
annotate_plot2()

dobs_diff = (simout_rich$beta_different)-dtrue
plot_delta_function(simout = simout_rich,dobs = dobs_diff, ylim = lm_beta, colu = collst[4], doplot = FALSE)
title("D. Richness, beta", adj = 0)

# shannon
dtrue = (simout_shannon$alpha_2_true-simout_shannon$alpha_1_true)
dobs_same = ((simout_shannon$alpha_2_same-simout_shannon$alpha_1_same)-dtrue)/simout_shannon$alpha_1_same
plot_delta_function(simout = simout_shannon,dobs = dobs_same, ylim = lm_alpha, colu = collst[3])
annotate_plot()

dobs_diff = ((simout_shannon$alpha_2_diff-simout_shannon$alpha_1_diff)-dtrue)/simout_shannon$alpha_1_diff
plot_delta_function(simout = simout_shannon,dobs = dobs_diff,ylim=lm_alpha, colu = collst[4], doplot = FALSE)
title("B. Shannon, alpha", adj = 0)

dtrue = (simout_shannon$beta_true)
dobs_same = (simout_shannon$beta_same)-dtrue
plot_delta_function(simout = simout_shannon,dobs = dobs_same, ylim = lm_beta, colu = collst[3])
annotate_plot2()

dobs_diff = (simout_shannon$beta_different)-dtrue
plot_delta_function(simout = simout_shannon,dobs = dobs_diff, ylim = lm_beta, colu = collst[4], doplot = FALSE)
title("E. Shannon, beta", adj = 0)

# simpson
dtrue = (simout_simpson$alpha_2_true-simout_simpson$alpha_1_true)
dobs_same = ((simout_simpson$alpha_2_same-simout_simpson$alpha_1_same)-dtrue)/simout_simpson$alpha_1_same
plot_delta_function(simout = simout_simpson,dobs = dobs_same, ylim = lm_alpha, colu = collst[3])
annotate_plot()

dobs_diff = ((simout_simpson$alpha_2_diff-simout_simpson$alpha_1_diff)-dtrue)/simout_simpson$alpha_1_diff
plot_delta_function(simout = simout_simpson,dobs = dobs_diff,ylim=lm_alpha, colu = collst[4], doplot = FALSE)
title("C. Simpson, alpha", adj = 0)

dtrue = (simout_simpson$beta_true)
dobs_same = (simout_simpson$beta_same)-dtrue
plot_delta_function(simout = simout_simpson,dobs = dobs_same, ylim = lm_beta, colu = collst[3])
annotate_plot2()

dobs_diff = (simout_simpson$beta_different)-dtrue
plot_delta_function(simout = simout_simpson,dobs = dobs_diff, ylim = lm_beta, colu = collst[4], doplot = FALSE)
title("F. Simpson, beta", adj = 0)

mtext("Relative bias", side = 2, outer = TRUE, line = 0.6)
mtext("Simulated change in richness", side = 1, outer = TRUE, line = 0.6)

legend(-10, .16, c("Same Surveyor", "Different Surveyors"),
       #fill = adjustcolor(collst[3:4], alpha_level),
       lty = 1, col = collst[3:4],
       bty = "n", cex = 1.2)

dev.off()


### now try MCMC fitting example
set.seed(23432)

## First, fix alpha diversity estimates
# sample data
deltaS = -2
true_div0 = divdat[divdat$SITE_CODE=="kzbg.at.dragnet" & divdat$plot == 5,]
true_div1 = simulate_change(true_div0, deltaS = deltaS)

obs_div0 = simulate_noise(divdat = true_div0)
obs_div1 = simulate_noise(divdat = true_div1)

# check diversity estimates
beta_test = get_beta(obs_div0, obs_div1)
beta_test[,grep("beta", names(beta_test))]

## run MCMC routine

# merge
obs_div = unique(rbind(obs_div0[,c("SITE_CODE", "plot", "species", "group")],
                       obs_div1[,c("SITE_CODE", "plot", "species", "group")]))
obs_div_merged1 = obs_div_merged0 = obs_div
obs_div_merged0$cover_resurvey_diff = obs_div_merged0$cover_resurvey_same = obs_div_merged0$cover = NA
obs_div_merged1$cover_resurvey_diff = obs_div_merged1$cover_resurvey_same = obs_div_merged1$cover = NA

ps = match(paste(obs_div_merged0$SITE_CODE, obs_div_merged0$plot, obs_div_merged0$species, obs_div_merged0$group),
           paste(obs_div0$SITE_CODE, obs_div0$plot, obs_div0$species, obs_div0$group))
obs_div_merged0$cover = obs_div0$cover[ps]
obs_div_merged0$cover_resurvey_same = obs_div0$cover_resurvey_same[ps]
obs_div_merged0$cover_resurvey_diff = obs_div0$cover_resurvey_diff[ps]
obs_div_merged0[is.na(obs_div_merged0)] = 0

ps = match(paste(obs_div_merged1$SITE_CODE, obs_div_merged1$plot, obs_div_merged1$species, obs_div_merged1$group),
           paste(obs_div1$SITE_CODE, obs_div1$plot, obs_div1$species, obs_div1$group))
obs_div_merged1$cover = obs_div1$cover[ps]
obs_div_merged1$cover_resurvey_same = obs_div1$cover_resurvey_same[ps]
obs_div_merged1$cover_resurvey_diff = obs_div1$cover_resurvey_diff[ps]
obs_div_merged1[is.na(obs_div_merged1)] = 0

# get likelihoods
obs_div_same0 = get_pq_est(obs_div_merged0, type = "same")
obs_div_diff0 = get_pq_est(obs_div_merged0, type = "diff")
obs_div_same1 = get_pq_est(obs_div_merged1, type = "same")
obs_div_diff1 = get_pq_est(obs_div_merged1, type = "diff")

## Next, fix beta diversity estimates
p_lose = abs(pmin(deltaS,0))/sum(!is.na(obs_div_same0$cover) & obs_div_same0$cover>0)

# losses
sps = 3 # 3 is obs error, 6 is lost
pi = obs_div_diff0$pi[sps]
qi = obs_div_diff0$qi[sps]
obs_div_diff0$cover[sps]
(1-pi)*p_lose # 1 found and 2 was lost
(1-pi)*pi # 1 found and 2 was missed or misided

(1-pi)*pi*qi # 1 found and 2 missed
(1-pi)*pi*(1-qi) # 1 found and 2 mis-ided

# mis-ids
