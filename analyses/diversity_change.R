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

divdat = d[,c("SITE_CODE", "plot", "species", "cover", "group", "nonzero_cov_mean")]
#divdat = divdat[!is.na(divdat$cover),]

# simulate diversity change
simulate_change = function(divdat, deltaS = -1) {
  divdat_new = divdat
  
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
            rnd = sample(x = which(dtmp$cover>0), size = 1)
            dtmp$cover[rnd] = pmax(dtmp$cover[rnd] + dS_percent,0)
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
        
        S_target = S+deltaS
        S_max = S + length(new_species)
        if(S_target <= S_max) {
          rnd = sample(x = length(new_species), size = deltaS)
          uplot = sort(unique(dtmp$plot))
          dnew = data.frame(SITE_CODE=site, plot = uplot,
                                  species = names(new_species[rnd]),
                                  cover = unname(new_species[rnd]))
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

# add noise to observations
simulate_noise = function(divdat) {
  # make estimates for a single resurvey
  divdat$Ntrials_same = 1
  divdat$Ntrials_diff = 1
  
  # get pq for self-resurvey
  divdat$nonzero_cov_mean_same =  divdat$nonzero_cov_mean
  regression_same = colMeans(posterior_epred(mod_pgs, newdata = divdat, re_formula = NA))
  divdat$pzero_same = regression_same
  divdat$pq_hat_same = (regression_same)/(2-regression_same)
  
  # get q for non-self survey
  divdat$nonzero_cov_mean_diff =  divdat$nonzero_cov_mean
  regression_diff =colMeans(posterior_epred(mod_pd, newdata = divdat, re_formula = NA))
  divdat$pzero_diff = regression_diff
  divdat$p_hat_diff = with(divdat, (sqrt(2) *sqrt(-pq_hat_same* regression_diff^2 + 3 *pq_hat_same* regression_diff - 2 *pq_hat_same + regression_diff^2 - 3 *regression_diff + 2) + regression_diff - 2)/(regression_diff - 2))
  divdat$q_hat_diff = with(divdat, (regression_diff*(1-p_hat_diff^2+2*p_hat_diff)-(-2*p_hat_diff^2+4*p_hat_diff))/(-2+regression_diff*2)/p_hat_diff)
  
  
  q_hat = (-(p_hat^2*(regression_diff - 2)) + 2*p_hat*(regression_diff - 2) + regression_diff)/(2*p_hat*(regression_diff - 1))
  
  hist(pmax(pmin(q_hat,1),0))
  
  # get p, shared, and alpha for each plot
  p = tapply(divdat$p_hat, paste(divdat$SITE_CODE, divdat$plot), mean)
  alpha_1 = tapply(covdat_diff[,1]>0, paste(divdat$SITE_CODE, divdat$plot), function(x) sum(x,na.rm=TRUE))
  alpha_2 = tapply(covdat_diff[,2]>0, paste(divdat$SITE_CODE, divdat$plot), function(x) sum(x,na.rm=TRUE))
  shared = tapply((covdat_diff[,1]>0) & (covdat_diff[,2]>0), paste(divdat$SITE_CODE, divdat$plot), function(x) sum(x,na.rm=TRUE))
  gamma = tapply((covdat_diff[,1]>0) | (covdat_diff[,2]>0), paste(divdat$SITE_CODE, divdat$plot), function(x) sum(x,na.rm=TRUE))
  ps = alpha_1 > 0 & alpha_2 > 0
  #plot(alpha_1[ps], alpha_2[ps])
  #abline(a=0,b=1,lty=3)
  #N_hat = shared[ps]/((1-p[ps])^2)
  #plot(N_hat, gamma[ps])
  #q_hat = (1-alpha_2[ps]/N_hat)/p[ps]
  R = tapply(regression_diff, paste(divdat$SITE_CODE, divdat$plot), mean)
  
  tmp = data.frame(shared = shared[ps],
             alpha_1 = alpha_1[ps],
             alpha_2 = alpha_2[ps],
             gamma = gamma[ps],
             p = round(p[ps],3),
             R = R[ps])
  
  #tmp$p_est_q1 = round((tmp$gamma-tmp$shared)/tmp$gamma,3)
  #tmp$p_est_q0 = round(((tmp$gamma-tmp$shared)/2)/(tmp$gamma-(tmp$gamma-tmp$shared)/2),3)
  tmp = colMeans(tmp)
  tmp
  # N_hat[2]
  
  # likelihood of shared = 7, gamma = 9
  p = 0.047
  q = seq(0,1,length=100)
  plot(q, log((p*q)^2+
              (1-p)*(p*(1-q))+
                (p*(1-q))*(p*(1-q))), type = "l")
  
  
  
  # aggregate by site
  divdat_ag = with(divdat, aggregate(
    list(p_hat=p_hat, pzero_same=pzero_same, pzero_diff=pzero_diff),
    list(SITE_CODE=SITE_CODE, plot=plot),
    FUN = sum
  ))
  
  # estimate overall q (for each site?)
  
  plot(divdat$cover, divdat$p_hat)
  plot(divdat$cover, regression_same)
  
  # get q for non-self-resurvey
  plot(regression_diff, regression_same)
  
  matplot(divdat$cover, cbind(regression_same, regression_diff), pch=1)
  
  
  # aggregate by site
  tmp = divdat
  R = predict(mod_pd, newdata = tmp, re_formula = ~0)[,1]
  p = tmp$p_hat
  q = pmin(pmax((-(p^2*(R-2))+2*p*(R-2)+R)/(2*p*(R-1)),0),1)
  plot(p,q)
  
  p = p_hat; R = regression_diff
  q_hat = (-(p^2*(R-2))+2*p*(R-2)+R)/(2*p*(R-1))
  
  R = regression_same/regression_diff
  p = p_hat
  q_hat = (-(p^2*(R - 1)) + p*(R - 2) + 2*R - 1)/
    (p*(R - 2) + R)
  
  # if error, then 50% chance of ID mistake, and 50% of missing observation
  # or just chance of missing?
}


# compare observations (both alpha change and turnover)


# loop: test observed change across deltaS values
deltaS_lvls = -10:10

# two samples: one for the initial, and one for the changed
# then compare them

## TODO: or just once? or different for cover vs. id error?
