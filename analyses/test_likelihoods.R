#### WORKS!
# TODO:
# come up with change algorithm that gives good scenarios of change



setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
require(brms)
require(nlme)
load("output/brms_coverid_models_small.rda")
source("analyses/diversity_change_functions.R")

# plotting parameters
collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.5
cex_level = 0.8
lwduse = 1.5

# model parameters
a = 0.05 # intercept
b = 1.2  # scaling exponent
pimm = 0.1 # probability of immigration

# likelihood function
get_LL = function() {
  # calculate LL
  LL = 0
  
  ## species observed, but estimated state absent (misID)
  ps_misID0 = true_pa_est$obs0==0 & obs_cover$obs0>0
  LL = LL + sum(log(obs_div_diff0$pi[ps_misID0]*(1-obs_div_diff0$qi[ps_misID0])))
  
  ps_misID1 = true_pa_est$obs1==0 & obs_cover$obs1>0
  LL = LL + sum(log(obs_div_diff1$pi[ps_misID1]*(1-obs_div_diff1$qi[ps_misID1])))
  
  # species not observed, but estimated state present (missed)
  # time0
  ps_miss0 = true_pa_est$obs0==1 & obs_cover$obs0==0
  # ignore misID species (as misses are pseudo-absences)
  if(sum(ps_misID0)>sum(ps_miss0)) {
    # if not enough missed species to account of all misIDs, add one
    delta_misID = sum(ps_misID0)-sum(ps_miss0)
    ps_tmp = sample(which(obs_cover$obs0==0 & true_pa_est$obs0==0 &!ps_misID0), delta_misID)
    true_pa_est$obs0[ps_tmp] = 1
    
    ps_miss0 = true_pa_est$obs0==1 & obs_cover$obs0==0
  }
  misID_match0 = sample(which(ps_miss0),sum(ps_misID0))
  ps_miss0[misID_match0] = FALSE
  LL = LL + sum(log(obs_div_diff0$pi[ps_miss0]*obs_div_diff0$qi[ps_miss0]))
  
  # time1
  ps_miss1 = true_pa_est$obs1==1 & obs_cover$obs1==0
  # ignore misID species (as misses are pseudo-absences)
  if(sum(ps_misID1)>sum(ps_miss1)) {
    # if not enough missed species to account of all misIDs, add one
    delta_misID = sum(ps_misID1)-sum(ps_miss1)
    ps_tmp = sample(which(obs_cover$obs1==0 & true_pa_est$obs1==0 &!ps_misID1), delta_misID)
    true_pa_est$obs1[ps_tmp] = 1
    
    ps_miss1 = true_pa_est$obs1==1 & obs_cover$obs1==0
  }
  misID_match1 = sample(which(ps_miss1),sum(ps_misID1))
  ps_miss1[misID_match1] = FALSE
  LL = LL + sum(log(obs_div_diff1$pi[ps_miss1]*obs_div_diff1$qi[ps_miss1]))
  
  ## colonizations and extinctions based on estimated states
  # extinction
  LL = LL + sum(log(p_lose[true_pa_est$obs0==1 & true_pa_est$obs1==0]))
  LL = LL + sum(log(1-p_lose[true_pa_est$obs0==1 & true_pa_est$obs1==1]))
  
  # colonisation
  ncol = sum(true_pa_est$obs0==0 & true_pa_est$obs1==1)
  LL = LL + pbinom(ncol, S-N, pimm, lower.tail = FALSE, log.p = TRUE)
  
  ## observation and estimated state both present
  ps = true_pa_est$obs0==1 & obs_cover$obs0>0
  LL = LL + sum(log((1-obs_div_diff0$pi[ps])))
  
  ps = true_pa_est$obs1==1 & obs_cover$obs1>0
  LL = LL + sum(log((1-obs_div_diff1$pi[ps])))
  
  return(list(LL=LL, true_pa_est=true_pa_est))
}


set.seed(234324)

## sample data
site = "kzbg.at.dragnet"
plot = 8

# div0 is initial community
# div1 is community post-extinction and colonisation
true_div1 = true_div0 = divdat[divdat$SITE_CODE==site & divdat$plot==plot,]

## simulate extinction
p_lose = pnorm(0, true_div0$cover, sqrt(a*true_div0$cover^b))
true_div1$cover = pmax(0,rnorm(nrow(true_div0), true_div0$cover, sqrt(a*true_div0$cover^b)))

## simulate colonization from rest of site
new_species = divdat[divdat$SITE_CODE==site & divdat$plot!=plot,]
new_species = new_species[!new_species$species%in%true_div1$species,]
new_species = aggregate(x = list(cover = new_species$cover),
                        by = list(species = new_species$species, group = new_species$group),
                        FUN = sum)
N = sum(true_div0$cover>0)
S = sum(new_species$cover>0)+N

n_immigrants = rbinom(1,S-N,pimm)

if(n_immigrants>0) {
  imm_ps = sample(1:nrow(new_species),n_immigrants)
  true_div1 = rbind(true_div1,
                    data.frame(SITE_CODE = site,
                               plot = plot,
                               species = new_species$species[imm_ps],
                               cover = new_species$cover[imm_ps],
                               group = new_species$group[imm_ps]))
}

## simulate sampling error
obs_div0 = simulate_noise(divdat = true_div0)
obs_div1 = simulate_noise(divdat = true_div1)

## div0 and div1 
obs_div = unique(rbind(obs_div0[,c("SITE_CODE", "plot", "species", "group")],
                       obs_div1[,c("SITE_CODE", "plot", "species", "group")]))
obs_div = obs_div[order(obs_div$species),]
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

# add in missing species from region with zero cover
ps = which(!new_species$species%in%obs_div$species)
obs_div_merged0 = rbind(obs_div_merged0,
                        data.frame(SITE_CODE = site,
                                   plot = plot,
                                   species = new_species$species[ps],
                                   group = new_species$group[ps],
                                   cover = new_species$cover[ps],
                                   cover_resurvey_same = 0,
                                   cover_resurvey_diff=0))

obs_div_merged1 = rbind(obs_div_merged1,
                        data.frame(SITE_CODE = site,
                                   plot = plot,
                                   species = new_species$species[ps],
                                   group = new_species$group[ps],
                                   cover = new_species$cover[ps],
                                   cover_resurvey_same = 0,
                                   cover_resurvey_diff=0))

## get p and q (based on observed cover in each plot and time period)
cover_est0 = obs_div_merged0$cover_resurvey_diff
cover_est0[cover_est0==0] = obs_div_merged1$cover_resurvey_diff[cover_est0==0]
cover_est0[cover_est0==0] = obs_div_merged0$cover[cover_est0==0]

cover_est1 = obs_div_merged1$cover_resurvey_diff
cover_est1[cover_est1==0] = obs_div_merged0$cover_resurvey_diff[cover_est1==0]
cover_est1[cover_est1==0] = obs_div_merged1$cover[cover_est1==0]

obs_div_diff0 = get_pq_est(obs_div_merged0, type = "diff", true_covers = cover_est0)
obs_div_diff1 = get_pq_est(obs_div_merged1, type = "diff", true_covers = cover_est1)


## extract observed and true covers
obs_cover = data.frame(species = obs_div_diff0$species, obs0 = obs_div_diff0$cover_resurvey_diff, obs1 = obs_div_diff1$cover_resurvey_diff)
true_cover = data.frame(species = obs_div_diff0$species,
                        true0 = true_div0$cover[match(obs_div_diff0$species,true_div0$species)],
                        true1 = true_div1$cover[match(obs_div_diff0$species,true_div1$species)])
true_cover[is.na(true_cover)] = 0

## estimate extinction probability based on observed cover
p_lose = pnorm(0, cover_est0, sqrt(a*cover_est0^b))


### get likelihood vs. change
# initialize true presence absence
true_pa_est = data.frame(species = obs_div_diff0$species,
                         obs0 = 0,
                         obs1 = 0)
true_pa_est$obs0[true_cover$true0>0] = 1
true_pa_est$obs1[true_cover$true1>0] = 1

# get initial likelihood at true community state
LL0 = get_LL()$LL

# get diversity metrics at true community state
alpha0 = colSums(true_pa_est[,2:3]>0)
gamma0 = sum(rowSums(true_pa_est[,2:3])>0)
delta_richness0 = alpha0[2]-alpha0[1]
beta0 = gamma0/mean(alpha0)

niter = 10 # number of random steps to take
datout = data.frame(
  LL = rep(NA, niter),
  LLrat = NA,
  alpha1 = NA,
  alpha2 = NA,
  gamma = NA,
  beta = NA,
  delta_richness = NA
)


for(i in 1:niter) {
  # choose species to remove or add
  if(i > 1) {
    rowps = sample(1:nrow(true_pa_est),1)
    colps = sample(2:3,1)
    true_pa_est[rowps,colps] = abs(true_pa_est[rowps,colps]-1)
    #ps = sample(which(true_pa_est[,-1]>0),1)
    #true_pa_est[,-1][(ps-1)%%nrow(true_pa_est)+1,floor((ps-1)/nrow(true_pa_est))+1] = 0
  }
  
  LLout = get_LL()
  LL = LLout$LL
  true_pa_est = LLout$true_pa_est
  
  # get changes
  alpha = colSums(true_pa_est[,2:3]>0)
  gamma = sum(rowSums(true_pa_est[,2:3])>0)
  delta_richness = alpha[2]-alpha[1]
  beta = gamma/mean(alpha)
  
  datout$LL[i] = LL
  datout$LLrat[i] = LL-LL0
  datout$alpha1[i] = alpha[1]
  datout$alpha2[i] = alpha[2]
  datout$gamma[i] = gamma
  datout$beta[i] = beta
  datout$delta_richness[i] = delta_richness
  
  if(i/10 == floor(i/10))
    round(print(i/niter),2)
}


plot(datout$delta_richness-delta_richness0, datout$LL-LL0)
abline(v=delta_richness0, lty=2)

plot(datout$beta-beta0, datout$LL-LL0)
abline(v=delta_richness0, lty=2)




###### Plotting
# extinction probability vs. abundance
cov_sq = seq(0.001, 1, by = 0.001)
proc_sd = sqrt(a*cov_sq^b)
prob_ext = pnorm(0, cov_sq, proc_sd)

par(mar=c(5,4,2,2))

plot(cov_sq, prob_ext, type = "l", log = "x",
     xaxs = "i", axes = FALSE,
     lwd = 1.5,
     xlab = "",
     ylab = "")
axis(2)
axis(1, at = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), las = 2)
box()
mtext("Species-level cover", 1, line = 3.5)
mtext("Prob. local extinction", 2, line = 2.5)
abline(h=c(0,1), lty=3)

# colonization probability vs. N and S
# extinction probability vs. abundance
N = sum(true_div0$cover>0)
S = sum(new_species$cover>0)+N

Nsq = seq(0,S, by = 1)
prob_imm = pbinom(0,S-Nsq,pimm,lower.tail = FALSE)

plot(Nsq, prob_imm, type = "s",
     xaxs = "i", axes = FALSE,
     lwd = 1.5,
     xlim=c(0,S+1),
     ylim = c(0,1),
     xlab = "",
     ylab = "")
axis(2)
axis(1, at = seq(0,S+1, by=2), las = 2)
box()
mtext("Plot-level richness", 1, line = 3.2)
mtext("Prob. ≥1 immigration", 2, line = 2.5)
abline(h=c(0,1), lty=3)
abline(v=c(N,S), lty=2)

# note: these are equal
#pbinom(0,3,pimm,lower.tail = FALSE)
#dbinom(1,3,pimm)+dbinom(2,3,pimm)+dbinom(3,3,pimm)

# and... this equals 1
#dbinom(1,3,pimm)+dbinom(2,3,pimm)+dbinom(3,3,pimm)+dbinom(0,3,pimm)


