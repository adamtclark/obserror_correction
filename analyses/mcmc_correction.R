setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
#source("util/save_model_outputs.R")
require(brms)
require(nlme)
load("output/brms_coverid_models_small.rda")
#load("output/diversity_change.rda")

set.seed(123432)
source("analyses/diversity_change_functions.R")

collst = c("purple", "forestgreen", "dodgerblue3", "firebrick")
alpha_level = 0.5
cex_level = 0.8
lwduse = 1.5

# parameters
a = 0.05 # intercept
b = 1.2  # scaling exponent
pimm = 0.1 # probability of immigration

### now try MCMC fitting example
set.seed(23432)

## First, fix alpha diversity estimates
# sample data
true_div1 = true_div0 = divdat[divdat$SITE_CODE=="kzbg.at.dragnet" & divdat$plot==8,]
#true_div1 = simulate_change(true_div0, deltaS = 3, process_noise = FALSE)
#true_div1 = simulate_change(true_div1, deltaS = -5, process_noise = FALSE)

# remove species
p_lose = pnorm(0, true_div0$cover, sqrt(a*true_div0$cover^b))
cover_new = rnorm(nrow(true_div0), true_div0$cover, sqrt(a*true_div0$cover^b))
cover_new[cover_new<=0] = 0
cover_new[cover_new>0.01] = round(cover_new[cover_new>0.01],2)
cover_new[cover_new>0 & cover_new <= 0.01] = 0.01
true_div1$cover = cover_new

# add species
new_species = divdat[divdat$SITE_CODE=="kzbg.at.dragnet" & divdat$plot!=8,]
new_species = new_species[!new_species$species%in%true_div1$species,]
new_species = aggregate(x = list(cover = new_species$cover),
                        by = list(species = new_species$species, group = new_species$group),
                        FUN = sum)
N = sum(true_div0$cover>0)
S = sum(new_species$cover>0)+N

n_immigrants = rbinom(1,S-N,pimm)
imm_ps = sample(1:nrow(new_species),n_immigrants)

true_div1 = rbind(true_div1,
                  data.frame(SITE_CODE = unique(true_div1$SITE_CODE),
                             plot = unique(true_div1$plot),
                             species = new_species$species[imm_ps],
                             cover = new_species$cover[imm_ps],
                             group = new_species$group[imm_ps]))

# get immigration likelihood
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
true_cover_diff = (obs_div_merged0$cover+obs_div_merged1$cover)/((obs_div_merged0$cover>0) + (obs_div_merged1$cover>0))
#tmpps = true_cover_diff0==0 & obs_div_merged0$cover_resurvey_diff != 0
#true_cover_diff0[tmpps] = obs_div_merged0$cover_resurvey_diff[tmpps]
obs_div_same0 = get_pq_est(obs_div_merged0, type = "same", true_covers = obs_div_merged0$cover)
obs_div_diff0 = get_pq_est(obs_div_merged0, type = "diff", true_covers = true_cover_diff)
obs_div_same1 = get_pq_est(obs_div_merged1, type = "same", true_covers = obs_div_merged1$cover)
obs_div_diff1 = get_pq_est(obs_div_merged1, type = "diff", true_covers = obs_div_merged1$cover)

## Next, fix beta diversity estimates
#abs(pmin(deltaS,0))/sum(!is.na(obs_div_same0$cover) & obs_div_same0$cover>0)

# losses
sps = 12 # 3 is obs error, 6 is lost
obs_div_diff0[sps,]
pi = obs_div_diff0$pi[sps]
qi = obs_div_diff0$qi[sps]
obs_div_diff0$cover[sps]
(1-pi)*p_lose # 1 found and 2 was lost
(1-pi)*pi # 1 found and 2 was missed or misided

(1-pi)*pi*qi # 1 found and 2 missed
(1-pi)*pi*(1-qi) # 1 found and 2 mis-ided

# mis-ids







#################################################
# Let's start easy - just use real data
p_lose = pnorm(0, true_cover_diff, sqrt(a*true_cover_diff^b))

transition_type = data.frame(
  both_found = obs_div_diff0$cover>0 & obs_div_diff1$cover>0,
  extinction = obs_div_diff0$cover>0 & obs_div_diff1$cover==0,
  immigration = obs_div_diff0$cover==0 & obs_div_diff1$cover>0
  )
#transition_type$extinction[10] = TRUE
#transition_type$both_found[10] = FALSE
sum(2*log(1-obs_div_diff0$pi[transition_type$both_found])) + # both found
      sum(log(1-obs_div_diff0$pi[transition_type$extinction]) + log(p_lose[transition_type$extinction])) + #1 found and 2 lots
      sum(log(1-obs_div_diff0$pi[transition_type$immigration])+log(dbinom(sum(transition_type$immigration),S-N,pimm))) # immigration and found

############### Now let's make it harder - observed data
# create estimated true community state
true_cover_est = (obs_div_merged0$cover_resurvey_diff+obs_div_merged1$cover_resurvey_diff)/((obs_div_merged0$cover_resurvey_diff>0) + (obs_div_merged1$cover_resurvey_diff>0))

obs_div_diff0 = get_pq_est(obs_div_merged0, type = "diff", true_covers = true_cover_est)
obs_div_diff1 = get_pq_est(obs_div_merged1, type = "diff", true_covers = true_cover_est)

true_pa_est = data.frame(species = obs_div_diff0$species,
                      obs0 = 0,
                      obs1 = 0)
true_pa_est$obs0[obs_div_diff0$cover>0] = 1
true_pa_est$obs1[obs_div_diff1$cover>0] = 1

obs_cover = data.frame(obs0 = obs_div_merged0$cover_resurvey_diff, obs1 = obs_div_merged1$cover_resurvey_diff)
p_lose = pnorm(0, true_cover_est, sqrt(a*true_cover_est^b))
true_pa_est[true_pa_est$obs0!=1 | true_pa_est$obs1!=1,]


# calculate LL
LL = 0

## colonizations and extinctions from "true" states
# extinction
LL = LL + sum(log(p_lose[true_pa_est$obs0==1 & true_pa_est$obs1==0]))
LL = LL + sum(log(1-p_lose[true_pa_est$obs0==1 & true_pa_est$obs1==1]))

# colonisation
ncol = sum(true_pa_est$obs0==0 & true_pa_est$obs1==1)
LL = LL + pbinom(ncol, S-N, pimm, lower.tail = FALSE, log.p = TRUE)

# 1 obs - 1 true observations
ps = true_pa_est$obs0==1 & obs_cover$obs0>0
LL = LL + sum(log((1-obs_div_diff0$pi[ps])))

ps = true_pa_est$obs1==1 & obs_cover$obs1>0
LL = LL + sum(log((1-obs_div_diff1$pi[ps])))


########### TODO - working from here
# 1 obs - 0 true observations (misID)
ps_misID = true_pa_est$obs0==0 & obs_cover$obs0>0
true_pa_est[ps,]

# 0 obs - 1 true observations (missed)
ps_miss = true_pa_est$obs0==1 & obs_cover$obs0==0

# for every misID, choose whether to associate to a missed
## WITHIN obs0/true0
# would need to link an obs > 0 / true == 0 TO an obs == 0 / true > 0
# NOTE that obs1/true1 has already been accounted for above.
# TO DO THIS: misID = p*(1-q)
## IF misID0 and miss1, THEN p*(1-q)*p*(q) (across both lines)
## IF misID0 and correct1, THEN ignore additional p*q on case with obs0/true1
##### OH! Basically, just determines whether we multiply by p(1-q) or by pq
## TO DO: think
# in case where obs1 misses the misided species, then no entry for real species exists (an no observation)
# do we need to add in a line for this?...

data.frame(true_pa_est, obs_cover)[ps_misID | ps_miss,]


# 0-0





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





###### old

if(FALSE) {
  #case with deltaS = 0
  p_lose = pnorm(0, obs_div_diff0$cover, sqrt(a*obs_div_diff0$cover^b))
  
  # misid's are more likely than extinctions?
  ps = which(obs_div_diff0$cover_resurvey_diff==0 & obs_div_diff0$cover != 0)
  obs_div_diff0$pi[ps] # missed or misided
  p_lose[ps]
  
  # misid's vs. immigration?
  ps = which(obs_div_diff0$cover_resurvey_diff!=0 & obs_div_diff0$cover == 0)
  obs_div_diff0[ps,]$pi*obs_div_diff0[ps,]$qi
  
  obs_div_diff0$cover_resurvey_diff[ps]
  #hist(quantile(td, probs = runif(1e3), type = 4))
  
  # simulate change
  divdat_small = divdat[divdat$SITE_CODE=="kzbg.at.dragnet",]
  divdat_small_new = simulate_change(divdat_small, deltaS = 3)
  divdat_small_new = simulate_change(divdat_small_new, deltaS = -5)
  
  divdat_0 = divdat_small[divdat_small$plot==8,]
  divdat_1 = divdat_small_new[divdat_small_new$plot==8,]
  
  divdat_0_obs = simulate_noise(divdat_0)
  divdat_1_obs = simulate_noise(divdat_1)
  
  # extract true change as array
  true_change = data.frame(
    species = divdat_0$species,
    group = divdat_0$group,
    cover_0 = divdat_0$cover,
    cover_1 = NA
  )
  
  true_change$cover_1 = divdat_1$cover[match(divdat_0$species, divdat_1$species)]
  
  ps = which(!(divdat_1$species %in% divdat_0$species))
  true_change=rbind(true_change,
                    data.frame(
                      species = divdat_1$species[ps],
                      group = divdat_1$group[ps],
                      cover_0 = 0,
                      cover_1 = divdat_1$cover[ps]
                    ))
  
  true_change$cover_1[is.na(true_change$cover_1)] = 0
  #true_change$cover_1 = round(true_change$cover_1,2)
  true_change = true_change[order(true_change$species),]
  
  
  # extract obs change as array
  obs_change = data.frame(
    species = divdat_0_obs$species,
    group = divdat_0_obs$group,
    cover_0_true = divdat_0_obs$cover,
    cover_1_true = NA,
    cover_0_obs = divdat_0_obs$cover_resurvey_diff,
    cover_1_obs = NA
  )
  
  obs_change$cover_1_obs = divdat_1_obs$cover_resurvey_diff[match(divdat_0_obs$species, divdat_1_obs$species)]
  obs_change$cover_1_true = divdat_1_obs$cover[match(divdat_0_obs$species, divdat_1_obs$species)]
  
  ps = which(!(divdat_1_obs$species %in% divdat_0_obs$species))
  obs_change=rbind(obs_change,
                   data.frame(
                     species = divdat_1_obs$species[ps],
                     group = divdat_1_obs$group[ps],
                     cover_0_true = 0,
                     cover_1_true = divdat_1_obs$cover[ps],
                     cover_0_obs = 0,
                     cover_1_obs = divdat_1_obs$cover_resurvey_diff[ps]
                   ))
  
  obs_change$cover_1_true[is.na(obs_change$cover_1_true)] = 0
  obs_change$cover_1_obs[is.na(obs_change$cover_1_obs)] = 0
  #obs_change$cover_1 = round(obs_change$cover_1,2)
  obs_change = obs_change[order(obs_change$species),]
  
  # true alpha diversity change
  alpha_true = c(sum(obs_change$cover_0_true>0),
                 sum(obs_change$cover_1_true>0))
  
  # true beta diversity
  gamma_true = sum(obs_change$cover_0_true+obs_change$cover_1_true>0)
  beta_true = gamma_true/mean(alpha_true)
  
  # round to 0.01 accuracy
  obs_change[,3:6][obs_change[,3:6]>0 & obs_change[,3:6] < 0.01] = 0.01
  obs_change[,3:6] = round(obs_change[,3:6],2)
  
  # make matrix for storing error types
  obs_change
  
  # Next step...
  #obs_change$cover_0_correct = obs_change$cover_0
  #obs_change$cover_0_misid = 
  #obs_change$cover_0_missed = 
}





