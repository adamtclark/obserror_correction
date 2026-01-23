setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())
require(brms)
require(nlme)
require(viridis)
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
  ## NOTE: this script does not account for cases where
  # the number of misIDed species (i.e. observation > 1, true value = 0)
  # exceeds the number of missed species(i.e. observation = 0, true value > 0)
  # For these cases, an additional "latent" species missed by both surveyors
  # would need to be added.
  # As this additional species would simply further reduce the relative likelihood
  # of an incorrect species assemblage, we do not include it here.
  
  
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
    ps_replace = which(obs_cover$obs0==0 & true_pa_est$obs0==0 &!ps_misID0)
    ps_tmp = sample(ps_replace, pmin(delta_misID, length(ps_replace)))
    # note: if delta_misID >  length(ps_replace), then this suggests
    # a species was missed by both observers - not included in this likelihood calc.
    true_pa_est$obs0[ps_tmp] = 1
    
    ps_miss0 = true_pa_est$obs0==1 & obs_cover$obs0==0
  }
  misID_match0 = sample(which(ps_miss0),pmin(sum(ps_miss0), sum(ps_misID0)))
  ps_miss0[misID_match0] = FALSE
  LL = LL + sum(log(obs_div_diff0$pi[ps_miss0]*obs_div_diff0$qi[ps_miss0]))
  
  # time1
  ps_miss1 = true_pa_est$obs1==1 & obs_cover$obs1==0
  # ignore misID species (as misses are pseudo-absences)
  if(sum(ps_misID1)>sum(ps_miss1)) {
    # if not enough missed species to account of all misIDs, add one
    delta_misID = sum(ps_misID1)-sum(ps_miss1)
    ps_replace = which(obs_cover$obs1==0 & true_pa_est$obs1==0 &!ps_misID1)
    ps_tmp = sample(ps_replace, pmin(delta_misID, length(ps_replace)))
    # note: if delta_misID >  length(ps_replace), then this suggests
    # a species was missed by both observers - not included in this likelihood calc.
    true_pa_est$obs1[ps_tmp] = 1
    
    ps_miss1 = true_pa_est$obs1==1 & obs_cover$obs1==0
  }
  misID_match1 = sample(which(ps_miss1),pmin(sum(ps_miss1), sum(ps_misID1)))
  ps_miss1[misID_match1] = FALSE
  LL = LL + sum(log(obs_div_diff1$pi[ps_miss1]*obs_div_diff1$qi[ps_miss1]))
  
  ## colonizations and extinctions based on estimated states
  # extinction
  LL = LL + sum(log(p_lose[true_pa_est$obs0==1 & true_pa_est$obs1==0]))
  LL = LL + sum(log(1-p_lose[true_pa_est$obs0==1 & true_pa_est$obs1==1]))
  
  # colonisation
  ncol = sum(true_pa_est$obs0==0 & true_pa_est$obs1==1)
  if(S-N>0) {
    LL = LL + pbinom(ncol, S-N, pimm, lower.tail = FALSE, log.p = TRUE)
  }
  
  ## observation and estimated state both present
  ps = true_pa_est$obs0==1 & obs_cover$obs0>0
  LL = LL + sum(log((1-obs_div_diff0$pi[ps])))
  
  ps = true_pa_est$obs1==1 & obs_cover$obs1>0
  LL = LL + sum(log((1-obs_div_diff1$pi[ps])))
  
  return(list(LL=LL, true_pa_est=true_pa_est))
}


if(FALSE) {
  set.seed(234324)
  datouttot = NULL
  plotlst = unique(divdat[,c("SITE_CODE", "plot")])
  for(k in 1:nrow(plotlst)) {
    ## sample data
    site = plotlst$SITE_CODE[k]
    plot = plotlst$plot[k]
    
    # div0 is initial community
    # div1 is community post-extinction and colonization
    true_div1 = true_div0 = divdat[divdat$SITE_CODE==site & divdat$plot==plot,]
    
    ## simulate extinction
    p_lose = pnorm(0, true_div0$cover, sqrt(a*true_div0$cover^b))
    true_div1$cover = pmax(0,rnorm(nrow(true_div0), true_div0$cover, sqrt(a*true_div0$cover^b)))
    
    ## simulate colonization from rest of site
    N = sum(true_div0$cover>0)
    new_species = divdat[divdat$SITE_CODE==site & divdat$plot!=plot,]
    new_species = new_species[!new_species$species%in%true_div1$species,]
    if(nrow(new_species)>0) {
      new_species = aggregate(x = list(cover = new_species$cover),
                              by = list(species = new_species$species, group = new_species$group),
                              FUN = sum)
      S = sum(new_species$cover>0)+N
    } else {
      S = N
    }
    
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
    if(length(ps)>0) {
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
    }
    
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
    true_pa_est0 = true_pa_est
    
    # get initial likelihood at true community state
    LL0 = get_LL()$LL
    
    # get diversity metrics at true community state
    alpha0 = colSums(true_pa_est[,2:3]>0)
    gamma0 = sum(rowSums(true_pa_est[,2:3])>0)
    delta_richness0 = alpha0[2]-alpha0[1]
    beta0 = gamma0/mean(alpha0)
    
    niter0 = 10 # number of random steps to take
    datout = data.frame(
      LL = rep(NA, 2*niter0+1),
      LLrat = NA,
      alpha0 = NA,
      alpha1 = NA,
      gamma = NA,
      beta = NA,
      delta_richness = NA
    )
    
    
    n=1
    for(j in 1:3) {
      true_pa_est = true_pa_est0
      if(j == 1) {
        simtype = "null"
        niter = 1
      }
      if(j == 2) {
        simtype = "loss"
        niter = pmin(niter0, sum(true_pa_est$obs1==1))
      } else if(j==3) {
        simtype = "gain"
        niter = pmin(niter0, sum(true_pa_est$obs1==0))
      }
      if(niter>0) {
        for(i in 1:niter) {
          # choose species to remove or add
          if(simtype == "loss") {
            # remove a species - bias towards survey 1 to cause loss
            rowps = sample(which(true_pa_est$obs1==1),1)
            colps = sample(2:3,1,prob=c(0.3,0.7))
            true_pa_est[rowps,colps] = 0
          } else if(simtype=="gain") {
            # add a species - bias towards survey 1 to cause gain
            rowps = sample(which(true_pa_est$obs1==0),1)
            colps = sample(2:3,1,prob=c(0.3,0.7))
            true_pa_est[rowps,colps] = 1
          }
          
          #rowps = sample(1:nrow(true_pa_est),1)
          #colps = sample(2:3,1)
          #true_pa_est[rowps,colps] = abs(true_pa_est[rowps,colps]-1)
          #ps = sample(which(true_pa_est[,-1]>0),1)
          #true_pa_est[,-1][(ps-1)%%nrow(true_pa_est)+1,floor((ps-1)/nrow(true_pa_est))+1] = 0
          
          LLout = get_LL()
          LL = LLout$LL
          
          # get changes
          alpha = colSums(true_pa_est[,2:3]>0)
          gamma = sum(rowSums(true_pa_est[,2:3])>0)
          delta_richness = alpha[2]-alpha[1]
          beta = gamma/mean(alpha)
          
          datout$LL[n] = LL
          datout$LLrat[n] = LL-LL0
          datout$alpha0[n] = alpha[1]
          datout$alpha1[n] = alpha[2]
          datout$gamma[n] = gamma
          datout$beta[n] = beta
          datout$delta_richness[n] = delta_richness
          
          n = n+1
        }
      }
    }
    datout = datout[!is.na(datout$LL) & is.finite(datout$LL),]
    datout$site = site; datout$plot = plot
    
    datouttot[[k]] = datout
    print(k/nrow(plotlst))
  }
  save(list=c("datouttot"), file = "~/Dropbox/Projects/117_ObservationError/src/output/likelihoods.rda")
} else {
  load("~/Dropbox/Projects/117_ObservationError/src/output/likelihoods.rda")
}

###### Plotting
# aggregate likelihoods across all samples
rngfn = function(x) {
  LL0 = x$LL[1]
  LLrat = x$LL[is.finite(x$LL)]-LL0
  range(LLrat,na.rm=TRUE)
}
tmp = sapply(datouttot, rngfn)
siten = sapply(datouttot, function(x) nrow(x))
minn = 6 # minimum sample size
LLrng = range(tmp[is.finite(tmp)])
collst = viridis(length(datouttot))

## Delta_richness
plot(c(-10,10),LLrng,type="n", xlab = "", ylab = "")
#plot(c(-0.2,1),LLrng,type="n", xlab = "", ylab = "")

xtot=NULL; ytot=NULL; fieldnum = NULL
for(i in 1:length(datouttot)) {
  datout = datouttot[[i]]
  delta_richness0 = datout$delta_richness[1]
  beta0 = datout$beta[1]
  gamma0 = datout$gamma[1]
  alpha0 = c(datout$alpha0[1],datout$alpha1[1])
  LL0 = datout$LL[1]
  
  x = datout$delta_richness-delta_richness0
  #x = datout$beta-beta0
  y = datout$LL-max(datout$LL,na.rm=TRUE)
  ps = which(is.finite(x) & is.finite(y))
  if(length(ps)>=minn) {
    x = x[ps]; y=y[ps]
    xrng = sort(c(0,seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length=100)))
    mod = loess(sqrt(-y)~x, enp.target = floor(sqrt(length(x))))
    prd = -(predict(mod, newdata = data.frame(x=xrng))^2)
    #xrng = sort(x)
    #prd = y[order(x)]
    lines(xrng, prd, lwd = 0.8, col = adjustcolor(collst[i], alpha.f = 0.3))
    points(x,y,pch=16,cex=0.2,col=adjustcolor(collst[i],alpha.f = 0.3))
    xtot=c(x,xtot)
    ytot=c(y,ytot)
    fieldnum=c(fieldnum,rep(i,length(x)))
  }
}

ps = which(is.finite(xtot) & is.finite(ytot))
x = xtot[ps]; y=ytot[ps]
xrng = sort(c(0,seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length=100)))
wts = (table(fieldnum)/length(fieldnum))[fieldnum] # weight loess by field-level sample size
wts[!is.finite(wts)]=NA
wts = wts/sum(wts,na.rm=TRUE)
wts[!is.finite(wts)] = 0
mod = loess(sqrt(-y)~x, enp.target = floor(sqrt(sum(siten>=minn))), weights = wts)
prd = -(predict(mod, newdata = data.frame(x=xrng))^2)
lines(xrng, prd, lwd = 3, col = adjustcolor(1, alpha.f = 0.8))
abline(h=0,v=0,lty=2)

## Beta
plot(c(-0.2,1),LLrng,type="n", xlab = "", ylab = "")

xtot=NULL; ytot=NULL; fieldnum = NULL
for(i in 1:length(datouttot)) {
  datout = datouttot[[i]]
  delta_richness0 = datout$delta_richness[1]
  beta0 = datout$beta[1]
  gamma0 = datout$gamma[1]
  alpha0 = c(datout$alpha0[1],datout$alpha1[1])
  LL0 = datout$LL[1]
  
  x = datout$beta-beta0
  y = datout$LL-max(datout$LL,na.rm=TRUE)
  ps = which(is.finite(x) & is.finite(y))
  if(length(ps)>=minn) {
    x = x[ps]; y=y[ps]
    xrng = sort(c(0,seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length=100)))
    mod = loess(sqrt(-y)~x, enp.target = floor(sqrt(length(x))))
    prd = -(predict(mod, newdata = data.frame(x=xrng))^2)
    #xrng = sort(x)
    #prd = y[order(x)]
    lines(xrng, prd, lwd = 0.8, col = adjustcolor(collst[i], alpha.f = 0.3))
    points(x,y,pch=16,cex=0.2,col=adjustcolor(collst[i],alpha.f = 0.3))
    xtot=c(x,xtot)
    ytot=c(y,ytot)
    fieldnum=c(fieldnum,rep(i,length(x)))
  }
}

ps = which(is.finite(xtot) & is.finite(ytot))
x = xtot[ps]; y=ytot[ps]
xrng = sort(c(0,seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length=100)))
wts = (table(fieldnum)/length(fieldnum))[fieldnum] # weight loess by field-level sample size
wts[!is.finite(wts)]=NA
wts = wts/sum(wts,na.rm=TRUE)
wts[!is.finite(wts)] = 0
mod = loess(sqrt(-y)~x, enp.target = floor(sqrt(sum(siten>=minn))), weights = wts)
prd = -(predict(mod, newdata = data.frame(x=xrng))^2)
lines(xrng, prd, lwd = 3, col = adjustcolor(1, alpha.f = 0.8))
abline(h=0,v=0,lty=2)


## extinction probability vs. abundance
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


