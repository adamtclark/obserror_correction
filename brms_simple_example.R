require(brms)
require(BayesianTools)
rm(list=ls())
set.seed(12322)

# parameters
N=200 # number of species
J=2   # number of repeated observations
mu=0.15 # mean cover across species
phi=1.5   # phi parameter for beta distribution
psi=0.8 # true occupation probability
sigma=0.2   # measurement variability
gamma0=-0.5 # logit intercept (detection prob)
gamma1=4    # logit slope (detection prob)

#generate the cover values for every iteration
y.mat <- matrix(rbinom(N, 1, psi)*rbeta(N, mu*phi, (1-mu)*phi), 
                nrow=N, ncol=1)

par(mfrow=c(1,1), mar=c(4,4,2,2))
hist(y.mat, breaks = 20)
xT = y.mat[,1]

### Step 1: fit model to true data
moddat = data.frame(rep = 1:nrow(y.mat), xT = xT)

#fit <- brm(
#  formula = bf(
#    xT ~ 1,
#    phi ~ 1,
#    zi ~ 1),
#  family = zero_inflated_beta(), data = moddat,
#  cores = 4, iter = 4000)
#summary(fit)
#c(qlogis(mu), log(phi), qlogis(1-psi))
# fits correctly

## likelihood function
nparam = 3
param = c(mu, phi, psi)
likelihood <- function(param){
  # parameters
  xT_mu = param[1]
  xT_phi = param[2]
  xT_zi = param[3]
  
  # likelihoods
  LL_beta = dbeta(xT[xT>0], xT_mu*xT_phi, (1-xT_mu)*xT_phi, log = TRUE)
  LL_logit = dbinom(xT>0, 1, xT_zi, log = TRUE)
  
  LL_tot = sum(LL_beta)+sum(LL_logit)
  
  return(LL_tot)
}
setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(0, 0, 0),
                             upper = c(1, 2, 1))

## run MCMC and save outputs
nmcmc = 5e4
settings <- list(iterations = nmcmc, consoleUpdates=1e1)
#res <- runMCMC(bayesianSetup = setup, settings = settings)
#smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
#apply(smpout[,1:3],2,function(x) quantile(x, c(0.025, 0.5, 0.957)))
#c(mu, phi, psi)

### step 2: add measurement error to non-zero cases
u.array <- array(dim = c(N, J))
for(i in 1:N){
  for(j in 1:J){
    u.array[i, j] <- ifelse(y.mat[i] > 0,
                            plogis(rnorm(1, qlogis(y.mat[i]), sigma)),
                            0)
  }
}

## likelihood function
xO = u.array
xmu = rowMeans(xO)
ps = xT>0

nparam = 4
param = c(mu, phi, psi, sigma)

likelihood <- function(param){
  # parameters
  xT_mu = param[1]
  xT_phi = param[2]
  xT_zi = param[3]
  xT_sigma = param[4]
  
  # likelihoods
  delta_x = (qlogis(xO[ps,1])-qlogis(xO[ps,2]))
  LL_covererror = dnorm(delta_x, 0, xT_sigma*sqrt(2), log = TRUE)
  LL_beta = dbeta(plogis(rowMeans(qlogis(xO[ps,]))), xT_mu*xT_phi, (1-xT_mu)*xT_phi, log = TRUE)
  #LL_beta = dbeta(xO[ps,], xT_mu*xT_phi, (1-xT_mu)*xT_phi, log = TRUE)
  LL_logit = dbinom(xmu>0, 1, xT_zi, log = TRUE)
  
  LL_tot = sum(LL_covererror) + sum(LL_beta) + sum(LL_logit)
  
  return(LL_tot)
}
setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(0, 0, 0, 0),
                             upper = c(1, 2, 1, 1))

## run MCMC and save outputs
nmcmc = 5e4
settings <- list(iterations = nmcmc, consoleUpdates=1e1)
#res <- runMCMC(bayesianSetup = setup, settings = settings)
#smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
#apply(smpout[,1:4],2,function(x) quantile(x, c(0.025, 0.5, 0.957)))
#c(mu, phi, psi, sigma)

#par(mar=c(4,4,2,2))
#plot(res, start = (nmcmc/3)/2, whichParameters = 1:nparam)


### add in false negatives
x.array <- array(dim = c(N, J))
for(i in 1:N){
  for(j in 1:J){
    x.array[i, j] <- ifelse(y.mat[i] > 0,
                            rbinom(1, 1, 
                                   plogis(gamma0 + gamma1*y.mat[i])),
                            0)
  }
}

cbind(x.array, y.mat>0)

## fit model
xT = y.mat[,1]
xO = u.array*x.array
xmu = rowMeans(xO)
ps = xT>0
Z = xT>0
y = xT

nparam = 6
param = c(mu, phi, psi, sigma, gamma0, gamma1)#,
#xT)

likelihood <- function(param){
  ## parameters
  xT_mu = param[1]   # mean (nonzero) cover across species
  xT_phi = param[2]  # phi for cover across species
  xT_zi = param[3]   # probability of true zero
  xT_sigma = param[4]  # cover observation error
  xT_gamma0 = param[5] # false absence intercept
  xT_gamma1 = param[6] # false absence slope
  
  # helping variables
  #y = param[-c(1:nparam)] # true abundance
  #Z = y>0 # true presence/absence
  X = xO > 0 # detection
  
  ## likelihoods
  # X|Z: observed presence given true presence
  #p_XgZy = plogis(xT_gamma0+xT_gamma1*y)
  #LL_XgZy = matrix(nrow = length(y), ncol = 2)
  #LL_XgZy[Z,] = dbinom(X[Z,], 1, p_XgZy[Z], log = TRUE)
  
  # Z
  #LL_Z = dbinom(Z*1,1,xT_zi, log = TRUE)
  
  LL_PA = matrix(nrow = length(y), ncol = J)
  p_XgZy = plogis(xT_gamma0+xT_gamma1*y)
  for(j in 1:J) {
    # observed presence (exists and was found)
    LL_PA[X[,j],j] = log(xT_zi) + log(p_XgZy[X[,j]])
    
    # observed absence (exists but was not found)
    LL_PA[!X[,j] & Z,j] = log(xT_zi) + log(1-p_XgZy[!X[,j] & Z])
    
    # observed absence (does not exist)
    LL_PA[!X[,j] & !Z,j] = log(1-xT_zi)
    
    #LL_PA[!X[,j],j] = log(exp(log(xT_zi) + log(1-p_XgZy[!X[,j]])) + (1-xT_zi))
  }
  
  # U|X: observed cover given observed presence
  #LL_UgXy = matrix(nrow = length(y), ncol = J)
  #for(j in 1:J) {
  #  # observed presence
  #  LL_UgXy[X[,j],j] = dnorm(qlogis(xO[X[,j],j]), qlogis(y[X[,j]]), xT_sigma, log = TRUE)
  #}
  ps = X[,1] & X[,2]
  delta_x = (qlogis(xO[ps,1])-qlogis(xO[ps,2]))
  LL_UgXy = dnorm(delta_x, 0, xT_sigma*sqrt(2), log = TRUE)
  
  # Y|Z
  LL_ygZ = rep(NA, length(y))
  LL_ygZ[Z] = dbeta(y[Z], xT_mu*xT_phi, (1-xT_mu)*xT_phi, log = TRUE)
  
  # total
  LL_tot = sum(LL_UgXy, na.rm=TRUE)+sum(LL_ygZ, na.rm=TRUE)+
  sum(LL_PA, na.rm=TRUE)
  #sum(LL_XgZy, na.rm=T)+sum(LL_Z, na.rm=TRUE)
  
  return(LL_tot)
}

setup <- createBayesianSetup(likelihood = likelihood,
                             lower = c(0, 0, 0, 0, -3, -3),
                             upper = c(1, 2, 1, 1, 3, 10))

## run MCMC and save outputs
nmcmc = 5e4
settings <- list(iterations = nmcmc, consoleUpdates=1e1)
res <- runMCMC(bayesianSetup = setup, settings = settings)
smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
apply(smpout[,1:6],2,function(x) quantile(x, c(0.025, 0.5, 0.957)))
c(mu, phi, psi, sigma, gamma0, gamma1)

par(mfrow=c(3,2), mar=c(4,4,2,2))
hist(smpout[,1], breaks = 20, main = "mu", xlim = c(0,1)); abline(v = mu, lty=2, col = 2, lwd = 2)
hist(smpout[,2], breaks = 20, main = "phi", xlim = c(0,2)); abline(v = phi, lty=2, col = 2, lwd = 2)
hist(smpout[,3], breaks = 20, main = "psi", xlim = c(0,1)); abline(v = psi, lty=2, col = 2, lwd = 2)
hist(smpout[,4], breaks = 20, main = "sigma", xlim = c(0,1)); abline(v = sigma, lty=2, col = 2, lwd = 2)
hist(smpout[,5], breaks = 20, main = "gamma0", xlim = c(-3,3)); abline(v = gamma0, lty=2, col = 2, lwd = 2)
hist(smpout[,6], breaks = 20, main = "gamma1", xlim = c(-3,10)); abline(v = gamma1, lty=2, col = 2, lwd = 2)

################## pick up from here - think about how to remove y from above?









minv_use = c(0, 0, 0, 0, -3, -3)
maxv_use = c(1, 2, 1, 1, 3, 8)
density_fun <- function(param, minv = minv_use, maxv = maxv_use) {
  # flat priors
  dsum <- 0
  dwidth <- maxv - minv
  dsum <- 0
  for (i in 1:length(minv)) {
    if (param[i] <= maxv[i] & param[i] >= minv[i]) {
      dsum <- dsum + log(1/dwidth[i])
    }
    else {
      dsum <- (-Inf)
    }
  }
  
  #beta prior
  #xT_mu = param[1]   # mean (nonzero) cover across species
  #xT_phi = param[2]  # phi for cover across species
  #dsum = dsum+sum(dbeta(param[-c(1:length(minv))], param[1]*param[2], (1-param[1])*param[2], log = TRUE))
  
  return(dsum)
}

sampler_fun <- function(n = 1, minv = minv_use, maxv = maxv_use) {
  dout = matrix(nrow = n, ncol = length(minv))#+N)
  for(i in 1:length(minv)) {
    dout[,i] = runif(n, minv[i], maxv[i])
  }
  #dout[,-c(1:length(minv))] = rbeta(n*N, mu*phi, (1-mu)*phi)
  return(dout)
}

prior <- createPrior(density = density_fun,
                     sampler = sampler_fun,
                     lower = c(minv_use),#, rep(0, N)),
                     upper = c(maxv_use))#, rep(1, N)))

setup <- createBayesianSetup(likelihood = likelihood,
                             prior = prior)
#                             lower = c(0, 0, 0, 0, -3, -3, rep(0, N)),
#                             upper = c(1, 2, 1, 1, 3, 8, rep(1, N)),
#                             best = c(mu, phi, psi, sigma, gamma0, gamma1, pmax(rowMeans(xO), 0.01)))

nmcmc = 5e4
settings <- list(iterations = nmcmc, consoleUpdates=1e1)

res <- runMCMC(bayesianSetup = setup, settings = settings)
plot(res, whichParameters = 1:6, start = (nmcmc/3)/2)
smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
apply(smpout[,1:nparam],2,function(x) quantile(x, c(0.025, 0.5, 0.957)))
c(mu, phi, psi, sigma, gamma0, gamma1)
## move Y into prior?





bayesianSetup_detfun0 <- createBayesianSetup(likelihood = likelihood_detfun0,
                                             prior = prior_USE)










## try brms on observation process
xO = u.array
xmu = rowMeans(u.array)
moddat = data.frame(rep = 1:nrow(y.mat), xO = c(xO[,1], xO[,2]), xmu = xmu)

fit <- brm(
  formula = bf(
    xO ~ (1|rep),
    phi ~ 1,
    zi ~ 1),
  family = zero_inflated_beta(), data = moddat,
  cores = 4, iter = 4000)
summary(fit)
c(sigma, qlogis(mu), log(phi), qlogis(1-psi))
# clearly not fitting well
