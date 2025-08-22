require(brms)
#require(BayesianTools)
#require(truncnorm)
#require(fitdistrplus)
rm(list=ls())
  
doplot = TRUE

source("~/Dropbox/Rfunctions/logit_funs.R")
#source("~/Dropbox/Rfunctions/truncnorm_funs.R")


### next steps
# try brms implementation
# try to add in probability of zeros

### notes?
# Forecast richness at a scale, and then estimate expected error effects and true changes. 
# Maybe add in prior from Emma? 
# Assumption: SAR, plus that abundance distribution is scale invariant 
# Add in distance matrix option (taxonomic, phylo, etc.)
# Use abundance distribution to estimate gains and losses

### make data
set.seed(23432)

# parameters
N=200 # number of species
J=2   # number of repeated observations
iter=1
mu=0.15 # mean cover across species
phi=1.5   # phi parameter for beta distribution
psi=0.8 # true occupation probability
gamma0=-1.5 # logit intercept (detection error)
gamma1=2    # logit slope (detection error)
sigma=0.5   # measurement variability

#generate the cover values for every iteration
y.mat <- matrix(rbinom(N*iter, 1, psi)*rbeta(N*iter, mu*phi, (1-mu)*phi), 
                nrow=N, ncol=iter)
hist(y.mat, breaks = 20)

## try out brms
moddat = data.frame(rep = 1:nrow(y.mat), xT = y.mat[,1])

fit <- brm(
  formula = bf(
    xT ~ 1,
    phi ~ 1,
    zi ~ 1),
  family = zero_inflated_beta(), data = moddat,
  cores = 4, iter = 4000)
summary(fit)
c(logit(mu), log(phi), logit(1-psi))


#make the measurement matrix for every iteration
u.array <- array(dim = c(N, J, iter))
for(k in 1:iter){
  for(i in 1:N){
    for(j in 1:J){
      u.array[i, j, k] <- ifelse(y.mat[i, k] > 0,
                                 plogis(rnorm(1, qlogis(y.mat[i, k]), sigma)),
                                 0)
    }
  }
}

xO = u.array[,,1]
xmu = rowMeans(u.array)

## try brms on observation process
moddat = data.frame(rep = 1:nrow(y.mat), xO = c(xO[,1], xO[,2]), xmu = xmu)

fit <- brm(
  formula = bf(xO ~ (1|rep),
               phi ~ (1|rep),
               zi ~ 1, family = zero_inflated_beta())+
    bf(logit(xO) ~ (1|rep), family = gaussian()),
  data = moddat,
  cores = 4, iter = 4000)
summary(fit)
c(sigma, logit(mu), log(phi), logit(1-psi))


# pick up from here?





#make a detection history matrix for every iteration
x.array <- array(dim = c(N, J, iter))
for(k in 1:iter){
  for(i in 1:N){
    for(j in 1:J){
      x.array[i, j, k] <- ifelse(y.mat[i, k] > 0,
                                 rbinom(1, 1, 
                                        plogis(gamma0 + gamma1*y.mat[i, k])),
                                 0)
    }
  }
}




















## parameters
b0 = 0.02
b1 = 0.3 # observation error parameters
n = 50  # sample size

## hyperparamters
xT_mu = 0.3
xT_phi = 1.5
a = xT_mu*xT_phi; b = (1-xT_mu)*xT_phi
#a = mu*phi and b = (1−mu)*phi
#mu = a/phi
#phi = b/(1-mu)
#mu = a/b/((1+a/b)) = a/(a+b)
#phi = b/(1-a/(a+b))

b0 = -2
b1 = 3

## true states
xT = rbeta(n, shape1 = a, shape2 = b)
x_sd = (b0 + b1*xT) # observation error sd

#x = rbeta(n, a, b)
#hist(x, breaks = 20)
#dat = data.frame(x=x)
#mout = brm(x~1, family = Beta(), data = dat,
#           cores = 4, iter = 4000)
#summary(mout)
#c(log(mu), phi)

## observed states
xO = ilogit(cbind(rnorm(n, logit(xT), x_sd),
           rnorm(n, logit(xT), x_sd)))
xmu = rowMeans(xO)

if(doplot) {
  #dev.new()
  par(mfrow=c(3,1), mar=c(4,4,2,2))
  hist(xT, breaks = 20)
  plot((xT), x_sd); abline(a=b0, b=b1)
  matplot(xT, xO)
}

moddat = data.frame(rep = 1:nrow(xO), xO1 = xO[,1], xO2 = xO[,2], xOmu = xmu)
medat = data.frame(rep = 1:nrow(moddat), obs = c(moddat$xO1, moddat$xO2), xmu = moddat$xOmu)

#prior_use = prior(normal(-2, 2), nlpar = "b0") + prior(normal(-2, 2), nlpar = "b1")
brm_out_0 <- brm(bf(obs ~ xT, nl = TRUE) +
                   lf(xT ~ (1|rep)) +
                   nlf(sigma ~ log(exp(b0) + exp(b1)*xT)) +
                   lf(b0 ~ 1) +
                   lf(b1 ~ 1),
                 data = medat, family = Beta(),
                 cores = 4, iter = 4000,
                 #save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.85))












### fit model
nparam = 4 # number of parameters
#param = c(b0, b1, xT_mu, xT_sd, xT)
dxest = xT/xmu-1
param = c(b0, b1, xT_mu, xT_sd, dxest)

## likelihood function
likelihood <- function(param){
  # estimated parameters
  b0 = param[1]
  b1 = param[2]
  
  # estimated hyperparameters
  xT_mu = param[3]#xT_mu#mean(xmu)#
  xT_sd = param[4]#xT_sd#sd(xmu)#
  
  # estimated states
  xE = xmu*(1+param[(nparam+1):length(param)])
  x_sd_est = b0+xE*b1
  
  LL1 = sum(log(truncnorm(x = xO[,1], a=0, b=1, mu = xE, sigma = x_sd_est, type = "density")$density))
  LL2 = sum(log(truncnorm(x = xO[,2], a=0, b=1, mu = xE, sigma = x_sd_est, type = "density")$density))
  
  llObservation =  LL1+LL2
  llRandomeffects = sum(log(truncnorm(x = xE, a=0, b=1, mu = xT_mu, sigma = xT_sd, type = "density")$density))

  return(llObservation+llRandomeffects)
}

# try out likelihoods
if(doplot) {
  par(mfrow=c(2,1))
  b0lst = seq(0, 1, length=1e3)
  out = numeric(100)
  for(i in 1:length(b0lst)) {
    out[i] = likelihood(param = c(b0lst[i], b1, xT_mu, xT_sd, dxest))
  }
  plot(b0lst, out, ylim=quantile(out, c(0.5, 1)), type = "l"); abline(v = b0, lty =2)
  abline(v = b0lst[which.max(out)], col = 2, lty = 2)
  
  b1lst = seq(0, 1, length=1e3)
  out = numeric(100)
  for(i in 1:length(b1lst)) {
    out[i] = likelihood(param = c(b0, b1lst[i], xT_mu, xT_sd, dxest))
  }
  plot(b1lst, out, ylim=quantile(out, c(0.5, 1)), type = "l"); abline(v = b1, lty =2)
  abline(v = b1lst[which.max(out)], col = 2, lty = 2)
}

## set up Bayesian run
#c(b0, b1, xT_mu, xT_sd, xT)
setup <- createBayesianSetup(likelihood = likelihood,
                             #lower = c(0, 0), upper = c(1,2))
                             lower = c(0, 0, 0, 0, rep(-1, n)),
                             upper = c(1, 1, 1, 1, rep(1, n)))

## run MCMC and save outputs
nmcmc = 5e4
settings <- list(iterations = nmcmc, consoleUpdates=1e1)
res <- runMCMC(bayesianSetup = setup, settings = settings)

smpout = getSample(res, start = (nmcmc/3)/2, parametersOnly=FALSE)
xE_mu = colMeans(smpout[,1:nparam])
xE_sd = apply(smpout[,1:nparam], 2, sd)
gelman_out = gelmanDiagnostics(res)
gelman_out

## plot diagnostics for a single MCMC run
plot(res, start = (nmcmc/3)/2, whichParameters = 1:nparam)
round(colMeans(smpout[,1:nparam]),3)
c(b0, b1, xT_mu, xT_sd)

par(mfrow=c(4,1), mar=c(4,4,2,2))
hist(smpout[,1], xlim = c(0,1), breaks = 20); abline(v = b0, lty = 2, col = 2, lwd = 2)
hist(smpout[,2], xlim = c(0,1), breaks = 20); abline(v = b1, lty = 2, col = 2, lwd = 2)
hist(smpout[,3], xlim = c(0,1), breaks = 20); abline(v = xT_mu, lty = 2, col = 2, lwd = 2)
hist(smpout[,4], xlim = c(0,1), breaks = 20); abline(v = xT_sd, lty = 2, col = 2, lwd = 2)

par(mfrow=c(1,1))
xTest = rowMeans(xmu*(t(smpout[,(nparam+1):(nparam+n)])+1))
xTest_sd = apply(xmu*(t(smpout[,(nparam+1):(nparam+n)])+1),1,function(x) quantile(x, pnorm(c(-1,1))))
plot(xT, xTest, ylim = range(c(xTest_sd[1,], xTest_sd[2,]))); abline(a = 0, b = 1, lty = 2)
segments(xT, xTest_sd[1,], xT, xTest_sd[2,])

cor(xT, xmu)
cor(xT, xTest)

# plot estimate for observation i
i = sample(n, 1)
hist(xmu[i]*(1+smpout[,(nparam+1):(nparam+n)][,i]), breaks = 20, main = "xT[i]", xlim = c(0,1)); abline(v=xT[i], col = 2, lwd=2)


################ Repeat same analysis in brms

require(brms)
moddat = data.frame(rep = 1:nrow(xO), xO1 = xO[,1], xO2 = xO[,2], xOmu = xmu)
medat = data.frame(rep = 1:nrow(moddat), obs = c(moddat$xO1, moddat$xO2), xmu = moddat$xOmu)

prior_use = prior(normal(-2, 2), nlpar = "b0") + prior(normal(-2, 2), nlpar = "b1")
brm_out_0 <- brm(bf(obs ~ xT, nl = TRUE) +
                   lf(xT | trunc(lb = 0, ub = 1) ~ (1|rep)) +
                   nlf(sigma ~ log(exp(b0) + exp(b1)*xT)) +
                   lf(b0 ~ 1) +
                   lf(b1 ~ 1),
                 data = medat, family = gaussian(), prior = prior_use,
                 cores = 4, iter = 4000,
                 #save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.85))

brm_out_0 <- brm(bf(obs ~ xT, nl = TRUE) +
                   lf(xT | trunc(lb = 0, ub = 1) ~ -1+(xmu+1|rep)) +
                   nlf(sigma ~ log(exp(b0) + exp(b1)*xT)) +
                   lf(b0 ~ 1) +
                   lf(b1 ~ 1),
                 data = medat, family = gaussian(), prior = prior_use,
                 cores = 4, iter = 4000,
                 #save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.85))

medat2 = medat
medat2$xT = xT
brm_out_0 <- brm(bf(obs ~ xT, nl = TRUE) +
                   nlf(sigma ~ log(exp(b0) + exp(b1)*xT)) +
                   lf(b0 ~ 1) +
                   lf(b1 ~ 1),
                 data = medat2, family = gaussian(), prior = prior_use,
                 cores = 4, iter = 4000,
                 #save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.85))


summary(brm_out_0)
exp(fixef(brm_out_0)[,c(1, 3:4)])
c(b0, b1)

fixef(brm_out_0)[1,c(1, 3:4)]
truncnorm(a=0, b=1, mu = xT_mu, sigma = xT_sd, type = "mean")$mean

plot(brm_out_0)

fit_out = fitted(brm_out_0)
#plot(fit_out[1:50,1], fit_out[51:100,1])
xTest = fit_out[1:50,1]
plot(xT, xTest)
points(xT, xmu, col = 2)
cor(xT, xmu)
cor(xT, xTest)

