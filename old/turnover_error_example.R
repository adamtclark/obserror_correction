niter = 100 # number of iterations
psisq = seq(1, 0.2, length = niter)

dataout = data.frame(iter = 1:niter, psi = NA,
                     Ntrue_1 = NA, Ntrue_2 = NA,
                     Nobs_1 = NA, Nobs_2 = NA,
                     Beta_true = NA, Beta_obs = NA)

# add in change (can be positive or negative
# still need turnover from LOSS
# parameters
N=200 # number of species
J=1   # number of repeated observations
mu=0.15 # mean cover across species
phi=1.5   # phi parameter for beta distribution
sigma=0.2   # measurement variability
gamma0=-0.5 # logit intercept (detection prob)
gamma1=4    # logit slope (detection prob)

for(ii in 1:niter) {
  psi=psisq[ii] # true occupation probability

  # generate the true abundance distribution for the global species pool
  y <- matrix(rbeta(N, mu*phi, (1-mu)*phi), 
                  nrow=N, ncol=1)
  
  #### TODO: change abundances too?
  
  
  # generate two true communities, including absences
  y.mat <- array(dim = c(N, 2))
  for(k in 1:2) {
    y.mat[,1] <- rbinom(N, 1, psi)*y
    y.mat[,2] <- rbinom(N, 1, mean(psisq))*y
  }
  
  
  # generate two messy measurements for the community (non-zero cases only)
  u.array <- array(dim = c(N, J, 2))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:2) {
        u.array[i, j, k] <- ifelse(y.mat[i,k] > 0,
                              plogis(rnorm(1, qlogis(y.mat[i,k]), sigma)),
                              0)
      }
    }
  }
  
  # add in false negatives for both measurement
  x.array <- array(dim = c(N, J, 2))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:2) {
        x.array[i, j, k] <- ifelse(y.mat[i,k] > 0,
                              rbinom(1, 1, 
                                     plogis(gamma0 + gamma1*y.mat[i,k])),
                              0)
      }
    }
  }
  
  # combine observation error and false negatives
  xO = u.array*x.array
  
  # get true Beta
  dataout$Beta_true[ii] = sum(rowSums(y.mat)>0)/mean(colSums(y.mat>0))
  
  # get observed Beta
  dataout$Beta_obs[ii] = sum(rowSums(xO[,1,])>0)/mean(colSums(xO[,1,]>0))
  
  dataout$psi[ii] = psi
  dataout[ii,c("Ntrue_1", "Ntrue_2")] = apply(y.mat>0, 2, sum)
  dataout[ii,c("Nobs_1", "Nobs_2")] = c(apply(xO>0, 2:3, sum))
}

# get richness change´
true_change = dataout$Ntrue_2-dataout$Ntrue_1
obs_change = dataout$Nobs_2-dataout$Nobs_1

plot(true_change, obs_change, type = "p",
     xlab = "true richness change",
     ylab = "observed richness change")
abline(a=0, b=1, lty=2)
abline(h=0, v=0, lty=3)

matplot(psisq/mean(psisq), cbind(true_change, obs_change),
        xlab = "diversity difference",
        ylab = "richness change",
        type = "p", pch = 1)
abline(h=0, v=1, lty=3)

# repeat for turnover
plot(dataout$Beta_true, dataout$Beta_obs, type = "p",
     xlab = "true beta",
     ylab = "observed beta")
abline(a=0, b=1, lty=2)
abline(h=0, v=0, lty=3)

matplot(psisq/mean(psisq), cbind(dataout$Beta_true, dataout$Beta_obs),
        xlab = "diversity difference",
        ylab = "beta",
        type = "p", pch = 1)
abline(v=1, lty=3)

# FINALLY: run regression algorithm to try to correct?
