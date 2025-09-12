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
covdat[!is.na(covdat) & covdat==0] = NA
covdat = covdat/100
d$cover = covdat[,1] # keep only first survey

divdat = d[,c("SITE_CODE", "plot", "species", "cover")]
divdat = divdat[!is.na(divdat$cover),]

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
  
  
  # if error, then 50% chance of ID mistake, and 50% of missing observation
  # or just chance of missing?
}


# compare observations (both alpha change and turnover)


# loop: test observed change across deltaS values
deltaS_lvls = -10:10

# two samples: one for the initial, and one for the changed
# then compare them

## TODO: or just once? or different for cover vs. id error?
