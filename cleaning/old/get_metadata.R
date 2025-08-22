setwd("~/Dropbox/Projects/075_DRAGNet/AddOns/resurvey/src/")
rm(list = ls())

d = read.csv("../data/2023/processed/NutNet_ReSurvey_nometadat_2023.csv")

nutnet_data = read.csv("~/Dropbox/SharedFolders/NutNet data/comb-by-plot_2024-05-31.csv")

#tmp = unique(nutnet_data[nutnet_data$site_code=="temple.us",c("plot", "block")])
#tmp[order(tmp$plot),]

# merge plot, block, and treatment
d$SITE_CODE[d$SITE_CODE=="msla.us-2"] = "msla.us"
site_index = sort(unique(d$SITE_CODE))

new_cols= c("continent", "country", "region",
            "managed", "burned", "grazed",
            "anthropogenic", "habitat", "elevation",
            "latitude", "longitude", "experiment_type",
            "block", "plot","trt",
            "first_nutrient_year", "first_fenced_year")

for(i in 1:length(new_cols)) {
  d[,new_cols[i]] = NA
}


for(i in 1:length(site_index)) {
  ps = which(nutnet_data$site_code == site_index[i])
  
  if(length(ps>0)) {
    ps2 = which(d$SITE_CODE==site_index[i])
    tmp = unique(nutnet_data[ps,c(3:14,15:16,22,23,24)])
    mtch = match(paste(d$PLOT, d$BLOCK)[ps2], paste(tmp$plot, tmp$block))
    
    d[ps2,new_cols] = tmp[mtch,]
  }
}

# missing info
d$other = NA
dragnet_data = read.csv("new_metadata.csv")

for(i in 1:length(site_index)) {
  ps = which(dragnet_data$site_index == site_index[i])
  
  if(length(ps>0)) {
    ps2 = which(d$SITE_CODE==site_index[i])
    tmp = unique(dragnet_data[ps,])
    
    # all cases with explicit bloc
    d[ps2,c(new_cols, "other")] = tmp[match(d$PLOT[ps2], tmp$plot),c(new_cols, "other")]
    
    
    # all cases without explicit bloc
  }
}

missing = sort(unique(d$SITE_CODE[is.na(d$continent)]))
missing

# set trt to control if no treatment started yet
d$YEAR = as.numeric(substr(d$DATE, 7,10))
d$trt[d$first_nutrient_year > d$YEAR] = "Control"

# remove zero cover
d = d[d$COVER>0,]

# remove disturbance treatments
d$other[is.na(d$other)] = "Control"
d = d[d$other=="Control" | d$SITE_CODE == "stzk.at.dragnet",]

# remove non-species cover (except "Litter" and "Bare ground")
d$SPECIES = gsub("*", "", d$SPECIES, fixed = TRUE)
d$SPECIES = gsub("^", "", d$SPECIES, fixed = TRUE)

sort(unique(d$SPECIES))
d$SPECIES[d$SPECIES%in%c("Achillea millefolium","Achillea_millefolium","ACHILLEA_MILLEFOLIUM","Achillea_millefolium ")] = "Achillea millefolium"
d$SPECIES[d$SPECIES%in%c("Alopecurus pratensis")] = "Alopecurus_pratensis"
d$SPECIES[d$SPECIES%in%c("Ambrosia artemisifolia", "Ambrosia_artemesiifolia", "Ambrosia_artemisiifolia")] = "Ambrosia artemisiifolia"
d$SPECIES[d$SPECIES%in%c("Anthoxanthum_odoratum")] = "Anthoxanthum odoratum"
d$SPECIES[d$SPECIES%in%c("Andropogon_gerardii", "ANDROPOGON_GERARDII", "Andropogon_gerardii")] = "Andropogon gerardii"
d$SPECIES[d$SPECIES%in%c("Alopecurus_pratensis")] = "Alopecurus pratensis"
d$SPECIES[d$SPECIES%in%c("Arrhenatherum_elatius")] = "Arrhenatherum elatius"
d$SPECIES[d$SPECIES%in%c("Bellis perennis", "Bellis_perennis")] = "Bellis perennis"
d$SPECIES[d$SPECIES%in%c("Cynodon dactylon", "Cynodon_dactylon")] = "Cynodon dactylon"
d$SPECIES[d$SPECIES%in%c("Festuca_rubra", "Festuca.rubra.agg.")] = "Festuca rubra"
d$SPECIES[d$SPECIES%in%c("Juncus spp.", "Juncus ")] = "Juncus spp."
d$SPECIES[d$SPECIES%in%c("Lespedeza_capitata", "LESPEDEZA_CAPITATA")] = "Lespedeza capitata"
d$SPECIES[d$SPECIES%in%c("Lotus corniculatus", "Lotus_corniculatus")] = "Lotus corniculatus"
d$SPECIES[d$SPECIES%in%c("Poa alpina", "Poa_alpina")] = "Poa alpina"
d$SPECIES[d$SPECIES%in%c("Poa_pratensis", "POA_PRATENSIS")] = "Poa pratensis"
d$SPECIES[d$SPECIES%in%c("Ranunculus acris", "Ranunculus_acris")] = "Ranunculus acris"
d$SPECIES[d$SPECIES%in%c("Salsola tragus", "Salsola_tragus")] = "Salsola tragus"
d$SPECIES[d$SPECIES%in%c("Schizachrium_scoparium", "SCHIZACHYRIUM_SCOPARIUM", "Schizachyrium_scoparium_scoparium")] = "Schizachrium scoparium"
d$SPECIES[d$SPECIES%in%c("Solanum_carolinense", "SOLANUM_CAROLINENSE")] = "Solanum carolinense"
d$SPECIES[d$SPECIES%in%c("Solidago_speciosa", "SOLIDAGO_SPECIOSA")] = "Solidago speciosa"
d$SPECIES[d$SPECIES%in%c("Taraxacum officinale", "Taraxacum_officinale")] = "Taraxacum officinale"
d$SPECIES[d$SPECIES%in%c("Trifolium pratense", "Trifolium_pratense", "TRIFOLIUM_PRATENSE")] = "Trifolium pratense"
d$SPECIES[d$SPECIES%in%c("Trifolium repens", "Trifolium_repens")] = "Trifolium repens"
d$SPECIES[d$SPECIES%in%c("Moss", "Unknown moss")] = "Moss"

d$SPECIES[d$SPECIES%in%c("Bare ground", "Bare Ground", "Bare_ground", "Bare_Ground", "Bare_soil", "Bare_Soil")] = "Bare ground"

d$SPECIES[d$SPECIES%in%c("Animal", "Animal_disturbance", "Animal_Disturbance", "Disturbance", "Gravel",
                         "Lichen", "Light_cyano_crust", "Overstory", "Rock", "Rocks", "Standing_dead", "Total", "Unknown_red_mushroom",
                         "Woody Overstory", "Woody_overstory", "Woody_Overstory")] = "NON-PLANT"

d$SPECIES = gsub("_", " ", d$SPECIES, fixed = TRUE)
d$SPECIES = tolower(d$SPECIES)

d$SPECIES[d$SPECIES%in%c("gopher mound", "animal disturbance", "basal vegetation", "brown", "logs/sticks", "unknown forb", "vegetation", "total cover", "total plant basal")] = "NON-PLANT"
d$SPECIES[d$SPECIES%in%c("brickellia eupatoriodes")] = "brickellia eupatorioides"
d$SPECIES[d$SPECIES%in%c("bryophytes %", "moss", "dark cyano crust", "lichen %", "unknown moss", "lycopodium.sp.")] = "cryptogam"
d$SPECIES[d$SPECIES%in%c("geranium disectum")] = "geranium dissectum"
d$SPECIES[d$SPECIES%in%c("litterexotic", "litternative")] = "litter"
d$SPECIES[d$SPECIES%in%c("schizachrium scoparium", "schizachyrium scoparius")] = "schizachrium scoparium"
d$SPECIES[d$SPECIES%in%c("sonches asper")] = "sonchus asper"
d$SPECIES[d$SPECIES%in%c("bare ground %", "bare soil", "bare")] = "bare ground"
d$SPECIES[d$SPECIES%in%c("carex spp.")] = "carex sp"
d$SPECIES[d$SPECIES%in%c("erigeron pumilis")] = "erigeron pumilus"
d$SPECIES = gsub("spp.", "sp", d$SPECIES)

d = d[d$SPECIES!="NON-PLANT",]
d = d[d$SPECIES!="non-plant",]

# reshape to multiple columns per survey
tmp = paste(d$SITE, d$SURVEY_ID, d$RESURVEY, d$SPECIES, d$plot)
max(table(tmp)) # should be 1
duplicates = duplicated(tmp)
d = d[!duplicates,]

d$BLOCK = NULL; d$PLOT = NULL; d$SUBPLOT = NULL; d$other = NULL

d_resurvey = NULL
ps = which(site_index=="kzbg.at.dragnet")
site_index = site_index[c(ps, c(1:length(site_index))[-ps])]
for(i in 1:length(site_index)) {
  ps = which(d$SITE_CODE == site_index[i])
  species_list = sort(unique(d$SPECIES[ps]))
  
  base_data = unique(d[ps,c(1:2,9:25)])
  base_data = base_data[order(base_data$plot),]
  
  resurvey_index = sort(unique(d$RESURVEY[ps]))
  
  pltlst = sort(unique(base_data$plot))
  
  
  tmp = base_data[rep(1:nrow(base_data), each = length(species_list)),]
  tmp$species = species_list
  
  tmp = data.frame(tmp, cover_1 = NA, surveyor_1 = NA, date_1 = NA,
                   cover_2 = NA, surveyor_2 = NA, date_2 = NA,
                   cover_3 = NA, surveyor_3 = NA, date_3 = NA,
                   cover_4 = NA, surveyor_4 = NA, date_4 = NA)
  
  for(j in resurvey_index) {
    ps2 = which(d$SITE_CODE == site_index[i] & d$RESURVEY == j)
    
    cover = d$COVER[ps2][match(paste(tmp$species, tmp$plot), paste(d$SPECIES[ps2], d$plot[ps2]))]
    surveyor = ((d$SURVEY_ID[ps2][match(paste(tmp$species, tmp$plot), paste(d$SPECIES[ps2], d$plot[ps2]))]))
    date = ((d$DATE[ps2][match(paste(tmp$species, tmp$plot), paste(d$SPECIES[ps2], d$plot[ps2]))]))
    
    # add in true zeros
    surveyed_plots = d$plot[ps2]
    cover[tmp$plot%in%surveyed_plots & is.na(cover)] = 0
    
    tmp[,paste("cover", j, sep = "_")] = cover
    tmp[,paste("surveyor", j, sep = "_")] = surveyor
    tmp[,paste("date", j, sep = "_")] = date
  }
  
  # remove species that were not found by anyone
  tmp = tmp[rowSums(!is.na(tmp[,grep("cover", colnames(tmp))]))>.1,] # keep only resurveyed plots
  tmp = tmp[rowSums(tmp[,grep("cover", colnames(tmp))], na.rm=TRUE)>0,] # throw out cases where a species was not found
  
  d_resurvey = rbind(d_resurvey, tmp)
}

#write.csv(d_resurvey, "~/Dropbox/Projects/117_ObservationError/src/data/NutNet_ReSurvey_processed_240814.csv", row.names = FALSE)

plot_sites_fun = function(agdat, nameslist = c("div_1", "div_2", "div_3"), alphalevel = 0.2) {
  stlst = sort(unique(agdat$SITE_CODE))
  require(lmodel2)
  
  colnames(agdat)[which(colnames(agdat) %in% nameslist)] = c("div_1", "div_2", "div_3")
  
  slope_lst = matrix(nrow = length(stlst), ncol = 3, data = NA)
  
  suppressWarnings({
    for(i in 1:length(stlst)) {
      ps = which(agdat$SITE_CODE == stlst[i])
      xsq = seq(min(agdat$div_1[ps], na.rm=TRUE), max(agdat$div_1[ps], na.rm=TRUE), length=100)
      collst = viridis(3)
      collst2 = adjustcolor(viridis(3), alpha.f = alphalevel)
      
      if(sum(!is.na(agdat$div_2[ps]+agdat$div_1[ps]))>2) {
        mod1 = suppressMessages(lmodel2(div_2~div_1, agdat[ps,], "relative", "relative"))
        pd1 = mod1$regression.results[4,2]+mod1$regression.results[4,3]*xsq
        #predict(mod1, newdata = data.frame(div_1=xsq))
        lines(xsq, pd1, col = collst2[1])
        slope_lst[i,1] = mod1$regression.results[4,3]
      }
      
      if(sum(!is.na(agdat$div_3[ps]+agdat$div_1[ps]))>2) {
        mod2 = suppressMessages(lmodel2(div_3~div_1, agdat[ps,], "relative", "relative"))
        pd2 = mod2$regression.results[4,2]+mod2$regression.results[4,3]*xsq
        #predict(mod1, newdata = data.frame(div_1=xsq))
        lines(xsq, pd2, col = collst2[2])
        slope_lst[i,2] = mod2$regression.results[4,3]
      }
      if(sum(!is.na(agdat$div_4[ps]+agdat$div_1[ps]))>2) {
        mod3 = suppressMessages(lmodel2(div_4~div_1, agdat[ps,], "relative", "relative"))
        pd3 = mod3$regression.results[4,2]+mod3$regression.results[4,3]*xsq
        #predict(mod1, newdata = data.frame(div_1=xsq))
        lines(xsq, pd3, col = collst2[3])
        slope_lst[i,3] = mod3$regression.results[4,3]
      }
    }})
  return(slope_lst)
}

#pdf("pilot_plots.pdf", width = 4, height = 10)
par(mar=c(4,4,2,2))#, mfrow = c(3,1))

# plot
require(viridis)
collst = viridis(3)
plot(d_resurvey$cover_1, d_resurvey$cover_2, xlab = "cover, first survey", ylab = "cover, repeat surveys", col = collst[1])
points(d_resurvey$cover_1, d_resurvey$cover_3, col = collst[2])
points(d_resurvey$cover_1, d_resurvey$cover_4, col = collst[3])
abline(a=0, b = 1, lty = 2)

slopes = plot_sites_fun(d_resurvey, c("cover_1", "cover_2", "cover_3"))
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")


# estimate sd
#sd_est=2.7
#A1 = rnorm(1e5,0,sd_est); A2 = rnorm(1e5,0,sd_est)
#sqrt(mean((A1-A2)^2)/2)
#sqrt(mean(((A1-A2)^2)/2))
#mean(((A1-A2)^2)/2)
#sd_est^2

mean_cover = rowMeans(d_resurvey[,grep("cover", colnames(d_resurvey))], na.rm=TRUE)
var_cover = apply(d_resurvey[,grep("cover", colnames(d_resurvey))],1,function(x) var(x, na.rm=TRUE))

mean_cover = mean_cover[var_cover>0]
var_cover = var_cover[var_cover>0]

plot(mean_cover, var_cover, log = "xy", xlab = "mean cover est.", ylab = "cover est. variance")
mnsq = exp(seq(log(min(mean_cover, na.rm=TRUE)), log(max(mean_cover, na.rm=TRUE)), length=1e3))
lines(mnsq, mnsq^2, lty = 2, lwd = 2, col = "red")

lines(mnsq, (mnsq*0.1)^2, lty = 3, lwd = 2, col = "red")

mod = loess(log(var_cover)~log(mean_cover))
pred = predict(mod, newdata = data.frame(mean_cover = mnsq))

lines(mnsq, exp(pred), col = "blue", lwd = 2)
legend("topleft", c("sd(error) = mean est.", "sd(error) = 10% of mean est.", "loess fit"), lty = c(2,3,1), lwd = 2, col = c("red", "red", "blue"), bty = "n")





# missing species
prob_zero = apply(d_resurvey[,grep("cover", colnames(d_resurvey))], 1, function(x) mean(x[!is.na(x)]==0))
n_obs = apply(d_resurvey[,grep("cover", colnames(d_resurvey))], 1, function(x) sum(!is.na(x)))
nonzero_mean = apply(d_resurvey[,grep("cover", colnames(d_resurvey))], 1, function(x) mean(x[!is.na(x) & x>0]))
zero_obs = n_obs*prob_zero

plot(nonzero_mean, jitter(prob_zero), xlab = "mean non-zero cover", ylab = "probability of zero cover")

mod = glm(cbind(zero_obs, n_obs)~nonzero_mean, family = "binomial")
summary(mod)

mnsq = (seq((min(nonzero_mean)), max((nonzero_mean)), length=1e3))
pred = predict(mod, newdata = data.frame(nonzero_mean = mnsq), type = "response")

lines(mnsq, pred, col = "blue", lwd = 2)
legend("topright", c("logit regression"), lty = c(1), lwd = 2, col = c("blue"), bty = "n")

#dev.off()



###################
### number of species
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

## richness
d_resurvey_ag = with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
                                div_2 = cover_2,
                                div_3 = cover_3,
                                div_4 = cover_4),
                          by = list(SITE_CODE = SITE_CODE, trt = trt,
                                    first_nutrient_year = first_nutrient_year,
                                    block = block, plot = plot),
                          FUN = function(x) hillfun(x, q = 0)))
agfun = function(q) {
  with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
                     div_2 = cover_2,
                     div_3 = cover_3,
                     div_4 = cover_4),
                by = list(SITE_CODE = SITE_CODE, trt = trt,
                          first_nutrient_year = first_nutrient_year,
                          block = block, plot = plot),
                FUN = function(x) hillfun(x, q = q)))
}

d_resurvey_ag_RI = agfun(0)

collst = viridis(3)
plot(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_2, xlab = "richness, first survey", ylab = "richness, repeat surveys", col = collst[1])
points(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_3, col = collst[2])
points(d_resurvey_ag_RI$div_1, d_resurvey_ag_RI$div_4, col = collst[3])
abline(a=0, b = 1, lty = 2)
slopes = plot_sites_fun(d_resurvey_ag_RI)
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
abline(v = mean(slopes, na.rm=TRUE), lty = 2)


## shannon
d_resurvey_ag_SH = agfun(1)

collst = viridis(3)
plot(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_2, xlab = "Shannon, first survey", ylab = "Shannon, repeat surveys", col = collst[1])
points(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_3, col = collst[2])
points(d_resurvey_ag_SH$div_1, d_resurvey_ag_SH$div_4, col = collst[3])
abline(a=0, b = 1, lty = 2)
slopes = plot_sites_fun(d_resurvey_ag_SH)
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
abline(v = mean(slopes, na.rm=TRUE), lty = 2)

## Simpson
d_resurvey_ag_SI = agfun(2)

collst = viridis(3)
plot(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_2, xlab = "Shannon, first survey", ylab = "Shannon, repeat surveys", col = collst[1])
points(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_3, col = collst[2])
points(d_resurvey_ag_SI$div_1, d_resurvey_ag_SI$div_4, col = collst[3])
abline(a=0, b = 1, lty = 2)
slopes = plot_sites_fun(d_resurvey_ag_SI)
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)
hist(c(slopes), breaks = 20, xlab = "site-level slope", main = "")
abline(v = mean(slopes, na.rm=TRUE), lty = 2)

# note - higher error, due to impact of cover?


# total cover
cover_resurvey_ag = with(d_resurvey[!d_resurvey$species%in%c("bare ground", "litter", "cryptogam"),], aggregate(list(div_1 = cover_1,
                           div_2 = cover_2,
                           div_3 = cover_3,
                           div_4 = cover_4),
                      by = list(SITE_CODE = SITE_CODE, trt = trt,
                                first_nutrient_year = first_nutrient_year,
                                block = block, plot = plot),
                      FUN = function(x) sum(x, na.rm=TRUE)))
cover_resurvey_ag[,6:9][cover_resurvey_ag[,6:9]==0] = NA

collst = viridis(3)
plot(cover_resurvey_ag$div_1, cover_resurvey_ag$div_2, xlab = "total cover, first survey", ylab = "total cover, repeat surveys", col = collst[1])
points(cover_resurvey_ag$div_1, cover_resurvey_ag$div_3, col = collst[2])
points(cover_resurvey_ag$div_1, cover_resurvey_ag$div_4, col = collst[3])
abline(a=0, b = 1, lty = 2)
legend("bottomright", c("survey 2", "survey 3", "survey 4"), pch = 1, col = collst)



# error same vs. different people


# error vs. day lag


# error vs. treatment?



## add in lines per site
# average site-level R2, average cross-site R2


