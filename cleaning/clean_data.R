setwd("~/Dropbox/Projects/117_ObservationError/src/cleaning/")
rm(list = ls())

d = read.csv("../data/cleaning/NutNet_ReSurvey_nometadat_2023.csv")
nutnet_data = read.csv("../data/cleaning/comb-by-plot_2025-08-20.csv")

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
dragnet_data = read.csv("../data/cleaning/new_metadata.csv")

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
d$SPECIES[d$SPECIES%in%c("Alchemilla.sp.")] = "Alchemilla sp"
d$SPECIES[d$SPECIES%in%c("Boechera spp")] = "Boechera sp"
d$SPECIES[d$SPECIES%in%c("Hippocrepis.sp.")] = "Hippocrepis sp"
d$SPECIES[d$SPECIES%in%c("Silene_latifolia_ssp_alba")] = "Silene latifolia"
d$SPECIES[d$SPECIES%in%c("Thymus.sp.")] = "Thymus sp"
d$SPECIES[d$SPECIES%in%c("Vicia")] = "Vicia sp"
d$SPECIES[d$SPECIES%in%c("Echinacea")] = "Echinacea sp"

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
d$SPECIES[d$SPECIES%in%c("sprobolus species")] = "sporobolus sp"

# spelling errors
d$SPECIES[d$SPECIES%in%c("agropyron smithiii")] = "agropyron smithii"
d$SPECIES[d$SPECIES%in%c("allyssum desertorum")] = "alyssum desertorum"
d$SPECIES[d$SPECIES%in%c("andropogon gerardii")] = "andropogon gerardi"
d$SPECIES[d$SPECIES%in%c("axanopus sp.")] = "axonopus sp"
d$SPECIES[d$SPECIES%in%c("boechera spp")] = "boechera sp"
d$SPECIES[d$SPECIES%in%c("carex semprevirens")] = "carex sempervirens"
d$SPECIES[d$SPECIES%in%c("cirsium englemannii")] = "cirsium engelmannii"
d$SPECIES[d$SPECIES%in%c("gallium anisophyllon")] = "galium anisophyllon"
d$SPECIES[d$SPECIES%in%c("hieracium sylvatica")] = "hieracium sylvaticum"
d$SPECIES[d$SPECIES%in%c("melilotus species")] = "melilotus sp"
d$SPECIES[d$SPECIES%in%c("microstis unifolia")] = "microtis unifolia"
d$SPECIES[d$SPECIES%in%c("morella cerifera pumila")] = "myrica cerifera"
d$SPECIES[d$SPECIES%in%c("rubus species")] = "rubus sp"
d$SPECIES[d$SPECIES%in%c("schizachrium scoparium")] = "schizachyrium scoparium"
d$SPECIES[d$SPECIES%in%c("symphyotrichum ericode")] = "symphyotrichum ericoides"
d$SPECIES[d$SPECIES%in%c("symphyotrichym ericoides")] = "symphyotrichum ericoides"
d$SPECIES[d$SPECIES%in%c("symphyotrichym oblongifolium")] = "symphyotrichum oblongifolium"
d$SPECIES[d$SPECIES%in%c("verbascum thaspis")] = "verbascum thapsus"
d$SPECIES[d$SPECIES%in%c("vernona baldwinii")] = "vernonia baldwinii"

d$SPECIES = gsub("spp.", "sp", d$SPECIES)

d = d[d$SPECIES!="NON-PLANT",]
d = d[d$SPECIES!="non-plant",]

d$SPECIES = paste(toupper(substr(d$SPECIES,1,1)), substr(d$SPECIES,2,nchar(d$SPECIES)), sep = "")

# trim trailing space
ps = 1
tmp = d$SPECIES
while(length(ps)>0) {
  ps = which(substr(tmp,nchar(tmp),nchar(tmp))==" ")
  tmp[ps] = substr(tmp,1,nchar(tmp)-1)[ps]
}
d$SPECIES=tmp

# standardise sp.
tmp = strsplit(d$SPECIES, " ")
tmp2 = sapply(tmp, function(x) x[2])
ps = which(!is.na(tmp2) & tmp2 == "sp")
for(i in 1:length(ps)) {
  tmp[[ps[i]]][2] = "sp."
}
tmp = sapply(tmp, function(x) paste(x, collapse = " "))
d$SPECIES=tmp

# remove non-plants
d = d[!d$SPECIES%in%c("Bare ground", "Litter"),]

# remove non-identified plants
d = d[!d$SPECIES%in%c("Cryptogam", "Unknown brassicaceae 2021"),]

# look at species list
splst = sort(unique(d$SPECIES))
splst

# get families
require(U.Taxonstand)
load("../data/cleaning/Plants_WCVP.rdata")
spp_list  <- nameClean(splst, author = F)$NameClean
all_ <- nameMatch(spp_list, database)
all_1     <- all_[match(spp_list, all_$Submitted_Name), ]
all_spp   <- all_1$Accepted_SPNAME
fams      <- all_1$Family

sp_families = data.frame(species = splst, corrected = all_1$Accepted_SPNAME, family = fams)
sp_families$family[sp_families$species=="Axanopus sp."] = "Poaceae"
sp_families$corrected[sp_families$species=="Axanopus sp."] = "Axanopus sp."
which(is.na(sp_families)) # check na's
splst[all_1$Fuzzy] # check spelling errors

write.csv(sp_families, "../data/cleaning/family_lookup.csv", row.names = FALSE)

# update with WCVP taxonomy
d$SPECIES=sp_families$corrected[match(d$SPECIES, sp_families$species)]

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

# match family and corrected species name
ps = match(d_resurvey$species, sp_families$corrected)
d_resurvey$family = sp_families$family[ps]
sort(unique(d_resurvey$family))

d_resurvey$group = rep("Non-Graminoid", nrow(d_resurvey))
d_resurvey$group[d_resurvey$family %in% c("Poaceae", "Cyperaceae", "Juncaceae")] = "Graminoid"

d_resurvey = d_resurvey[,c(1:20, 33, 34, 21:32)]

write.csv(d_resurvey, "../data/NutNet_ReSurvey_processed_250821.csv", row.names = FALSE)
