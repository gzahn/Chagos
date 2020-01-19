# Analyses of Chagos Islands coral-associated bacteria 16S Amplicons
# Data cleaning and preparation - building a phylogeny
# Author: Geoffrey Zahn

# Packages and functions ####
library(tidyverse)
library(phyloseq)
library(vegan)
theme_set(theme_bw())

# Load phyloseq object from 01_Process_Raw_Data.R ####
ps <- readRDS("./output/phyloseq_object_16S_noncontam.RDS")

# Remove non-bacteria/archaea ####
ps <- subset_taxa(ps, Kingdom %in% c("Bacteria","Archaea"))

# Inspect and clean taxonomy
ps <- subset_taxa(ps, Phylum != "NA")

# Tidy up metadata ... proper date format, temp and salinity ranges, etc #### 

      # List sample variables
      glimpse(ps@sam_data)

# Dates   
ps@sam_data$date <- as.POSIXct(ps@sam_data$date,format = '%m/%d/%Y')

# Temp
ps@sam_data$`Temp Range (oC)`[is.na(ps@sam_data$`Temp Range (oC)`)] <- "NA-NA"
ps@sam_data$`Temp Range (oC)`[ps@sam_data$`Temp Range (oC)`=="Blank"] <- "Blank-Blank"
min <- as.numeric(unlist(map(str_split(ps@sam_data$`Temp Range (oC)`,"-"),1)))
max <- as.numeric(unlist(map(str_split(ps@sam_data$`Temp Range (oC)`,"-"),2)))
ps@sam_data$TempRangeMin <- min
ps@sam_data$TempRangeMax <- max

# Salinity
ps@sam_data$`Salinity Range (ppt)`[is.na(ps@sam_data$`Salinity Range (ppt)`)] <- "NA-NA"
ps@sam_data$`Salinity Range (ppt)`[ps@sam_data$`Salinity Range (ppt)`=="Blank"] <- "Blank-Blank"
min <- as.numeric(unlist(map(str_split(ps@sam_data$`Salinity Range (ppt)`,"-"),1)))
max <- as.numeric(unlist(map(str_split(ps@sam_data$`Salinity Range (ppt)`,"-"),2)))
ps@sam_data$SalinityRangeMin <- min
ps@sam_data$SalinityRangeMax <- max

# Depth
ps@sam_data$`Depth Range (m)` <- str_remove(ps@sam_data$`Depth Range (m)`,pattern = "m")
ps@sam_data$`Depth Range (m)`[is.na(ps@sam_data$`Depth Range (m)`)] <- "NA-NA"
ps@sam_data$`Depth Range (m)`[ps@sam_data$`Depth Range (m)`=="Blank"] <- "Blank-Blank"
min <- as.numeric(unlist(map(str_split(ps@sam_data$`Depth Range (m)`,"-"),1)))
max <- as.numeric(unlist(map(str_split(ps@sam_data$`Depth Range (m)`,"-"),2)))
ps@sam_data$DepthRangeMin <- min
ps@sam_data$DepthRangeMax <- max


# Remove Blank Controls ####
ps <- subset_samples(ps, Island != "Blank")


# Clean up Sample Data names ####
names(sample_data(ps)) <- 
  c("LibraryID","Barcode","SampleID","Date","Island","Site","Lat","Lon","ReefType","Exposure",
  "SpeciesTentative","SpeciesConfirmed","ColonyColor","AvgSiteTemp","TempRange_C","Salinity",
  "SalinityRange_ppt","Depth_m","DepthRange_m","TempRangeMin","TempRangeMax","SalinityRangeMin",
  "SalinityRangeMax","DepthRangeMin","DepthRangeMax")

# Convert Sample Data classes ####
glimpse(sample_data(ps))

ps@sam_data$LibraryID <- as.character(ps@sam_data$LibraryID)
ps@sam_data$Lat <- as.numeric(ps@sam_data$Lat)
ps@sam_data$Lon <- as.numeric(ps@sam_data$Lon)
ps@sam_data$AvgSiteTemp <- as.numeric(ps@sam_data$AvgSiteTemp)
ps@sam_data$Salinity <- as.numeric(ps@sam_data$Salinity)
ps@sam_data$Depth_m <- as.numeric(ps@sam_data$Depth_m)

# Change ColonyColor to non-ordered factor (Healthy Coral as reference category)

ps@sam_data$ColonyColor <- plyr::mapvalues(ps@sam_data$ColonyColor,from = unique(ps@sam_data$ColonyColor),
                to=c("Healthy","Pale","Bleached","Very pale"))

ps@sam_data$ColonyColor <- factor(ps@sam_data$ColonyColor,
                                   levels = c("Healthy","Pale","Very pale","Bleached"))

# create 4 groups of avg temperature
clus <- kmeans(na.omit(ps@sam_data$AvgSiteTemp), 3)
clus$cluster[which(!is.na(ps@sam_data$AvgSiteTemp))] <- clus$cluster
ps@sam_data$AvgSiteTempGroup <- clus$cluster

plot(ps@sam_data$AvgSiteTemp,ps@sam_data$AvgSiteTempGroup)

ps@sam_data$AvgSiteTempGroup <- plyr::mapvalues(ps@sam_data$AvgSiteTempGroup,
                                                from=as.character(unique(ps@sam_data$AvgSiteTempGroup)),
                                                to=c("<30.5","30.5-31",">31"))

ps@sam_data$AvgSiteTempGroup[which(is.na(ps@sam_data$AvgSiteTemp))] <- NA

# convert to non-ordered factor
ps@sam_data$AvgSiteTempGroup <- factor(ps@sam_data$AvgSiteTempGroup,
                                       levels = c("<30.5","30.5-31",">31"))


# Remove unneeded columns
ps@sam_data[,grep("Range",names(ps@sam_data))] <- NULL
ps@sam_data$Barcode <- NULL
ps@sam_data$SpeciesTentative <- NULL

ps@sam_data 
skimr::skim(ps@sam_data)

unique(ps@tax_table[,2])

# Clean up non-bacterial taxa ####
bact <- subset_taxa(ps, Kingdom == "Bacteria")
bact <- subset_taxa(bact, Phylum != "Cyanobacteria")
bact <- subset_taxa(bact, Family != "Mitochondria")
ps <- bact

ntaxa(ps)

# Save cleaned phyloseq object ####
saveRDS(ps, "./output/cleaned_ps_object.RDS")






# Site
# ReefType
# Exposure
# Species of Coral
# Mixture of above?

