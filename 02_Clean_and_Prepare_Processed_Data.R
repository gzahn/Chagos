# Packages and functions ####
library(tidyverse)
library(phyloseq)
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

# Remove "Barcode" column
ps@sam_data$Barcode <- NULL

# Save cleaned phyloseq object
saveRDS(ps, "./output/cleaned_ps_object.RDS")


# Merge by various factors ####

# Island
ps_island = merge_samples(ps,"Island")
  # Repair metadata
  ps_island@sam_data$Island <- row.names(ps_island@sam_data)
  # Rel-abund
  ps_island_ra <- transform_sample_counts(ps_island, function(x) x/sum(x))
saveRDS(ps_island_ra,"output/ps_island_ra.RDS")

# Colony Color
ps_color <- merge_samples(ps,"ColonyColor")  
  # Repair metadata
  ps_color@sam_data$ColonyColor <- row.names(ps_color@sam_data)
  # Rel-abund
  ps_color_ra <- transform_sample_counts(ps_color, function(x) x/sum(x))
saveRDS(ps_color_ra,"./output/ps_color_ra.RDS")




# Site
# ReefType
# Exposure
# Species of Coral
# Mixture of above?

