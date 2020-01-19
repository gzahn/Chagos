# Analyses of Chagos Islands coral-associated bacteria 16S Amplicons
# Initial data explorations
# Author: Geoffrey Zahn

# Packages and functions ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(MASS)
library(car)
library(ggmap)
library(corncob)
library(patchwork)
library(fantaxtic)
source("./R/plot_bar2.R")


# Register Google Maps Key
register_google(" ???hidden??? ")

# Custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000",
        "#493829","#816C5B","#A9A18C","#613318","#855723","#B99C6B",
        "#8F3B1B","#D57500","#DBCA69","#404F24","#668D3C","#BDD09F",
        "#4E6172","#83929F","#A3ADB8")

palplot <- colorblindr::palette_plot(pal)

# Load data ####
ps <- readRDS("./output/cleaned_ps_object.RDS")
psm_island_ra <- readRDS("./output/ps_island_ra.RDS")
psm_color_ra <- readRDS("./output/ps_color_ra.RDS")


# Transform full phyloseq object to relative abundance ####
ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))

# clean taxa names
ps <- clean_taxa_names(ps)

# Look at metadata ####
names(sample_data(ps))
skimr::skim(sample_data(ps_ra))
apply(ps_ra@sam_data,2,table)

# make data frame of sample data ####
meta = as(sample_data(ps_ra),"data.frame")

# Map of sites ####
chagosmap <- get_map(location = c(lon=mean(ps_ra@sam_data$Lon),
                                  lat=mean(ps_ra@sam_data$Lat)),
                     source="google", maptype="terrain", crop=FALSE, zoom = 9)

ggmap(chagosmap) + 
  geom_point(data = meta,aes(y=Lat,x=Lon))

ggplot(meta, aes(x=Lon,y=Lat)) + geom_point()


# Initial Plots
plot_bar2(psm_color_ra,fill = "Phylum") + scale_fill_manual(values=pal)

# Plot richness based on ColonyColor ####
plot_richness(ps,x="ColonyColor",measures=c("Shannon","Observed"))

richness = specnumber(otu_table(ps_ra))
shannon = diversity(otu_table(ps_ra),index = "shannon")
plot(ps_ra@sam_data$AvgSiteTemp,shannon)


# Model basic alpha diversity based on all factors ####
mod.shannon = aov(data=meta, shannon ~ (Island + ColonyColor + AvgSiteTemp)* SpeciesConfirmed)
mod.richness = aov(data=meta, richness ~ (Island + ColonyColor + AvgSiteTemp)* SpeciesConfirmed)
summary(mod.shannon)
summary(mod.richness)

# Stepwise AIC
step = stepAIC(mod.shannon)
summary(step)

step2 = stepAIC(mod.richness)
summary(step)
summary(step2)

# Shannon model: Island and AvgSiteTemp are significant
# Richness model: Island, AvgSiteTemp, AvgSiteTemp:SpeciesConfirmed are significant
#                 ColonyColor is 'not quite significant'


# What is correlation between ColonyColor and AvgSiteTemp
plot(meta$ColonyColor,meta$AvgSiteTemp)
ggplot(meta, aes(x=ColonyColor,y=AvgSiteTemp,fill=ColonyColor)) + 
  geom_boxplot(alpha=.5) + geom_point(alpha=.5) +
  labs(y="Mean Site Temperature",x="Colony Color",fill="Colony Color")
ggsave("./output/figs/Temp_vs_ColonyColor.png",dpi=300)


# Ordinations ####

# Transform full phyloseq object to relative abundance (again) 
ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))

dca <- ordinate(ps_ra,method="DCA")
nmds <- ordinate(ps_ra,method = "NMDS")

plot_ordination(ps,dca,color = "ColonyColor") + scale_color_manual(values = pal[c(11,4,6,1)])
ggsave("./output/figs/DCA_Ordination_by_ColonyColor.png",dpi=300)

plot_ordination(ps,dca,color = "Island")
ggsave("./output/figs/DCA_Ordination_by_Island.png",dpi=300)
















#













# x-axis = samples, y-axis = relative abundance of WHATEVER group

# Find potential taxa to plot
plot_bar2(ps_ra,fill="Phylum") + scale_fill_manual(values = pal)

cyanobacteria = subset_taxa(ps_ra, Phylum = "Bacteroidetes")
plot_bar2(cyanobacteria, fill="ColonyColor")










