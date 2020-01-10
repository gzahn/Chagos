# Packages and functions ####
library(tidyverse)
library(phyloseq)
source("./R/plot_bar2.R")

# Custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")
colorblindr::palette_plot(pal)

# Load data ####
ps <- readRDS("./output/cleaned_ps_object.RDS")
psm_island_ra <- readRDS("./output/ps_island_ra.RDS")
psm_color_ra <- readRDS("./output/ps_color_ra.RDS")

# Look at metadata ####
names(sample_data(ps))
glimpse(sample_data(ps))

# Transform full phyloseq object to relative abundance ####
ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))

plot_bar2(psm_color_ra,fill = "Phylum") #+ scale_fill_manual(values=pal)

# Plot richness
plot_richness(ps,x="ColonyColor",measures=c("Shannon","Observed"))




plot_bar(ps_ra)



