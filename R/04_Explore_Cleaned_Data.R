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
library(doParallel)
library(foreach)
library(igraph) 
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
ps <- readRDS("./output/cleaned_ps_object_w_tree.RDS")

# Merge by various factors ####

# Island ##
ps_island = merge_samples(ps,"Island")
# Repair metadata
ps_island@sam_data$Island <- row.names(ps_island@sam_data)
# Rel-abund
ps_island_ra <- transform_sample_counts(ps_island, function(x) x/sum(x))
saveRDS(ps_island_ra,"output/ps_island_ra.RDS")

# Colony Color ##
ps_color <- merge_samples(ps,"ColonyColor")  
# Repair metadata
ps_color@sam_data$ColonyColor <- row.names(ps_color@sam_data)
# Rel-abund
ps_color_ra <- transform_sample_counts(ps_color, function(x) x/sum(x))
saveRDS(ps_color_ra,"./output/ps_color_ra.RDS")


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


# Bar Plots, merged by various sample data ####
plot_bar2(ps_color_ra,fill = "Phylum") + scale_fill_manual(values=pal) + labs(y="Relative abundance") +
  theme(axis.text = element_text(face = "bold",size = 12),
        axis.title = element_text(face="bold",size=16),
        legend.title = element_text(face="bold",size=16))
ggsave("./output/figs/Barplot_Phylum_by_ColonyColor.png",dpi=300)

plot_bar2(ps_island_ra,fill = "Phylum") + scale_fill_manual(values=pal)  + labs(y="Relative abundance") +
  theme(axis.text = element_text(face = "bold",size = 12),
        axis.title = element_text(face="bold",size=16),
        legend.title = element_text(face="bold",size=16))
ggsave("./output/figs/Barplot_Phylum_by_Island.png",dpi=300)


# Plot richness based on ColonyColor ####
plot_richness(ps,x="ColonyColor",measures=c("Shannon","Observed"))

richness = specnumber(otu_table(ps_ra))
shannon = diversity(otu_table(ps_ra),index = "shannon")
plot(ps_ra@sam_data$AvgSiteTemp,shannon)


# Model basic alpha diversity based on many factors ####
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
  labs(y="Mean Site Temperature",x="Colony Color",fill="Colony Color") +
  scale_fill_manual(values = pal[c(11,4,6,1)]) + theme_bw()
ggsave("./output/figs/Temp_vs_ColonyColor.png",dpi=300)


# UniFrac distance calculation ####
set.seed(123)
unifrac <- UniFrac(ps_ra, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

# Ordinations ####

# Transform full phyloseq object to relative abundance (again) 
ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))

# Ordinations with various methods ####
dca <- ordinate(ps_ra,method="DCA")
nmds <- ordinate(ps_ra,method = "NMDS")
ordu <-  ordinate(ps_ra, "PCoA", "unifrac", weighted=TRUE)

# Plot ordinations
plot_ordination(ps_ra,dca,color = "ColonyColor") + scale_color_manual(values = pal[c(11,4,6,1)])
ggsave("./output/figs/DCA_Ordination_by_ColonyColor.png",dpi=300)

plot_ordination(ps_ra,dca,color = "Island")
ggsave("./output/figs/DCA_Ordination_by_Island.png",dpi=300)

plot_ordination(ps_ra, ordu, color="ColonyColor")
ggsave("./output/figs/DCA_Ordination_by_Island.png",dpi=300)

# Heatmaps ####
names(sample_data(ps_ra))

# order by ColonyColor
colonycolor_order <- c(grep("Healthy",ps_ra@sam_data$ColonyColor),
grep("Pale",meta$ColonyColor),
grep("Very pale",meta$ColonyColor),
grep("Bleached",meta$ColonyColor))
colonycolor_order <- ps_ra@sam_data$LibraryID[colonycolor_order]

hm_plot <- plot_heatmap(ps_ra, 
             sample.label="ColonyColor", sample.order = colonycolor_order,
             taxa.label = NULL, taxa.order = "Phylum",
             low="#000033", high="#66CCFF",) +
  labs(sample.label="Colony color",y="ESVs") +
  theme(axis.text.y = element_blank())

hm_plot$scales$scales[[2]]$name <- "Colony color"
hm_plot$scales$scales[[1]]$name <- "ESVs"
hm_plot$labels$fill <- "Relative abundance"

hm_plot
ggsave("./output/figs/Heatmap_ColonyColor.png",dpi=300)


# Phylogenetic tree plots ####
plot_tree(ps_ra,color="ColonyColor",ladderize = "left") +
  scale_color_manual(values = pal[c(11,4,6,1)]) + 
  theme(legend.position = "bottom")
ggsave("./output/figs/GTR_Phylogenetic_Tree_Colored_by_ColonyColor.png",dpi=300,height = 8,width = 6)



# Network analyses ####
ig <- make_network(ps_ra, max.dist = .6)

# by Coral Species
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "SpeciesConfirmed",label = NULL) + 
  scale_color_manual(values = pal[c(1,5)]) + labs(color = "Coral Species")
ggsave("./output/figs/Network_Plot_Species.png",dpi=300)

# by Temperature
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "AvgSiteTemp",label = NULL) +
  scale_color_gradient(low=pal[16],high = pal[4]) +
  labs(color="Mean Site\nTemperature")
ggsave("./output/figs/Network_Plot_Temperature.png",dpi=300)

# by ColonyColor
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "ColonyColor",label = NULL) +
  scale_color_manual(values = pal[c(11,4,6,1)]) + labs(color="Colony color")
ggsave("./output/figs/Network_Plot_ColonyColor.png",dpi=300)

# by Island
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "Island",label = NULL) +
  scale_color_manual(values = pal)
ggsave("./output/figs/Network_Plot_Island.png",dpi=300)

# by Depth_m
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "Depth_m",label = NULL) +
  scale_color_gradient(low="#37b4c4",high="#2a40a1") +
  labs(color="Depth (m)")
ggsave("./output/figs/Network_Plot_Depth.png",dpi=300)








#













# x-axis = samples, y-axis = relative abundance of WHATEVER group

# Find potential taxa to plot
plot_bar2(ps_ra,fill="Phylum") + scale_fill_manual(values = pal)

cyanobacteria = subset_taxa(ps_ra, Phylum = "Bacteroidetes")
plot_bar2(cyanobacteria, fill="ColonyColor")










