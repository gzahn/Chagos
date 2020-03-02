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

# Custom color palettes
source("./R/palettes.R")
colorblindr::palette_plot(pal.discrete)
colorblindr::palette_plot(pal.continuous)


# Register Google Maps Key
register_google(" ???hidden??? ")


# Load data ####
ps <- readRDS("./output/cleaned_ps_object_w_tree.RDS")
names(sample_data(ps))

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

# Coral Species ##
ps_species <- merge_samples(ps,"SpeciesConfirmed")  
# Repair metadata
ps_species@sam_data$SpeciesConfirmed <- row.names(ps_species@sam_data)
# Rel-abund
ps_species_ra <- transform_sample_counts(ps_species, function(x) x/sum(x))
saveRDS(ps_species_ra,"./output/ps_species_ra.RDS")

# Avg Site Temp Group ##
ps_temp <- merge_samples(ps,"AvgSiteTempGroup")
# Repair metadata
ps_temp@sam_data$AvgSiteTempGroup <- row.names(ps_temp@sam_data)
# Rel-abund
ps_temp_ra <- transform_sample_counts(ps_temp, function(x) x/sum(x))
saveRDS(ps_temp_ra,"./output/ps_temp_ra.RDS")


# BarPlots of relabundance for the merged data ####

# phylum by island
plot_bar2(ps_island_ra, fill = "Phylum") + labs(x="Island",y="Relative abundance") +  # island
  theme_bw() + 
  scale_fill_manual(values = pal.discrete) + 
  theme(axis.text.x = element_text(angle=60,hjust=1,face="bold"),
        axis.title = element_text(face="bold",size=16),
        legend.title = element_text(face="bold",size=16))
ggsave("./output/figs/Barplot_Phylum_relabund_by_Island.png",dpi=300)

# by coral colony color
plot_bar2(ps_color_ra, fill = "Phylum") + labs(x="Coral color",y="Relative abundance") +  # color
  theme_bw() + 
  theme(axis.title = element_text(face="bold",size=16),
        legend.title = element_text(face="bold",size=16)) +
        axis.text.x = element_text(angle=60,hjust=1,face="bold") +
  scale_fill_manual(values = pal.discrete) + theme(axis.text.x = element_text(angle=60,hjust=1))
ggsave("./output/figs/Barplot_Phylum_relabund_by_Coral_Color.png",dpi=300)

# by coral species
plot_bar2(ps_species_ra, fill = "Phylum") + labs(x="Coral species",y="Relative abundance") +  # species
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold.italic"),
        axis.title = element_text(face="bold",size=16),
        legend.title = element_text(face="bold",size=16)) + 
  scale_fill_manual(values = pal.discrete) 
ggsave("./output/figs/Barplot_Phylum_relabund_by_Coral_Species.png",dpi=300)

# by temp group
plot_bar2(ps_temp_ra, fill = "Phylum") + labs(x="Temperature group",y="Relative abundance") +  # temp
  theme_bw() + 
  theme(axis.text.x = element_text(face="bold",size=12),
        axis.title = element_text(face="bold",size=16),
        legend.title = element_text(face="bold",size=16)) +
  scale_fill_manual(values = pal.discrete) 
ggsave("./output/figs/Barplot_Phylum_relabund_by_Temp.png",dpi=300)


# Transform full phyloseq object to relative abundance ####
ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))

# clean taxa names
ps <- clean_taxa_names(ps)

# Look at metadata ####
names(sample_data(ps))
skimr::skim(sample_data(ps_ra))
apply(ps_ra@sam_data,2,table)

# make data frame of sample data #
meta = as(sample_data(ps_ra),"data.frame")

# Map of sites ####
chagosmap <- get_map(location = c(lon=mean(ps_ra@sam_data$Lon),
                                  lat=mean(ps_ra@sam_data$Lat)),
                     source="google", maptype="terrain", crop=FALSE, zoom = 9)

ggmap(chagosmap) + 
  geom_point(data = meta,aes(y=Lat,x=Lon))

ggplot(meta, aes(x=Lon,y=Lat)) + geom_point()



# Plot richness based on ColonyColor ####
plot_richness(ps,x="ColonyColor",measures=c("Shannon","Observed"))

richness = specnumber(otu_table(ps_ra))
shannon = vegan::diversity(otu_table(ps_ra),index = "shannon")
plot(ps_ra@sam_data$AvgSiteTemp,shannon)


# Model basic alpha diversity based on many factors ####
mod.shannon = aov(data=meta, shannon ~ (Island + ColonyColor + AvgSiteTemp)* SpeciesConfirmed)
mod.richness = aov(data=meta, richness ~ (Island + ColonyColor + AvgSiteTemp)* SpeciesConfirmed)
summary(mod.shannon)
summary(mod.richness)

glm.shannon = glm(data=meta, shannon ~ (Island + ColonyColor + AvgSiteTemp)* SpeciesConfirmed)
car::Anova(glm.shannon)
glm.richness = glm(data=meta, richness ~ (Island + ColonyColor + AvgSiteTemp)* SpeciesConfirmed)
car::Anova(glm.richness)


# Stepwise AIC
step = stepAIC(mod.shannon)
summary(step)

step2 = stepAIC(mod.richness)
summary(step)
summary(step2)

# Shannon model: Island and AvgSiteTemp are significant
# Richness model: Island, AvgSiteTemp, AvgSiteTemp:SpeciesConfirmed are significant
#                 ColonyColor is 'not quite significant'


# Make output file of stat tables for alpha diversity ##
sink("./output/stats/alpha_diversity_models.txt")
print("Shannon Diversity Stats Table")
print(formula(glm.shannon))
car::Anova(glm.shannon)
print("Richness Stats Table")
print(formula(glm.richness))
car::Anova(glm.richness)
sink(NULL)

# What is correlation between ColonyColor and AvgSiteTemp
plot(meta$ColonyColor,meta$AvgSiteTemp)
ggplot(meta, aes(x=ColonyColor,y=AvgSiteTemp,fill=ColonyColor)) + 
  geom_boxplot(alpha=.5) + geom_point(alpha=.5) +
  labs(y="Mean Site Temperature",x="Colony Color",fill="Colony Color") +
  scale_fill_manual(values = pal.discrete[c(11,4,6,1)]) + theme_bw() +
  theme(axis.text.x = element_text(face="bold",size=12),
        axis.title = element_text(face="bold",size=16),
        legend.title = element_text(face="bold",size=16))
ggsave("./output/figs/Temp_vs_ColonyColor.png",dpi=300)


# UniFrac distance calculation ####
set.seed(123)
unifrac <- UniFrac(ps_ra, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

# Ordinations ####

# Transform full phyloseq object to relative abundance (again) 
ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))

# Ordinations with various methods #
dca <- ordinate(ps_ra,method="DCA")
nmds <- ordinate(ps_ra,method = "NMDS")
ordu <-  ordinate(ps_ra, "PCoA", "unifrac", weighted=TRUE)

# Plot ordinations
plot_ordination(ps_ra,dca,color = "ColonyColor") + scale_color_manual(values = pal.discrete[c(11,4,6,1)]) +
  theme_bw()
ggsave("./output/figs/DCA_Ordination_by_ColonyColor.png",dpi=300)

plot_ordination(ps_ra,dca,color = "Island") + theme_bw() + scale_color_manual(values=pal.discrete)
ggsave("./output/figs/DCA_Ordination_by_Island.png",dpi=300)

plot_ordination(ps_ra, ordu, color="ColonyColor") + theme_bw() + scale_color_manual(values = pal.discrete[c(11,4,6,1)])
ggsave("./output/figs/PCoA_Ordination_by_ColonyColor.png",dpi=300)

plot_ordination(ps_ra, ordu, color="SpeciesConfirmed") + theme_bw() + scale_color_manual(values = pal.discrete[c(1,3)],name="Coral species")+
  theme(legend.text = element_text(face="italic"))
ggsave("./output/figs/PCoA_Ordination_by_CoralSpecies.png",dpi=300)

# permanova

veganmodel <- vegan::adonis(data = ,formula = otu_table(ps_ra) ~ ps_ra@sam_data$ColonyColor * ps_ra@sam_data$SpeciesConfirmed)

sink("./output/stats/permanova_results.txt")
veganmodel
sink(NULL)


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
             low="#000033", high="#66CCFF") +
  labs(sample.label="Colony color",y="ESVs") +
  theme(axis.text.y = element_blank())

hm_plot$scales$scales[[2]]$name <- "Colony color"
hm_plot$scales$scales[[1]]$name <- "ESVs"
hm_plot$labels$fill <- "Relative abundance"

hm_plot
ggsave("./output/figs/Heatmap_ColonyColor.png",dpi=300)


# Phylogenetic tree plots ####
plot_tree(ps_ra,color="ColonyColor",ladderize = "left",base.spacing = .04) +
  scale_color_manual(values = pal.discrete[c(11,4,6,1)]) + 
  theme(legend.position = "bottom",
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=12,face="bold"))
ggsave("./output/figs/GTR_Phylogenetic_Tree_Colored_by_ColonyColor.png",dpi=300,height = 14,width = 6)

# Network analyses ####
ig <- make_network(ps_ra, max.dist = .6)

# by Coral Species
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "SpeciesConfirmed",label = NULL) + 
  scale_color_manual(values = pal.discrete[c(1,5)]) + labs(color = "Coral Species")
ggsave("./output/figs/Network_Plot_Species.png",dpi=300)

# by Temperature
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "AvgSiteTemp",label = NULL) +
  scale_color_gradient(low=pal.discrete[16],high = pal.discrete[4]) +
  labs(color="Mean Site\nTemperature")
ggsave("./output/figs/Network_Plot_Temperature.png",dpi=300)

# by ColonyColor
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "ColonyColor",label = NULL) +
  scale_color_manual(values = pal.discrete[c(11,4,6,1)]) + labs(color="Colony color")
ggsave("./output/figs/Network_Plot_ColonyColor.png",dpi=300)

# by Island
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "Island",label = NULL) +
  scale_color_manual(values = pal.discrete)
ggsave("./output/figs/Network_Plot_Island.png",dpi=300)

# by Depth_m
set.seed(123)
plot_network(ig, physeq = ps_ra, color = "Depth_m",label = NULL) +
  scale_color_gradient(low="#37b4c4",high="#2a40a1") +
  labs(color="Depth (m)")
ggsave("./output/figs/Network_Plot_Depth.png",dpi=300)








#


















