# Analyses of Chagos Islands coral-associated bacteria 16S Amplicons
# Differential abundance analyses
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

# clean taxa names
ps <- clean_taxa_names(ps)

# Differential Abundance Based on Colony Color ####

# Find differentially-abundant groups...plot them as in http://statisticaldiversitylab.com/blog/167093

names(sample_data(ps))
set.seed(123)
da_analysis <- differentialTest(formula = ~ ColonyColor, #abundance
                                phi.formula = ~ 1, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps,
                                fdr_cutoff = 0.05)

plot(da_analysis)
ggsave("./output/figs/colonycolor_differential_abundance_plot.png",dpi = 300,width = 14,height = 6)

sigs_colcolor <- otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = ps)
length(sigs_colcolor)



# Plot differentially-abundant taxa (Colony Color) ####

# glm models with ColonyColor and no intercept

sigs_colcolor
colorblindr::palette_plot(pal)

set.seed(123)
corncob_da1 <- bbdml(formula = OTU8 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p1 <- plot(corncob_da1,color="ColonyColor") + ggtitle(sigs_colcolor[1]) + 
  scale_color_manual(values = pal[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

set.seed(123)
corncob_da2 <- bbdml(formula = OTU14 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p2 <- plot(corncob_da2,color="ColonyColor") + ggtitle(sigs_colcolor[2]) + 
  scale_color_manual(values = pal[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank()) 

set.seed(123)
corncob_da3 <- bbdml(formula = OTU24 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p3 <- plot(corncob_da3,color="ColonyColor") + ggtitle(sigs_colcolor[3]) + 
  scale_color_manual(values = pal[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

set.seed(123)
corncob_da4 <- bbdml(formula = OTU36 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p4 <- plot(corncob_da4,color="ColonyColor") + ggtitle(sigs_colcolor[4]) + 
  scale_color_manual(values = pal[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

set.seed(123)
corncob_da5 <- bbdml(formula = OTU110 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p5 <- plot(corncob_da5,color="ColonyColor") + ggtitle(sigs_colcolor[5]) + 
  scale_color_manual(values = pal[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

# combine figs with patchwork
combined_plots <- p1/p2/p3/p4/p5
combined_plots

ggsave(combined_plots,filename = "./output/figs/ColonyColor_diff_abund_taxa.png",
       height = 18,width = 12,dpi = 300,device = "png")

# Export stats tables

sink("./output/stats/corncob_ColonyColor_tables.txt")
print(sigs_colcolor[1])
corncob_da1
print(sigs_colcolor[2])
corncob_da2
print(sigs_colcolor[3])
corncob_da3
print(sigs_colcolor[4])
corncob_da4
print(sigs_colcolor[5])
corncob_da5
sink(NULL)


# Plot differentially-abundant taxa by avg site temp ####

# Differential abundance test by temperature group 
set.seed(123)
da_analysis_temp <- differentialTest(formula = ~ AvgSiteTempGroup, #abundance
                                     phi.formula = ~ AvgSiteTempGroup, #dispersion
                                     formula_null = ~ 1, #mean
                                     phi.formula_null = ~ 1,
                                     test = "Wald", boot = FALSE,
                                     data = ps,
                                     fdr_cutoff = 0.05)

plot(da_analysis_temp)
ggsave("./output/figs/tempgroup_differential_abundance_plot.png",width = 16,height = 6)

sigs_tempgroup <- otu_to_taxonomy(OTU = da_analysis_temp$significant_taxa, data = ps)
length(sigs_tempgroup)
length(ps@sam_data$AvgSiteTempGroup)
levels(ps@sam_data$AvgSiteTempGroup)

set.seed(123)
corncob_da7 <- bbdml(formula = OTU1 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p7 <- plot(corncob_da7,color = "AvgSiteTempGroup") + ggtitle(sigs_tempgroup[1]) + 
  scale_color_manual(values = pal[c(11,3,6)]) + 
  theme(axis.text.x = element_blank())

set.seed(123)
corncob_da8 <- bbdml(formula = OTU3 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p8 <- plot(corncob_da8,color="AvgSiteTempGroup") + ggtitle(sigs_tempgroup[2]) + 
  scale_color_manual(values = pal[c(11,3,6)]) + 
  theme(axis.text.x = element_blank())

temp_group_plot <- p7 / p8
temp_group_plot
ggsave(temp_group_plot,filename = "./output/figs/TempGroup_diff_abund_taxa.png",
       dpi=300,width = 12, height = 6,device = "png")



# Relative abundance plots ####

# Top 12 taxa with highest relative abundance
top12 <- get_top_taxa(ps, 11, relative = TRUE, discard_other = FALSE,
                      other_label = "Other")

# Transform top-12 phyloseq object to relative abundance #
top12_ra <- transform_sample_counts(top12, function(x) x/sum(x))

ps_family <- phyloseq::tax_glom(top12_ra, "Family")
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]

# Melt and plot
phyloseq::psmelt(ps_family) %>%
  ggplot(data = ., aes(x = ColonyColor, y = Abundance,color=ColonyColor)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(height = 0, width = .2,alpha=.5) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  theme(legend.position = "bottom", axis.text.x = element_blank()) +
  scale_color_manual(values = pal[c(11,4,6,1)])

ggsave("./output/figs/most_abundant_families_by_colonycolor.png",dpi=300,width = 10,height = 6)

