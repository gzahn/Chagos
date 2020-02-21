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
library(dada2)
library(Biostrings)
source("./R/plot_bar2.R")

# Custom color palettes
source("./R/palettes.R")
colorblindr::palette_plot(pal.discrete)
colorblindr::palette_plot(pal.continuous)

# Load data ####
ps <- readRDS("./output/cleaned_ps_object_w_tree.RDS")
ps@refseq

# Extract sequences and add to phlyloseq object ####
seqs <- rownames(tax_table(ps))

if(identical(seqs,rownames(tax_table(ps)))){
  ps@refseq <- DNAStringSet(x=seqs)
}

# clean taxa names
ps <- clean_taxa_names(ps)

# Differential Abundance Based on Colony Color ####

# Find differentially-abundant groups ##
set.seed(123)
da_analysis <- differentialTest(formula = ~ ColonyColor, #abundance
                                phi.formula = ~ 1, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps,
                                fdr_cutoff = 0.05)

plot(da_analysis) + labs(y="Differentially abundant taxa\n(relative abundance)") +
  theme(axis.text.y = element_text(face="bold.italic"),
        axis.title.y = element_text(face="bold",size=16),
        strip.text = element_text(face="bold",size=12))
ggsave("./output/figs/colonycolor_differential_abundance_plot.png",dpi = 300,width = 14,height = 6)

sigs_colcolor <- otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = ps)


# Extract sequences from taxa with significant differential abundance
sigs_colcolor_seqs <- which(taxa_names(ps) %in% names(sigs_colcolor))
x <- ps@refseq[sigs_colcolor_seqs]
names(x) <- paste0(names(x),"_", otu_to_taxonomy(names(x),ps))
writeXStringSet(x, "./output/ColonyColor_Significant_Taxa_Seqs.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")


# Plot differentially-abundant taxa (Colony Color) ####
da_analysis$significant_taxa
# glm models with ColonyColor and no intercept

sigs_colcolor
colorblindr::palette_plot(pal.discrete)

set.seed(123)
corncob_da1 <- bbdml(formula = OTU8 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p1 <- plot(corncob_da1,color="ColonyColor") + ggtitle(sigs_colcolor[1]) + 
  scale_color_manual(values = pal.discrete[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

set.seed(123)
corncob_da2 <- bbdml(formula = OTU14 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p2 <- plot(corncob_da2,color="ColonyColor") + ggtitle(sigs_colcolor[2]) + 
  scale_color_manual(values = pal.discrete[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank()) 

set.seed(123)
corncob_da3 <- bbdml(formula = OTU24 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p3 <- plot(corncob_da3,color="ColonyColor") + ggtitle(sigs_colcolor[3]) + 
  scale_color_manual(values = pal.discrete[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

set.seed(123)
corncob_da4 <- bbdml(formula = OTU36 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p4 <- plot(corncob_da4,color="ColonyColor") + ggtitle(sigs_colcolor[4]) + 
  scale_color_manual(values = pal.discrete[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

set.seed(123)
corncob_da5 <- bbdml(formula = OTU110 ~ ColonyColor,
                     phi.formula = ~ ColonyColor,
                     data = ps)

p5 <- plot(corncob_da5,color="ColonyColor") + ggtitle(sigs_colcolor[5]) + 
  scale_color_manual(values = pal.discrete[c(11,4,6,1)]) + 
  theme(axis.text.x = element_blank())

# combine figs with patchwork
combined_plots <- p1/p2/p3/p4/p5
combined_plots

ggsave(combined_plots,filename = "./output/figs/ColonyColor_diff_abund_taxa.png",
       height = 18,width = 12,dpi = 300,device = "png")


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

# Extract sequences from taxa with significant differential abundance
sigs_tempgroup_seqs <- which(taxa_names(ps) %in% names(sigs_tempgroup))
x <- ps@refseq[sigs_tempgroup_seqs]
names(x) <- paste0(names(x),"_", otu_to_taxonomy(names(x),ps))
writeXStringSet(x, "./output/TempGroup_Significant_Taxa_Seqs2.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

temp_diff_otus <- unlist(purrr::map(str_split(names(x),"_"),1))
temp_diff_otus <- noquote(temp_diff_otus)



# Corncob plots of differential abundance
set.seed(123)
corncob_da7 <- bbdml(formula = OTU11 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p7 <- plot(corncob_da7,color = "AvgSiteTempGroup") + ggtitle(sigs_tempgroup[1]) + 
  scale_color_manual(values = pal.discrete[c(11,3,6)]) + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

set.seed(123)
corncob_da8 <- bbdml(formula = OTU12 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p8 <- plot(corncob_da8,color="AvgSiteTempGroup") + ggtitle(sigs_tempgroup[2]) + 
  scale_color_manual(values = pal.discrete[c(11,3,6)]) + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

set.seed(123)
corncob_da9 <- bbdml(formula = OTU18 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p9 <- plot(corncob_da9,color="AvgSiteTempGroup") + ggtitle(sigs_tempgroup[3]) + 
  scale_color_manual(values = pal.discrete[c(11,3,6)]) + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

set.seed(123)
corncob_da10 <- bbdml(formula = OTU19 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p10 <- plot(corncob_da10,color="AvgSiteTempGroup") + ggtitle(sigs_tempgroup[4]) + 
  scale_color_manual(values = pal.discrete[c(11,3,6)]) + 
  theme(axis.text.x = element_blank(),
        axis.title = element_text(face="bold",size=16))


set.seed(123)
corncob_da11 <- bbdml(formula = OTU21 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p11 <- plot(corncob_da11,color="AvgSiteTempGroup") + ggtitle(sigs_tempgroup[5]) + 
  scale_color_manual(values = pal.discrete[c(11,3,6)]) + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

set.seed(123)
corncob_da12 <- bbdml(formula = OTU27 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p12 <- plot(corncob_da12,color="AvgSiteTempGroup") + ggtitle(sigs_tempgroup[6]) + 
  scale_color_manual(values = pal.discrete[c(11,3,6)]) + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

set.seed(123)
corncob_da13 <- bbdml(formula = OTU28 ~ AvgSiteTempGroup,
                     phi.formula = ~ AvgSiteTempGroup,
                     data = subset_samples(ps,AvgSiteTempGroup %in% levels(ps@sam_data$AvgSiteTempGroup)))

p13 <- plot(corncob_da13,color="AvgSiteTempGroup") + ggtitle(sigs_tempgroup[7]) + 
  scale_color_manual(values = pal.discrete[c(11,3,6)]) + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")



temp_group_plot <- p7 / p8 / p9 / p10 / p11 / p12 / p13
temp_group_plot

ggsave(temp_group_plot,filename = "./output/figs/TempGroup_diff_abund_taxa.png",
       dpi=300,width = 12, height = 16,device = "png")


# Stat tables for differential abundance of OTUS by colony color and temperature ####
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
ps@sam_data$AvgSiteTempGroup

sink("./output/stats/corncob_TempGroup_tables.txt")
print(sigs_tempgroup[1])
corncob_da7
print(sigs_tempgroup[2])
corncob_da8
print(sigs_tempgroup[3])
corncob_da9
print(sigs_tempgroup[4])
corncob_da10
print(sigs_tempgroup[5])
corncob_da11
print(sigs_tempgroup[6])
corncob_da12
print(sigs_tempgroup[7])
corncob_da13
sink(NULL)

# Differentially abundant taxa, split by coral species ####

# Split data by Coral Species
Pacuta <- subset_samples(ps,SpeciesConfirmed == "P. acuta")
Pdamicornis <- subset_samples(ps,SpeciesConfirmed == "P. damicornis")


# Differential abundance test in P acuta by Colony Color 
set.seed(123)
Pacuta_da_analysis <- differentialTest(formula = ~ ColonyColor, #abundance
                                phi.formula = ~ 1, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = Pacuta,
                                fdr_cutoff = 0.05)

pacuplot <- plot(Pacuta_da_analysis) + theme(axis.title.y = element_text(size=12,face="bold"))
pacuplot
ggsave("./output/figs/Pacuta_differential_abundance_by_ColonyColor_plot.png",width = 16,height = 6)

# Differential abundance test in P damicornis by Colony Color
set.seed(123)
Pdamicornis_da_analysis <- differentialTest(formula = ~ ColonyColor, #abundance
                                       phi.formula = ~ ColonyColor, #dispersion
                                       formula_null = ~ 1, #mean
                                       phi.formula_null = ~ ColonyColor,
                                       test = "Wald", boot = FALSE,
                                       data = Pdamicornis,
                                       fdr_cutoff = 0.05)

pdamiplot <- plot(Pdamicornis_da_analysis) + theme(axis.title.y = element_text(size=12,face="bold"))
pdamiplot
ggsave("./output/figs/Pdamicornis_differential_abundance_by_ColonyColor_plot.png",width = 16,height = 6)

pacuplot / pdamiplot
ggsave("./output/figs/acuta_and_damicornis_taxa_differential_abundance_combined_plot.png",width = 16,height = 8)


# Full differential abund and dispersion test ####
# Differential abundance test, controlling for effect of colony color on dispersion and controlling
# for effect of Coral Species on abundance

#rename columns in metadata for prettier plot
names(ps@sam_data)[10] <- "CoralSpecies"

set.seed(123)
Full_da_analysis <- differentialTest(formula = ~ ColonyColor + CoralSpecies, #abundance
                                            phi.formula = ~ ColonyColor, #dispersion
                                            formula_null = ~ CoralSpecies, #mean
                                            phi.formula_null = ~ 1,
                                            test = "Wald", boot = FALSE,
                                            data = ps,
                                            fdr_cutoff = 0.05)

plot(Full_da_analysis) + theme(axis.title.y = element_text(size=12,face="bold")) +
  labs(caption = "Joint test for differential abundance and variability across Colony Color,\n
       controlling for the effect of Coral Species on abundance")
ggsave("./output/figs/Full_da_test.png",width = 16,height = 8,dpi=300)
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
  geom_jitter(height = 0, width = .2,alpha=.5) + theme_bw() +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU) +
  theme(legend.position = "bottom", axis.text.x = element_blank()) +
  scale_color_manual(values = pal.discrete[c(11,4,6,1)])

ggsave("./output/figs/most_abundant_families_by_colonycolor.png",dpi=300,width = 10,height = 6)



