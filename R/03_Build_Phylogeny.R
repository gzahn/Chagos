# Analyses of Chagos Islands coral-associated bacteria 16S Amplicons
# Building and adding a phylogeny to the cleaned phyloseq object
# Author: Geoffrey Zahn

# Packages and functions ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(phangorn)
library(msa)
library(ape)
theme_set(theme_bw())


# Read in phyloseq object ####
ps <- readRDS("./output/cleaned_ps_object.RDS")


seqs <- rownames(tax_table(ps))
names(seqs) <- seqs # This propagates to the tip labels of the tree

# Multiple sequence alignment  ####
alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 10)
# save progress 
saveRDS(alignment,"./output/trees/Chagos_16S_dna_alignment_muscle.RDS")
# re-load point
# alignment <- readRDS("./output/trees/Chagos_16S_dna_alignment_muscle.RDS")

# Convert to phangorn format
phang.align = as.phyDat(alignment, type = "DNA")

# distance max likelihood
dm <- dist.ml(phang.align)

#save
saveRDS(dm,"./output/trees/ML_Distance.RDS")

# Initial neighbor-joining tree ####
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeNJ$tip.label <- seqs
#save
saveRDS(treeNJ, "./output/trees/treeNJ.RDS")
# re-load
# treeNJ <- readRDS("./output/trees/treeNJ.RDS")

# Estimate model parameters ####
fit = pml(treeNJ, data=phang.align)
#save
saveRDS(fit,"./output/trees/fit_treeNJ.RDS")
# re-load
# fit <- readRDS("./output/trees/fit_treeNJ.RDS")

# Likelihood of tree ####
fitJC <- optim.pml(fit, TRUE)
# save
saveRDS(fitJC, "./output/trees/Chagos_16S_tree_fitJC.RDS") # This is the new tree using optim.pml
write.tree(fitJC$tree, file = "./output/trees/Chagos_16S_tree_JC.nwk")
# reload point
# fitJC <- readRDS("./output/trees/Chagos_16S_tree_fitJC.RDS")

# Bootstrap on JC model tree ####
bs = bootstrap.pml(fitJC, bs=100, control = pml.control(trace = 0),multicore = TRUE, mc.cores = 4)
saveRDS(bs,"./output/trees/Chagos_16S_tree_fitJC_bootstrap.RDS")
# re-load point
# bs <- readRDS("./output/trees/Chagos_16S_tree_fitJC_bootstrap.RDS")

# Re-fit tree with GTR model ####
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                              control = phangorn::pml.control(trace = 1L),rearrangement = "stochastic")
saveRDS(fitGTR, "./output/trees/Chagos_16S_fitGTR2.RDS") # This is the new tree using optim.pml
write.tree(fitGTR$tree, file = "./output/trees/Chagos_16S_fitGTR2_tree.nwk")
# re-load
# fitGTR = readRDS("./output/trees/Chagos_16S_fitGTR2.RDS") # loading new tree. does it work??

# Bootstrap with GTR model tree
bs = bootstrap.pml(fitGTR, bs=100, optNni=TRUE, control = pml.control(trace = 0),
                   multicore=TRUE, mc.cores = 4)
saveRDS(bs,"./output/trees/Chagos_16S_tree_fitGTR_bootstrap.RDS")
# re-load point
# bs <- readRDS("./output/trees/Chagos_16S_tree_fitGTR_bootstrap.RDS")

# Compare trees ####

# look at GTR tree
plotBS(midpoint(fitGTR$tree), bs, p = 50, type="p")

# Compare trees (JC vs GTR)
anova(fitJC, fitGTR)

detach("package:phangorn", unload=TRUE)

# add tree to phyloseq object ####
ps2 = phyloseq(tax_table(tax_table(ps)),
               otu_table(otu_table(ps)),
               sample_data(sample_data(ps)),
               phy_tree(fitGTR$tree))

ps <- ps2

# Save updated phyloseq object with tree
saveRDS(ps, "./output/cleaned_ps_object_w_tree.RDS")
