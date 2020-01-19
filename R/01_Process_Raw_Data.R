# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr)
library(tidyr)
library(ggplot2)
library(readxl)
library(decontam)
library(phyloseq)

# metadata
meta <- readxl::read_xlsx("./metadata/Chagos_PlateMpas_MetaData_1_Geoff.xlsx")

# Find raw fastq files and prepare workspace ####
path <- "./raw_data" 
list.files(path, full.names = TRUE)

# Parse fwd and rev reads
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Get Sample Names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3)

# Peek at quality profiles
 plotQualityProfile(fnFs[c(12,59)]) # fwd reads
 plotQualityProfile(fnRs[c(12,59)]) # rev reads

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# make new directory for filtered files
dir.create(file.path(path,"filtered"))

# Filter and trim ####
# clean out any files from previous iterations
file.remove(filtFs); file.remove(filtRs)

# cut fwd reads at 250 and rev reads at 160
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, matchIDs = TRUE)

# reassign filts for any potentially lost samples
filtFs <- list.files("raw_data/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("raw_data/filtered", pattern = "_R_filt", full.names = TRUE)

# Get Sample Names
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)


# which ones were lost?
start = unlist(map(strsplit(fnFs, split = "_"),4))
end = unlist(map(strsplit(unlist(map(strsplit(filtFs, split = "/"),3)),"_"),1))
failed <- as.numeric(sample.names[which(start %in% end == FALSE)])
failures <- meta[meta$`Pooled library ID (from GIS)` %in% failed,]
write.csv(failures,"./metadata/failures.csv",row.names = FALSE, quote = FALSE)



# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+10)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e+10)

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

# Sample Inferrence ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# save progress
dir.create("./output")
saveRDS(dadaFs,"./output/dadaFs.RDS")
saveRDS(dadaRs,"./output/dadaRs.RDS")

# Merge fwd and rev reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, trimOverhang = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./output/seqtab.nochim.RDS")

# remove all ASVs that don't have at least 10 hits
seqtab.nochim <- seqtab.nochim[,(colSums(seqtab.nochim) > 9)]

# Assign taxonomy - Silva v132 exact match / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.nochim,"./taxonomy/silva_nr_v132_train_set.fa.gz", minBoot = 80,multithread = TRUE)
taxa <- addSpecies(taxa, "./taxonomy/silva_species_assignment_v132.fa.gz")
saveRDS(taxa,"./output/taxa.RDS")

# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df) <- unlist(map(strsplit(row.names(seqtab.df), "_"),1))

# Create phyloseq object ####

# subset to remove missing samples
meta = readxl::read_xlsx("./metadata/Chagos_PlateMpas_MetaData_1_Geoff.xlsx")
in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Pooled library ID (from GIS)` == TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]
dim(seqtab.nochim)
in.seqtab = which(meta$`Pooled library ID (from GIS)` %in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

# re-order
meta = meta[order(meta$`Pooled library ID (from GIS)`),]
row.names(meta) <- meta$`Pooled library ID (from GIS)`
seqtab.nochim <- (seqtab.nochim[row.names(meta),])
identical(row.names(seqtab.nochim), as.character(meta$`Pooled library ID (from GIS)`))

# check dimensions
dim(taxa)
dim(seqtab.nochim)
dim(meta)

# make phyloseq object ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))
# save it
saveRDS(ps, "./output/phyloseq_object_16S.RDS")

# Identify and remove contaminants
blanks = which(ps@sam_data$date == "Blank")
contamdf.prev <- isContaminant(ps, neg=blanks, threshold = 0.1)
table(contamdf.prev$contaminant)

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
identical(ps, ps.noncontam)
    # no contaminants from blanks found via prevalence method in other samples!

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./output/phyloseq_object_16S_noncontam.RDS")

