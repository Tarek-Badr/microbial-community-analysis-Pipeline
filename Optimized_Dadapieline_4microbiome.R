if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

library(dada2); packageVersion("dada2")
library(knitr)
library(BiocStyle)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(phangorn)

C:/Users/ngs-adm/Desktop/INTEGRATE
setwd("C:/Users/ngs-adm/Desktop/INTEGRATE")

path <- "C:/Users/ngs-adm/Desktop/INTEGRATE"

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

#Filter and trim
#Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and  maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging quality scores.
#watch out with trunclen, reads have to overlap at the end, standar script 250,200, you havr to try out, maxEE can be eased maxEE=c(2,5) if too many read are lost because of low quality

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,200),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE


head(out)

plotQualityProfile(filtFs)
plotQualityProfile(filtRs)

#Learn the Error Rates

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# apply the core sample inference algorithm to the dereplicated data.

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)

dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

#Inspecting the returned dada-class object:

dadaFs[[1]]

#Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) #min overlap is 12 as standard, but can be adjusted

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#CMost of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?

#Construct sequence table

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

#Considerations for your own data: Sequences that are much longer or shorter than expected may be the result of non-specific priming. 
#You can remove non-target-length sequences from your sequence table 
#(eg.  seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). This is analogous to "cutting a band" in-silico to get amplicons of the targeted lengt

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 288:460]  #our v3-v4 amplicon should be 460 bp

table(nchar(getSequences(seqtab2)))

#Remove chimeras
#repeat weíth seqtab2!!!!!!!!!!!!!!!!

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)

Pipeline_Track = head(track)
 write.csv(Pipeline_Track, file = "Pipeline_Track.csv")
################################################ repeating with seqtab2


seqtab.nochim2 <- removeBimeraDenovo(seqtab2, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim2)

sum(seqtab.nochim2)/sum(seqtab)

#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track2 <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim2")
rownames(track2) <- sample.names
head(track2)

################################
#Assign taxonomy
################################
ref_fasta = "C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/silva_nr_v132_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

taxa_species = addSpecies(taxa,"C:/Users/Lenovo R2G/Desktop/DADA2/Silva/Silva.nr_v132/silva_species_assignment_v132.fa")
taxa.print_s <- taxa_species # Removing sequence rownames for display only
rownames(taxa.print_s) <- NULL
head(taxa.print_s)


###############Assign Taxa with RDP
MiSeq_SOP/rdp_train_set_16.fa.gz
MiSeq_SOP/rdp_species_assignment_16.fa.gz


taxa_RDP <- assignTaxonomy(seqtab.nochim, "C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/rdp_train_set_16.fa.gz")
taxa.print_RDP  <- taxa_RDP # Removing sequence rownames for display only
rownames(taxa.print_RDP ) <- NULL
head(taxa.print_RDP)

taxa_species_RDP = addSpecies(taxa_RDP,"C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/rdp_species_assignment_16.fa.gz")
taxa.print_spp_RDP <- taxa_species_RDP # Removing sequence rownames for display only
rownames(taxa.print_spp_RDP) <- NULL
head(taxa.print_spp_RDP)


genus.species_RDP <- assignSpecies(taxa_RDP, "C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/rdp_species_assignment_16.fa.gz")
unname(genus.species_RDP)

##################################Assign Taxa with GTDB
C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/
/GTDB_bac-arc_ssu_r86.fa.gz
/GTDB_dada2_assignment_species.fa.gz

taxa_GTDB <- assignTaxonomy(seqtab.nochim, "C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_bac-arc_ssu_r86.fa.gz")

unname(taxa_GTDB)

taxa.print_GTDB  <- taxa_GTDB # Removing sequence rownames for display only

rownames(taxa.print_GTDB ) <- NULL

head(taxa.print_GTDB)


taxa_species_GTDB = addSpecies(taxa_GTDB,"C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_dada2_assignment_species.fa.gz", allowMultiple=TRUE)
taxa.print_spp_GTDB  <- taxa_species_GTDB  # Removing sequence rownames for display only
rownames(taxa.print_spp_GTDB ) <- NULL
head(taxa.print_spp_GTDB )


##################################
#Extracting the standard goods from DADA2
#The typical standard outputs from amplicon processing are a fasta file, a count table, and a taxonomy table. So here's one way we can generate those files from your DADA2 objects in R:
# giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)

asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:

asv_fasta_230200_GTDB <- c(rbind(asv_headers, asv_seqs))

write(asv_fasta_230200_GTDB, "asv_fasta_230200_GTDB.fa")

# count table:

asv_fasta_230200_GTDB <- t(seqtab.nochim)

row.names(asv_fasta_230200_GTDB) <- sub(">", "", asv_headers)

write.table(asv_fasta_230200_GTDB, "asv_fasta_230200_GTDB.tsv", sep="\t", quote=F, col.names=NA)

# tax table:

asv_tax_230200_GTDB <- taxa_GTDB

row.names(asv_tax_230200_GTDB) <- sub(">", "", asv_headers)

write.table(asv_tax_230200_GTDB, "ASVs_taxonomy_230200_GTDB.tsv", sep="\t", quote=F, col.names=NA)


######################### tax table with taxa print for better resolution

asv_tax_2 <- taxa.print_GTDB
row.names(asv_tax_2) <- sub(">", "", asv_headers)
write.table(asv_tax_2, "ASVs_taxonomy_2.tsv", sep="\t", quote=F, col.names=NA)

############ASV file with RDP (RDP is much better than SILVA)

asv_tax_RDP <- taxa.print_spp_RDP
row.names(asv_tax_RDP) <- sub(">", "", asv_headers)
write.table(asv_tax_RDP, "ASVs_taxonomy_RDP.tsv", sep="\t", quote=F, col.names=NA)

############ASV file with GTDB

asv_tax_GTDB <- taxa.print_GTDB
row.names(asv_tax_GTDB) <- sub(">", "", asv_headers)
write.table(asv_tax_GTDB, "ASVs_taxonomy_GTDB.tsv", sep="\t", quote=F, col.names=NA)


#######################################################
#Considerations for your own data: If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) and see if this fixes the assignments. If using DECIPHER for taxonomy, try IdTaxa (..., strand="both").

######################################################
#Alternatives: The recently developed IdTaxa taxonomic classification method  (i think the prior Method to assign Taxonomy is better)
######################################################

library(DECIPHER); packageVersion("DECIPHER")
dna_2 <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
#Load a training object
#SILVA_SSU_r138_2019.RData  #SILVA_SSU_r132_March2018.RData  #RDP_v16-mod_March2018.RData

load("C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_r89-mod_June2019.RData")

#load("C:/Users/Lenovo R2G/Desktop/DADA2/Silva/SILVA_SSU_r138_2019.RData") 

ids_GTDB <- IdTaxa(dna_2, trainingSet, strand="top", processors=NULL) # use all processors

print(ids_GTDB)
plot(ids_GTDB)
plot(ids_GTDB, trainingSet)
print(ids_r138)
plot(ids_r138)

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy

taxid <- t(sapply(ids_GTDB, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa_GTDB, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

taxa <- taxid 


################################
#Self remarks : some how with the training set i reach the genus level (taxa.print), but not with the full silva (taxa)??? How??
#################################
