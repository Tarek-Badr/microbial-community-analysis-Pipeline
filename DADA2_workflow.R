if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")

library(dada2); packageVersion("dada2")
library(knitr)
library(BiocStyle)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(dendextend)

setwd("C:/Users/ngs-adm/Desktop/INTeGRATE-ALL/0-INTEGRATE-Analysis")
path <- "C:/Users/ngs-adm/Desktop/INTeGRATE-ALL/00-Final Ilumina Analysis/0-Samples_1-50_merged_Fastq"

list.files(path)


# assign Forward and reverse fastq files
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq (samplename doesn't include "_")

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

#Filter and trim
#Assign the filenames for the filtered fastq.gz files and Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and  maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging quality scores.
#!!!!watch out with trunclen, reads have to overlap at the end, standard script was 250,200 you havr to try out, for our primer pair V3-V4 (Amplicon size of 460bp i went with 280,230 after many trials) 
#!!!!maxEE can be eased maxEE=c(2,5) if too many read are lost because of low quality

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
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

# apply the DADA2 algorithm to the dereplicated data.

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

#Inspecting the returned dada-class object:

dadaFs[[1]]

#Merge paired reads
#min overlap is 20 as standard, but can be adjusted, if you have perfect overlapping primers then you can increase it (minOverlap = 100) or so

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) 

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited in case you trimed away the overlap between your reads?
#Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:387]
table(nchar(getSequences(seqtab2)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
Pipeline_Track = head(track)
write.csv(track, file = "Pipeline_Track_final.csv")
write.csv(sample.names, file = "sample.names.csv")

################################
#Assign taxonomy
################################

ref_fasta = "C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz"

taxa_GTDB <- assignTaxonomy(seqtab.nochim, ref_fasta,tryRC = TRUE )

unname(taxa_GTDB)

taxa.print_GTDB  <- taxa_GTDB 

rownames(taxa.print_GTDB ) <- NULL

head(taxa.print_GTDB)

ref_fasta_species = "C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_bac120_arc122_ssu_r202_Species.fa.gz"

taxa_species_GTDB = addSpecies(taxa_GTDB,ref_fasta_species, allowMultiple=TRUE)

taxa.print_spp_GTDB  <- taxa_species_GTDB 

rownames(taxa.print_spp_GTDB ) <- NULL

head(taxa.print_spp_GTDB )

##################################
#Creating Output files from DADA2
##################################
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_230_GTDB <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta_230_GTDB, "asv_fasta_230_GTDB.fa")
# count table:
asv_tab_230_GTDB <- t(seqtab.nochim)
row.names(asv_tab_230_GTDB) <- sub(">", "", asv_headers)
write.table(asv_tab_230_GTDB, "ASVs_counts_230_GTDB.tsv", sep="\t", quote=F, col.names=NA)
# tax table:
asv_tax_230_GTDB <- taxa_GTDB
row.names(asv_tax_230_GTDB) <- sub(">", "", asv_headers)
write.table(asv_tax_230_GTDB, "ASVs_taxonomy_230_GTDB.tsv", sep="\t", quote=F, col.names=NA)

###################Analysis with DESEq2
#We read the output files from DADA2 and the sample information file
###################

count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("ASVs_taxonomy_GTDB.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
tax_tab_RDP <- as.matrix(read.table("ASVs_taxonomy_RDP.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

sample_info_tab <- samdf

###################################
#Make the Deseq Tables
###################################

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~Protocol) 
#deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
#############################
#clustering
#############################
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust) 

################################
#Phyloseq
################################
#Rarefaction curves
################################
library(vegan)
rarecurve(t(count_tab), step=100, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(t(count_tab)))))
########################################
#Richness and diversity estimates
########################################

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)  
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

plot_richness(ASV_physeq, color="Protocol", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "Simpson")) 
ASV_physeq.prop <- transform_sample_counts(ASV_physeq, function(otu) otu/sum(otu))
ASV_physeq.nmds.bray <- ordinate(ASV_physeq.prop, method="NMDS", distance="bray")     
plot_ordination(ASV_physeq.prop, ord.nmds.bray, label = "Protocol" , title="Bray NMDS")

#######################################################################################
#Taxonomic summaries
#######################################################################################

#We now construct a phyloseq object directly from the dada2 outputs.

dna <- Biostrings::DNAStringSet(taxa_names(ASV_physeq))
names(dna) <- taxa_names(ASV_physeq)
ps_GTDB <- merge_phyloseq(ASV_physeq, dna)
taxa_names(ps_GTDB) <- paste0("ASV", seq(ntaxa(ps_GTDB)))
ps_GTDB

#########################################################################################


                     ###########################################################################
                                            #Hope this was useful to you#
                                                  ##Tarek Badr##
                     ###########################################################################               
 
