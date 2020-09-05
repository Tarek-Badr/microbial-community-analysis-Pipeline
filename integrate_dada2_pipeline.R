if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

BiocManager::install("ShortRead")

library(dada2); packageVersion("dada2")
library(knitr)
library(BiocStyle)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(phangorn)


setwd("C:/Users/ngs-adm/Desktop/INTeGRATE-ALL/0-INTEGRATE-Analysis")

path <- "C:/Users/ngs-adm/Desktop/INTeGRATE-ALL/0-INTEGRATE-Analysis"

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)


#primer check
fwd_primer <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAGAGTTTGATCCTGGCTCAG"
rev_primer <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGCTGCCTCCCGTAGGAGT"
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))

# This function counts number of reads in which the primer is found

count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  return(sum(num_hits > 0))
}

count_primers(fwd_primer, fnFs[[1]])

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
#maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
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

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:387]

dim(seqtab)
dim(seqtab2)
# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

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
write.csv(Pipeline_Track, file = "Pipeline_Track_2.csv")


################################
#Assign taxonomy with GTDB
################################
C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/
/GTDB_bac-arc_ssu_r86.fa.gz
/GTDB_dada2_assignment_species.fa.gz

ref_fasta = "C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_bac-arc_ssu_r86.fa.gz"

taxa_GTDB <- assignTaxonomy(seqtab.nochim, ref_fasta )

unname(taxa_GTDB)

taxa.print_GTDB  <- taxa_GTDB # Removing sequence rownames for display only

rownames(taxa.print_GTDB ) <- NULL

head(taxa.print_GTDB)

ref_fasta_species = "C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_dada2_assignment_species.fa.gz"
taxa_species_GTDB = addSpecies(taxa_GTDB,ref_fasta_species, allowMultiple=TRUE)
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

asv_fasta_200200_GTDB <- c(rbind(asv_headers, asv_seqs))

write(asv_fasta_200200_GTDB, "asv_fasta_200200_GTDB.fa")

# count table:

asv_tab_200200_GTDB <- t(seqtab.nochim)
row.names(asv_tab_200200_GTDB) <- sub(">", "", asv_headers)
write.table(asv_tab_200200_GTDB, "ASVs_counts_200200_GTDB.tsv", sep="\t", quote=F, col.names=NA)


# tax table:

asv_tax_200200_GTDB <- taxa_GTDB
row.names(asv_tax_200200_GTDB) <- sub(">", "", asv_headers)
write.table(asv_tax_200200_GTDB, "ASVs_taxonomy_200200_GTDB.tsv", sep="\t", quote=F, col.names=NA)


######################### tax table with taxa print for better resolution

asv_tax_2 <- taxa.print_GTDB
row.names(asv_tax_2) <- sub(">", "", asv_headers)
write.table(asv_tax_2, "ASVs_taxonomy_2.tsv", sep="\t", quote=F, col.names=NA)

#######################################################
#Considerations for your own data: If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) and see if this fixes the assignments. If using DECIPHER for taxonomy, try IdTaxa (..., strand="both").

#######################################################
#Removing likely contaminants
#######################################################
library(decontam)
packageVersion("decontam") # 1.1.2 when this was put together

vector_for_decontam <- c(rep(FALSE, 18), rep(TRUE, 1), rep(FALSE, 1))

contam_df <- isContaminant(t(asv_tab_230200_GTDB), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 6 as contaminants

# getting vector holding the identified contaminant IDs

contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

asv_tax_230200_GTDB[row.names(asv_tax_230200_GTDB) %in% contam_asvs, ]

# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
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



library(dada2); packageVersion("dada2")
library(knitr)
library(BiocStyle)
library(gridExtra)
library(DECIPHER)
library(phangorn)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

install.packages('dendextend')

library(dendextend)

###################Analysis with DESEq2
#We read the output files from DADA2 and the sample information file
###################
C:/Users/ngs-adm/Desktop/INTEGRATE

count_tab <- read.table("C:/Users/ngs-adm/Desktop/INTEGRATE/ASVs_counts_200200_GTDB.tsv", header=T, row.names=1,
                          check.names=F, sep="\t")

tax_tab <- as.matrix(read.table("C:/Users/ngs-adm/Desktop/INTEGRATE/ASVs_taxonomy_200200_GTDB.tsv", header=T,
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

########clustering
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust) 



########
#If we have different colors in a column in our sample file then use this clustering to apply it
#####
#euc_dend <- as.dendrogram(euc_clust, hang=0.1)
#dend_cols <- as.character(sample_info_tab$Protocol[order.dendrogram(euc_dend)])
#labels_colors(euc_dend) <- dend_cols
#plot(euc_dend, ylab="VST Euc. dist.")


#############################
#Ordination
############################
# making our phyloseq object with transformed table

vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="char") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info_tab$Protocol[order(sample_info_tab$Protocol)])) + 
  theme(legend.position="none")

################################
#Phyloseq
###############################
#Rarefaction curves
############################

library(vegan)

rarecurve(t(count_tab), step=100, lwd=2, ylab="ASVs", label=F)

# and adding a vertical line at the fewest seqs in any sample

abline(v=(min(rowSums(t(count_tab)))))

########################################
#Richness and diversity estimates
########################################

# first we need to create a phyloseq object using our un-transformed count table

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)

tax_tab_phy <- tax_table(tax_tab)    #here based on the GTDB database (we can retry it for RDP or SILVA)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

# and now we can call the plot_richness() function on our phyloseq object

plot_richness(ASV_physeq, color="Protocol", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "Simpson")) 

# Transform data to proportions as appropriate for Bray-Curtis distances

ASV_physeq.prop <- transform_sample_counts(ASV_physeq, function(otu) otu/sum(otu))
ASV_physeq.nmds.bray <- ordinate(ASV_physeq.prop, method="NMDS", distance="bray")     

#plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

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
#Bar plot:
##############################################################################


top20 <- names(sort(taxa_sums(ps_GTDB), decreasing=TRUE))[1:20]  # adjust number to wished top ones
ps.top20 <- transform_sample_counts(ps_GTDB, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
#plot_bar(ps.top20, x="Day", fill="family") + facet_wrap(~When, scales="free_x")  #Family in SILVA, family in Decipher

BarPlot_Top20SPP_GTDB = plot_bar(ps.top20, fill="Species")  #Family in SILVA, family in Decipher
BarPlot_Top20Genus_GTDB = plot_bar(ps.top20, fill="Genus")
BarPlot_Top20Family_GTDB = plot_bar(ps.top20, fill="Family")

top50 <- names(sort(taxa_sums(ps_GTDB), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps_GTDB, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)

#plot_bar(ps.top20, x="Day", fill="family") + facet_wrap(~When, scales="free_x")  #Family/ Genus / Species etc,....

plot_bar(ps.top50, fill="Genus")  #Family in SILVA, family in Decipher
plot_bar(ps.top50, fill="Species")
plot_bar(ps.top20, fill="Genus")
############################################################################
#Heat Maps
############################################################################

theme_set(theme_bw())

rank_names(ps.top50)

plot_heatmap(ps.top50, sample.label="Protocol","Family")

plot_heatmap(ps.top50, method = "NMDS", distance = "bray",
             sample.label = "Protocol", taxa.label = "Genus", low = "#000033",
             high = "#FF3300", na.value = "black",
             max.label = 250, title = "Heatmap of top 50 Genus (GTDB)", sample.order = NULL, taxa.order = "Genus",
             first.sample = NULL, first.taxa = NULL)

plot_heatmap(ps.top50, sample.label = "Protocol", taxa.label = "Genus", low = "#000033",
             high = "#FF3300", na.value = "black",
             max.label = 250, title = "Heatmap of top 50 Genus (GTDB)", sample.order = NULL, taxa.order = "Genus",
             first.sample = NULL, first.taxa = NULL)

plot_heatmap(ps.top50, taxa.label = "Species", low = "#000033",
             high = "#FF3300", na.value = "black",
             max.label = 250, title = "Heatmap of top 50 Species (GTDB)")

heatmap(otu_table(ps.top50))

p_HM <- plot_heatmap(ps_GTDB, "NMDS", "bray", "Protocol", "Family")


















