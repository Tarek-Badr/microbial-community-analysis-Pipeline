library(dada2)
library(knitr)
library(BiocStyle)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ShortRead)
library(decontam)
library(dendextend)
library(vegan)

setwd("C:/Users/ngs-adm/Desktop/INTeGRATE-ALL/0-INTEGRATE-Analysis")
path <- "C:/Users/ngs-adm/Desktop/INTeGRATE-ALL/0-INTEGRATE-Analysis"

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

#primer check
fwd_primer <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAGAGTTTGATCCTGGCTCAG"
rev_primer <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGCTGCCTCCCGTAGGAGT"
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))

count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  return(sum(num_hits > 0))
}

count_primers(fwd_primer, fnFs[[1]])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 

head(out)

plotQualityProfile(filtFs)
plotQualityProfile(filtRs)

#Learn Error Rates
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

#Merge reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) #min overlap is 12 as standard, but can be adjusted

head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)  
dim(seqtab)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:387]
dim(seqtab2)

table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
Pipeline_Track = (track)
write.csv(Pipeline_Track, file = "Pipeline_Track_2.csv")

################################
#Assign taxonomy with GTDB
################################

ref_fasta = "C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_bac-arc_ssu_r86.fa.gz"

taxa_GTDB <- assignTaxonomy(seqtab.nochim, ref_fasta )
unname(taxa_GTDB)
taxa.print_GTDB  <- taxa_GTDB # Removing sequence rownames for display only
rownames(taxa.print_GTDB ) <- NULL
head(taxa.print_GTDB)

ref_fasta_species = "C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_dada2_assignment_species.fa.gz"

taxa_species_GTDB = addSpecies(taxa_GTDB,ref_fasta_species, allowMultiple=TRUE)
taxa.print_spp_GTDB  <- taxa_species_GTDB 
rownames(taxa.print_spp_GTDB ) <- NULL
head(taxa.print_spp_GTDB )

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta_GTDB <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta_GTDB, "asv_fasta_GTDB.fa")

# count table:
asv_tab_200200_GTDB <- t(seqtab.nochim)
row.names(asv_tab_200200_GTDB) <- sub(">", "", asv_headers)
write.table(asv_tab_200200_GTDB, "ASVs_counts_GTDB.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax_GTDB <- taxa_GTDB
row.names(asv_tax_GTDB) <- sub(">", "", asv_headers)
write.table(asv_tax_GTDB, "ASVs_taxonomy_GTDB.tsv", sep="\t", quote=F, col.names=NA)

#######################################################
#Removing contaminantions
#######################################################

vector_for_decontam <- c(rep(FALSE, 23), rep(TRUE, 1), rep(FALSE, 1))
contam_df <- isContaminant(t(asv_tab_GTDB), neg=vector_for_decontam)

table(contam_df$contaminant) 

contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
asv_tax_GTDB[row.names(asv_tax_GTDB) %in% contam_asvs, ]

# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]
# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]
# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)

######################################################

dna_2 <- DNAStringSet(getSequences(seqtab.nochim)) 
load("C:/Users/ngs-adm/Desktop/R/DADA2/MiSeq_SOP/GTDB_r89-mod_June2019.RData")
ids_GTDB <- IdTaxa(dna_2, trainingSet, strand="top", processors=NULL) 
print(ids_GTDB)
plot(ids_GTDB)
plot(ids_GTDB, trainingSet)

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

taxid <- t(sapply(ids_GTDB, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa_GTDB, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa <- taxid 


###################
Analysis with DESEq2
###################

count_tab <- read.table("C:/Users/ngs-adm/Desktop/INTEGRATE/ASVs_counts_GTDB.tsv", header=T, row.names=1,
                          check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("C:/Users/ngs-adm/Desktop/INTEGRATE/ASVs_taxonomy_GTDB.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
sample_info_tab <- samdf

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~Protocol) 
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
########clustering
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust) 

#############################
#Ordination
############################
#phyloseq object

vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

#visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="char") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info_tab$Protocol[order(sample_info_tab$Protocol)])) + 
  theme(legend.position="none")

################################
#Rarefaction curves
############################

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

dna <- Biostrings::DNAStringSet(taxa_names(ASV_physeq))
names(dna) <- taxa_names(ASV_physeq)
ps_GTDB <- merge_phyloseq(ASV_physeq, dna)
taxa_names(ps_GTDB) <- paste0("ASV", seq(ntaxa(ps_GTDB)))
ps_GTDB

##############################################################################
#Bar plot:
##############################################################################


top20 <- names(sort(taxa_sums(ps_GTDB), decreasing=TRUE))[1:20]  
ps.top20 <- transform_sample_counts(ps_GTDB, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
BarPlot_Top20SPP_GTDB = plot_bar(ps.top20, fill="Species")  
BarPlot_Top20Genus_GTDB = plot_bar(ps.top20, fill="Genus")
BarPlot_Top20Family_GTDB = plot_bar(ps.top20, fill="Family")
top50 <- names(sort(taxa_sums(ps_GTDB), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps_GTDB, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
plot_bar(ps.top50, fill="Genus") 
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
