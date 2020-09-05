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

C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/


count_tab <- read.table("C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")

tax_tab <- as.matrix(read.table("C:/Users/Lenovo R2G/Desktop/DADA2/MiSeq_SOP/ASVs_taxonomy_GTDB.tsv", header=T,
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







