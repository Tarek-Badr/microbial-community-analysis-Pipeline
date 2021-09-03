library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(tidyverse)
library(reshape2)

getwd()

ASV_species = read.delim("0-TopASVs_filtered_merged.txt",sep='\t')
ASV_family = read.delim("ASVs_200_NCfilt_merged_Family.txt",sep='\t')
ASV_phylum = read.delim("ASVs_200_NCfilt_merged_phylum.txt",sep='\t')



ASV_species = ASV_species[1:40]
ASV_family = ASV_family[1:40]
ASV_phylum = ASV_phylum[1:40]


#ASV_species$row <- seq_len(nrow(ASV_species))

ASV_species <- melt(ASV_species, id.vars = "Blast_ID")
ASV_family <- melt(ASV_family, id.vars = "Family")
ASV_phylum <- melt(ASV_phylum, id.vars = "Phylum")



p_Phylum = ggplot(ASV_phylum, aes(x = variable, y = value, fill = Phylum)) + 
  geom_bar(position="fill", stat = "identity") +
  xlab("\nSample") +
  ylab("Phylum\n") +
  theme_bw()

p_Phylum + theme(legend.position="top", legend.background = element_rect(fill="lightblue", 
                                                                  size=0.5, linetype="solid"))


p_fam = ggplot(ASV_family, aes(x = variable, y = value, fill = Family)) + 
  geom_bar(position="fill", stat = "identity") +
  xlab("\nSample") +
  ylab("Family\n") +
  theme_bw()

p_fam + theme(legend.position="top", legend.background = element_rect(fill="lightblue", 
                                                                      size=0.5, linetype="solid"))

p_sp = ggplot(ASV_species, aes(x = variable, y = value, fill = Blast_ID)) + 
  geom_bar(position="fill", stat = "identity") +
  xlab("\nSample") +
  ylab("Species\n") +
  theme_bw()

p_sp + theme(legend.position="top", legend.background = element_rect(fill="lightblue", 
                                                                     size=0.5, linetype="solid"))


#####################################
#Making Heatmaps
#####################################

Top_ASVs = read.delim("0-TopASVs_filtered_merged.txt",sep='\t',row.names = 1)
Top_ASVs

pheatmap(Top_ASVs)

pheatmap(Top_ASVs, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = TRUE, cluster_cols = FALSE,
) 

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

tasv <- t(apply(Top_ASVs, 1, cal_z_score))
pheatmap(wxs_crc_norm)


pheatmap(tasv, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = drows1, 
         clustering_distance_cols = dcols1, angle_col = 45) 
