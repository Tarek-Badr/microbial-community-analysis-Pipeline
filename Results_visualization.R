getwd()

ASV_species = read.delim("0-TopASVs_filtered_merged.txt",sep='\t')
ASV_family = read.delim("ASVs_200_NCfilt_merged_Family.txt",sep='\t')
ASV_phylum = read.delim("ASVs_200_NCfilt_merged_phylum.txt",sep='\t')


library(reshape2)

ASV_species = ASV_species[1:40]
ASV_family = ASV_family[1:40]
ASV_phylum = ASV_phylum[1:40]


#ASV_species$row <- seq_len(nrow(ASV_species))

ASV_species <- melt(ASV_species, id.vars = "Blast_ID")
ASV_family <- melt(ASV_family, id.vars = "Family")
ASV_phylum <- melt(ASV_phylum, id.vars = "Phylum")


#data$row <- seq_len(nrow(data))
#data2 <- melt(data2, id.vars = "Blast_ID")



library(tidyverse)

#data(iris)
#dummy <- iris %>% count(Species)
#ggplot(data = dummy, aes(Species, y = n, fill = Species)) +
#  geom_bar(stat = "identity")

#ggplot(data2, aes(x = variable, y = value, fill = Blast_ID)) + 
#  geom_bar(position="fill", stat = "identity") +
#  xlab("\nSample") +
#  ylab("Species\n") +
#  theme_bw()

#ggplot(ASV_species, aes(x = variable, y = value, fill = Blast_ID)) + 
#  geom_bar(position="fill", stat = "identity") +
#  xlab("\nSample") +
#  ylab("Species\n") +
#  theme_bw()

library(ggplot2)


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




ggplot(data2, aes(x = variable, y = value, fill = row)) + 
  geom_bar(position="fill", stat = "identity") +
  xlab("\nSample") +
  ylab("Species\n") +
  theme_bw()


p = ggplot(data2, aes(x = variable, y = value, fill = Blast_ID)) + 
  geom_bar(position="fill", stat = "identity") +
  xlab("\nSample") +
  ylab("Species\n") +
  theme_bw()

p = ggplot(ASV_species, aes(x = variable, y = value, fill = Blast_ID)) + 
  geom_bar(position="fill", stat = "identity") +
  xlab("\nSample") +
  ylab("Species\n") +
  theme_bw()


p + theme(legend.position="top", legend.background = element_rect(fill="lightblue", 
                                                                  size=0.5, linetype="solid"))





#####################################
#Making Heatmaps for Mutations per CLC gene
#####################################

Top_ASVs = read.delim("0-TopASVs_filtered_merged.txt",sep='\t',row.names = 1)
Top_ASVs


Mutec_Variants = read.table("Mutec_Variants.txt",sep='\t', header = T)
Mutec_Variants


#normalization (as samples size < 30, rld is prefered than vst)

library("pheatmap")
library("RColorBrewer")
library(scales)
library(ggplot2)


################Heatmap of mutations per CRC highly mutated genes
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

wxs_crc_norm <- t(apply(wxs_crc, 1, cal_z_score))
pheatmap(wxs_crc_norm)


pheatmap(wxs_crc_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = drows1, 
         clustering_distance_cols = dcols1, angle_col = 45) 


################Heatmap of mutations according to mutation type

Mutations_HM = read.delim("Mutec_Variants.txt",sep='\t',row.names = 1)

Mutations_HM

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

Mutations_HM_norm <- t(apply(Mutations_HM, 1, cal_z_score))
pheatmap(Mutations_HM_norm)


pheatmap(Mutations_HM_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = TRUE, cluster_cols = TRUE) 

####################################
#Adjust Table for plotting and Statistical significance
####################################

Mutations = read.delim("Mutations.txt",sep='\t', header = T)
Mutations

head(Mutations)

df.long <- pivot_longer(Mutations, cols=-1, names_to = "Variants", values_to = "Count")
df.long

df.long$Count <- as.numeric(as.vector(df.long$Count))
ggplot(data=df.long, aes(x=Variants, y=Count, fill=Mouse)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

df.long$Count <- as.numeric(as.vector(df.long$Count))

##################plotting

P = ggplot(data = df.long, aes(Variants, Count)) +
  geom_boxplot(aes(colour = Mouse))

sp = P + theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                          size = 12, angle = 45, hjust = 1),
               axis.text.y = element_text(face = "bold", color = "blue", 
                                          size = 12, angle = 45))

###best one##
SD = sp + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))
SD
#####box Plots
P_box = ggboxplot(df.long, x = "Variants", y = "Count", color = "Mouse", palette = "jco",
                  add = "jitter")

P_box
sp_box = P_box + theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                                  size = 12, angle = 45, hjust = 1),
                       axis.text.y = element_text(face = "bold", color = "blue", 
                                                  size = 12, angle = 45))
sp_box

SD_box = sp_box + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))
SD_box





#############################################statistical testing
sp  + annotation_logticks() 


SD_box + stat_compare_means(aes(group = Mouse))
SD_box + stat_compare_means(aes(group = Mouse), label =  "p.signif", label.x = 1.5) #wilcoxon test

##############################################################################################################
#Plotting the vcf files in waterfall plots and mutation load plots
#############################################################################################################














