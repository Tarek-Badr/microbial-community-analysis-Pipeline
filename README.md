# microbial-community-analysis-Pipeline
This repository details the analysis of microbial communities in clinical microbiological samples
The analysis starts from raw fastq files where used primers (and there reverse compliment if applicable) get trimmed using cutadapt 
and further analysis is carried out using DADA2 to end up with  amplicon sequence variant (ASV) table to be further used for
DESeq2 analysis. 
