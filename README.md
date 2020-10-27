# microbial-community-analysis-Pipeline

Copyright Mohamed Tarek Badr , Date: 27 oct 2020
University of Freiberg Germany 

This code was used for the analysis of 16srDNA sequencing samples from the INTeGRATE clinical study.
This repository details the analysis of microbial communities in clinical microbiological samples
The analysis starts from raw fastq files where used primers (and there reverse compliment if applicable) get trimmed using cutadapt 
and further analysis is carried out using DADA2 to end up with  amplicon sequence variant (ASV) table to be further used for
DESeq2 and phyloseq analysis. 
