if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings", force = TRUE)
BiocManager::install("GenomicAlignments", force = TRUE)
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
install.packages("seqinr")
library(seqinr)
library(ape)
install.packages("UniprotR")
install.packages("protti")
install.packages("r3dmol")
library(UniprotR)
library(protti)
library(r3dmol)
library(phangorn)


#Set Working directory for Data and Homework folder ####
setwd("GitHub/Bioinformatics/Data/Homework06/")


#Checked that directory was in place ####
getwd()

mySequences01 <- readDNAStringSet("Homework06/sequence01.fasta")


#### Translate DNA sequence into AA sequence ####
dna_sequences <- readDNAStringSet("Homework06/sequence01.fasta")
aa_sequence <- translate(dna_sequences)
aa_sequence
as.character(aa_sequence)

#### Writing a aa sequence into a fasta file ####
?writeXStringSet()
writeXStringSet(aa_sequence)

readAAStringSet(filepath = "Homework06", format="fasta",
                nrec=-1L, skip=0L, seek.first.rec=FALSE,
                use.names=TRUE, with.qualities=FALSE)

output_file <- "Homework06/amino_acid_sequence.fasta"

# Write the amino acid sequence to a .fasta file
writeXStringSet(aa_sequence, file = output_file,
                format = "fasta", width = 60)





