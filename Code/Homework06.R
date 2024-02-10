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
output_file <- "Homework06/amino_acid_sequence.fasta"
writeXStringSet(aa_sequence, file = output_file,
                format = "fasta", width = 60)

#### Read this file into R using the appropriate function ####
accession_numbers<- read.table("Homework06/AccNumbers.txt")

# Sample list of accession numbers
#accession_numbers <- c("AYN72248", "AYN72247", "AYN72246", "AYN72242", "KAF0875111")

# Convert the list to a character string
#accession_string <- paste(accession_numbers, collapse = ",")

# Print the formatted string
#print(accession_string)


AccessionNumbersGO <- GetProteinGOInfo(accession_numbers)
PlotGoInfo(AccessionNumbersGO)
?PlotGoInfo()

# Connect to UniProt database
up <- UniProt.ws(accession_numbers)

# Query UniProt database to retrieve protein entries
protein_entries <- getBM(attributes = c('acc', 'go_id'),
                         filters = 'acc',
                         values = accession_numbers,
                         mart = up)

# Show the retrieved GO terms
print(protein_entries)





