#Loading packages ####
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
library(seqinr)
library(ape)


#Set Working directory for Data and Midterm folder ####
setwd("Bioinformatics/")
setwd("Data/Midterm01/")
#Checking to see if working directory is correct
getwd()
#---

#1. Reading fasta file and print to view the file ####
mySequences <- readDNAStringSet("sequences.fasta")
print(mySequences, show = "complete")
View(mySequences)

#Aligning the sequences
msa_alignment <- msa(mySequences)
print(msa_alignment, show = "complete")
#---

#2. Checking for differences between sequences ####
#Homo_sapiens_6 has one less bp than the other sequences
#Homo_sapiens_6 also has a mutation which is AAG instead of AAA, but this is a silent mutation, that is all I see so far
#
#---

#3. Search to see which gene these sequences represent ####
#I used BLAST to determine that this is the hemoglobin subunit beta gene
#The accession number that best matched my search for sequence 20 is LC121775, found on GenBank
#---

#4. Finding the outlier in the sequences ####
HGS <- msaConvert(msa_alignment, type = "seqinr::alignment")

#Compute a distance matrix using the 'dist.alignment' function
distance_matrix <- dist.alignment(HGS)

#Print the distance matrix
print(distance_matrix)
#Homo_sapiens_6 is the most distantly related

#Translating 6 and writing a new fasta file for seqeunce 6
# Homo_sapiens_6 <-"AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG"
Homo_sapiens_6 <- mySequences$Homo_sapiens_6
#printing to view Homo_sapiens_6
print(Homo_sapiens_6)
#Conflict between translate in seqinr and biostrings
detach("package:seqinr", unload = TRUE) # Unload seqinr if it's loaded

#Translate the DNA sequence to amino acids
aa_sequence <- Biostrings::translate(Homo_sapiens_6)
#Error doesn't impact the translation

#Showing the translated amino acid sequence completely with gaps
print(aa_sequence)
print(toString(aa_sequence))

#Writing a aa sequence into a fasta file
aa_xstring <- AAStringSet(aa_sequence)

#Specify the output file name
output_file <- "amino_acid_sequence.fasta"

#Write the XStringSet object to the FASTA file
writeXStringSet(aa_xstring, file = output_file, format = "fasta", width = 60)
#Done and is in "C:/Users/18638/OneDrive/Desktop/GitHub/Bioinformatics/Data/Midterm01"
#---

#5. Find a match to your protein sequence ####
#The accession number is KAI2558340, 100% identity, hemoglobin subunit beta
#---

#6. Finding the diseases associated with the gene #####
#According to OMIM, the diseases associated with this gene are sickle cell, erythrocytosis, and Methemoglobinemia, just to name a few.
#Homo_sapiens_6 has one of these diseases. Probably sickle cell.
#---

#7. 3D model ####
#The 3D model is in "C:/Users/18638/OneDrive/Desktop/GitHub/Bioinformatics/Output/Midterm01"





















