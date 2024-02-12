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


#Set Working directory for Data and Homework folder ####
setwd("Bioinformatics/")
setwd("Data/")
setwd("Homework06/")

#Checked that directory was in place ####
getwd()

mySequences01 <- readDNAStringSet("sequence01.fasta")
mySequences01

#### Translate DNA sequence into AA sequence ####
aa_sequence <- translate(mySequences01)
aa_sequence
as.character(aa_sequence)

#### Writing a aa sequence into a fasta file ####
output_file <- "amino_acid_sequence.fasta"
writeXStringSet(aa_sequence, file = output_file,
                format = "fasta", width = 60)

#### Read this file into R using the appropriate function ####
accession_numbers<- read.table("AccNumbers.txt")

# Sample list of accession numbers
accession_numbers <- c("A0A6G1AHE8", "A0A7J7SV63", "P21445", "Q2F7I8", "A0A7J7U5J2")

# Convert the list to a character string
accession_string <- paste(accession_numbers, collapse = ",")

# Print the formatted string
print(accession_string)

#Reading accession numbers into GetProteinGOInfo ####
AccessionNumbersGO <- GetProteinGOInfo(accession_numbers)
str(AccessionNumbersGO)
#write.csv(AccessionNumbersGO, "AccessionNumbersGO.csv", row.names = FALSE)

PlotGoInfo(AccessionNumbersGO)
?PlotGoInfo()



# Assuming AccessionNumbersGO is your original data frame
# Extract GO terms and their counts from AccessionNumbersGO
go_terms <- unlist(strsplit(AccessionNumbersGO$Gene.Ontology..GO., ";"))
go_terms <- gsub("\\[.*?\\]", "", go_terms)  # Remove GO IDs from GO terms

# Create a data frame with GO terms and their counts
go_counts <- data.frame(GoTerm = go_terms, Count = rep(1, length(go_terms)))

# Summarize the counts for each GO term
go_counts <- aggregate(Count ~ GoTerm, go_counts, sum)

# Plot the GO information
barplot(go_counts$Count, names.arg = go_counts$GoTerm,
        xlab = "GO Terms", ylab = "Count", main = "GO Term Distribution")


