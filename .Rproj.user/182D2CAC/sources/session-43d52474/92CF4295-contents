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

#3. Translate DNA sequence into AA sequence ####
aa_sequence <- translate(mySequences01)
aa_sequence
as.character(aa_sequence)

# Writing a aa sequence into a fasta file
output_file <- "amino_acid_sequence.fasta"
writeXStringSet(aa_sequence, file = output_file,
                format = "fasta", width = 60)
?writeXStringSet
#4. Read this file into R using the appropriate function ####
accession_numbers<- read.table("AccNumbers.txt")

#5. Sample list of accession numbers ####
accession_numbers <- c("A0A6G1AHE8", "A0A7J7SV63", "P21445", "Q2F7I8", "A0A7J7U5J2")

# Convert the list to a character string
accession_string <- paste(accession_numbers, collapse = ",")

# Print the formatted string
print(accession_string)

#6. Reading accession numbers into GetProteinGOInfo ####
AccessionNumbersGO <- GetProteinGOInfo(accession_numbers)
str(AccessionNumbersGO)
#write.csv(AccessionNumbersGO, "AccessionNumbersGO.csv", row.names = FALSE)

#7. Extract GO terms and their counts from AccessionNumbersGO ####
View(AccessionNumbersGO)
df <-- read.csv("Data/Homework06/AccessionNumbersGO.csv")
#PlotGoInfo(AccessionNumbersGO) #--> Did not work
go_terms <- unlist(strsplit(AccessionNumbersGO$Gene.Ontology..GO., ";"))
go_terms <- gsub("\\[.*?\\]", "", go_terms)  # Remove GO IDs from GO terms

# Create a data frame with GO terms and their counts
go_counts <- data.frame(GoTerm = go_terms, Count = rep(1, length(go_terms)))

# Summarize the counts for each GO term
go_counts <- aggregate(Count ~ GoTerm, go_counts, sum)

# Plot the GO information
barplot(go_counts$Count, names.arg = go_counts$GoTerm,
        xlab = "GO Terms", ylab = "Count", main = "GO Term Distribution")

#9. What are some interesting GO terms for your gene?
#virion membrane
#membrane
#viral envelope
#host cell plasma membrane
#fusion of virus membrane with host plasma membrane




#10. Use GetPathology_Biotech() and Get.diseases() to find information on any diseases or pathologies associated with your gene ####
GetPathology_Biotech(accession_numbers)
#NA on all counts
summary(AccessionNumbersGO)

Get.diseases(accession_numbers)
#Error


#11. We are going to access structural information using the protti package ####

viewtibble <- fetch_uniprot(accession_numbers)
View(viewtibble)

#12. Pull any available structural information from the Protein DataBase
fetch_pdb("1ZMR")
#pdb_ids auth_asym_id label_asym_id reference_database_accession protein_name 
#<chr>   <chr>        <chr>       <chr>       <chr>        
#1ZMR    A            A           P0A799   Phosphoglyce…
fetch_pdb("2HWG")
#pdb_ids auth_asym_id label_asym_id reference_database_accession protein_name 
#<chr>   <chr>        <chr>           <chr>                        <chr>        
#1 2HWG    A            A             P08839                       Phosphoenolp…
#2 2HWG    A            A             P08839                       Phosphoenolp…
#3 2HWG    B            B             P08839                       Phosphoenolp…
#4 2HWG    B            B             P08839                       Phosphoenolp…


#13. Get information on any available 3D structures for your gene
fetch_alphafold_prediction(accession_numbers)
#Is in file Bioinformatics\Output\Homework_06








