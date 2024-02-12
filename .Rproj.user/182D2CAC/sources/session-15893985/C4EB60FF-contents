if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
install.packages("seqinr")
library(seqinr)
library(ape)

setwd("Data/Homework05/")
#Set Working directory for Data and Homework folder

getwd()
#Checked that directory was in place

mySequences01 <- readDNAStringSet("sequence01.fasta")
mySequences02 <- readDNAStringSet("sequence02.fasta")
mySequences03 <- readDNAStringSet("sequence03.fasta")
mySequences04 <- readDNAStringSet("sequence04.fasta")
mySequences05 <- readDNAStringSet("sequence05.fasta")
#Read each sequence and gave each a name

mySequenceFile <- c(mySequences01, mySequences02, mySequences03, mySequences04, mySequences05)
#Combined each file into one to complete msa

mySequenceFile
#printed sequence file to check that combine worked

print(mySequenceFile, show="complete")
#Print the alignment with show="complete"



####Run clustal ####
clustalw_alignment <- msa(mySequenceFile, "ClustalW")
print(clustalw_alignment)




####Function to count gaps in a sequence ####
count_gaps <- function(sequence) {
#Search for the presence of "-" in the sequence
  gap_positions <- grepl("-", sequence)
  
#Count the number of TRUE values (i.e., gaps)
  num_gaps <- sum(gap_positions)
  
  return(num_gaps)
}

#Count gaps in each sequence of the alignment
num_gaps_in_sequences <- sapply(mySequenceFile, count_gaps)

#Total gaps in the alignment
total_gaps <- sum(num_gaps_in_sequences)

#Print the total number of gaps
print(total_gaps)
#0 gaps total in consensus



#####Calculate the width of the alignments ####
alignment_length <- width(mySequenceFile)

#Print the length of the alignments
print(alignment_length)
#[1] 1950 1950 1950 1497 1950



####Calculate GC content ####
# Convert the alignment to a DNAStringSet object
alignment <- DNAStringSet(mySequenceFile)

# Calculate the GC content for each position in the alignment
gc_content <- vcountPattern("GC", alignment) / width(alignment)

# Print the GC content for each position
print(gc_content)
#[1] 0.05179487 0.05230769 0.05230769 0.04542418 0.05179487
#    5.18%      5.23%      5.23%      4.54%      5.18%





####Convert alignment to SeqinR format ####
msa(mySequenceFile)





# Example MSA object (replace this with your actual MSA object)






HyenaEnv2 <- msaConvert(mySequenceFile, type = "seqinr::alignment")


myAlignment <- msa(mySequenceFile)
converted_alignment <- msaConvert(myAlignment, type = "seqinr::alignment")
# Compute a distance matrix using the 'dist.alignment' function
distance_matrix <- dist.alignment(myAlignment)



# Print the distance matrix
print(distance_matrix)

















# CLEAN UP #####
# Clear environment
rm(list = ls()) 
# Clear packages
# requires the package pacman to work# CLEAN UP #####
# Clear environment
p_unload(all)  # Remove all add-ons
# Clear console
cat("\014")  # ctrl+L
# Clear mind :)