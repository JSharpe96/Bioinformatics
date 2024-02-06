if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
install.packages("seqinr")
library(seqinr)
library(ape)
install.packages("phangorn")
library(phangorn)


setwd("Data/Homework05/")
#Set Working directory for Data and Homework folder ####

getwd()
#Checked that directory was in place ####

mySequences01 <- readDNAStringSet("sequence01.fasta")
mySequences02 <- readDNAStringSet("sequence02.fasta")
mySequences03 <- readDNAStringSet("sequence03.fasta")
mySequences04 <- readDNAStringSet("sequence04.fasta")
mySequences05 <- readDNAStringSet("sequence05.fasta")
mycombinedSeq <- readDNAStringSet("mycombinedSeq.txt")

mySequenceFile <- c(mySequences01, mySequences02, mySequences03, mySequences04, mySequences05)
mySequenceFile 

#### Creating a consensus ####
alignment_set <- DNAStringSet(mySequenceFile)
consensus <- consensusString(alignment_set)
print(consensus)


#### Creating an msa alignment ####
myFirstAlignment <- msa(mySequenceFile)
myFirstAlignment


print(myFirstAlignment, show="complete")

####Function to count gaps in a consensus ####
count_gaps <- function(sequence) {
  #Search for the presence of "-" in the sequence
  gap_positions <- grepl("-", sequence)
  
  #Count the number of TRUE values (i.e., gaps)
  num_gaps <- sum(gap_positions)
  
  return(num_gaps)
}

#Count gaps in each sequence of the consensus
num_gaps_in_sequences <- sapply(consensus, count_gaps)

#Total gaps in the alignment
total_gaps <- sum(num_gaps_in_sequences)

#Print the total number of gaps
print(total_gaps)
#1 gap



#####Calculate the width of the alignments ####
alignment_length <- width(consensus)

#Print the length of the alignments
print(alignment_length)
#ALignmnet length is 1950



####Calculate GC content ####
# Convert the alignment to a DNAStringSet object
alignment <- DNAStringSet(consensus)

# Calculate the GC content for each position in the alignment
c_content <- vcountPattern("C", alignment) / width(alignment)
G_content <- vcountPattern("G", alignment) / width(alignment)

# Print the GC content for each position
print(c_content)
print(G_content)
#[1] 0.2758974
#[1] 0.2071795
#GC content is .482 = 48.2%

####Convert alignment to SeqinR format ####
HyenaEnv2 <- msa(mycombinedSeq)
HyenaEnv2Com


HyenaEnv2Com <- msaConvert(HyenaEnv2, type="seqinr::alignment")
# Compute distance matrix using SeqinR functions
# Assuming you want to compute identity distance
d <- dist.alignment(HyenaEnv2Com)

# View the distance matrix
print(d)

####Creat Phylogenetic Tree ####
Env2Tree <- nj(d)
plot(Env2Tree, main="Phylogenetic Tree of HyenaEnv2 Gene Sequences")
#MG805960.1 Hyaena hyaena and MG805961.1 Parahyaena brunnae are the most closely related
#MG805960.1 Hyaena hyaena and MG805959.1 Crocuta crocuta are the least closely related


#### Translate DNA sequence into AA sequence ####
dna_sequences <- readDNAStringSet("Data/Homework05/sequence01.fasta")
amino_acid_sequences <- translate(dna_sequences)
amino_acid_sequences

####Convert alignment to phangorn ####
Alignment_phyDat <- msaConvert(myFirstAlignment, type="phangorn::phyDat")
Alignment_phyDat
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")



