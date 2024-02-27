# load in all of the libraries that you might need
# this should always be at the start of your script
library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)

# set the working directory to the folder containing all of your scripts and data
# filepaths and files should always be in quotes. Variables in R should not.
setwd("/Users/ojohnson/Documents/GitHub/Bioinformatics_Spring2024")

# read in albatross Cytochrome B sequences 
# assign each one to a variable
# Note that these fasta files are contained in a folder called 'Diomedea_exulans'
seq_1 <- readDNAStringSet("Diomedea_exulans/sequence_1.fasta")
seq_2 <- readDNAStringSet("Diomedea_exulans/sequence_2.fasta")
seq_3 <- readDNAStringSet("Diomedea_exulans/sequence_3.fasta")
seq_4 <- readDNAStringSet("Diomedea_exulans/sequence_4.fasta")
seq_5 <- readDNAStringSet("Diomedea_exulans/sequence_5.fasta")
seq_6 <- readDNAStringSet("Diomedea_exulans/sequence_6.fasta")
seq_7 <- readDNAStringSet("Diomedea_exulans/sequence_7.fasta")
seq_8 <- readDNAStringSet("Diomedea_exulans/sequence_8.fasta")

# combine samples into a single variable using the combine ('c') function
seqs <- c(seq_1, seq_2, seq_3, seq_4, seq_5, seq_6, seq_7, seq_8)

# The default GenBank names are very very long, so lets
# rename the samples to something shorter and more legible
# we do this by assigning a list of characters (using the same 'c' function)
# to the 'names' of the combined seqs variable
# check what these names are by first running just the names() function
names(seqs)
names(seqs) <- c("sanfordi_U48946.1", "chionoptera_AF076048.1", "epomophora_AF076049.1",
                 "gibsoni_AF076050.1", "antipodensis_MH330008.1", "antipodensis_MH330010.1",
                 "exulans_U48947.1", "amsterdamensis_U48948.1")

# run the MSA! Assign it to a new variable
albatrossAln <- msa(seqs)

# check the alignment length, two different ways
nchar(albatrossAln)
print(albatrossAln, show="complete") # here, you can also calculate the number of gaps by hand, or use the next step

# Calculate the GC content. First, calculate the frequency of each nucleotide
alFreq <- alphabetFrequency(albatrossAln)
alFreq # here, it also gives you the number of dashes (-) in the alignment, which is the number of gaps

# now pull out the total number of G's and C's
# here, the 'sum' function takes the sum, as you might expect
# the square brackets are for accessing rows and columns of a matrix
# values before the comma access rows, those after the comma access columns
# we want the columns
GC <- sum(alFreq[,"C"]) + sum(alFreq[,"G"]) 
AT <- sum(alFreq[,"A"]) + sum(alFreq[,"T"]) 
# and calculate the percentage that are G or C (out of the total nucleotides)
GC / (GC + AT )


# calculate the GC content a different way, using the 'GC' function in the seqinr package
# we can only run this on one sample at a time, so let's read in the first sample using the 
# read.fasta() function in the seqinr package
# note that because there are multiple 'read.fasta()' functions, we need to specify that
# want to use the one in the seqinr package using the double colon
seq_1.seqinr <- seqinr::read.fasta("Diomedea_exulans/sequence_1.fasta")
# then select the sequence data from the variable. Use the '$' to access it 
seqinr::GC(seq_1.seqinr$U48946.1)

# calculate the identity matrix
# first, convert the alignment to the seqinr format using msaConvert
# because the dist.alignment() function is part of the seqinr package
albatrossAln2 <- msaConvert(albatrossAln, type="seqinr::alignment")
d <- dist.alignment(albatrossAln2, "identity")
d
# this is a fancy way to compare only my 'epomophora' sample to the other samples in the matrix
# and convert the numbers to a percentage 
100 - (round(as.matrix(d)[, "epomophora_AF076049.1", drop=FALSE], digits = 2) * 100)

# translate one sample to an amino acid sequence
# we again need to specify which package to use because the 'translate()' function 
# exists in both the Biostrings and seqinr packages
seq_1_AA <- Biostrings::translate(seq_1)
print(seq_1_AA)


# write the alignment to a fasta file (harder than I expected to figure this out) 
# there is a write function in the phangorn package, but not one that I could find in seqinr or Biostrings
# Biostrings has a write function, but not for fasta-formatted files
albatrossAln_phyDat <- msaConvert(albatrossAln, type="phangorn::phyDat")
write.phyDat(albatrossAln_phyDat, "Diomedea_exulans/albatross_alignment.fasta", format = "fasta")


