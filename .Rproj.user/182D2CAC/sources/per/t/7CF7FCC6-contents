# load in all of the libraries that you might need
# this should always be at the start of your script
library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)
install.packages("beast")

# set the working directory to the folder containing all of your scripts and data
# filepaths and files should always be in quotes. Variables in R should not.
setwd("/Users/18638/OneDrive/Desktop/GitHub/Bioinformatics/Data/Final_Project/")
#----

# read in albatross Cytochrome B sequences 
# assign each one to a variable
# Note that these fasta files are contained in a folder called 'Diomedea_exulans'
seq_1 <- readDNAStringSet("sequence_539.fasta")
seq_2 <- readDNAStringSet("sequence_540.fasta")
seq_3 <- readDNAStringSet("sequence_544.fasta")
seq_4 <- readDNAStringSet("sequence_901.fasta")
seq_5 <- readDNAStringSet("sequence_902.fasta")
seq_6 <- readDNAStringSet("sequence_903.fasta")
seq_7 <- readDNAStringSet("sequence_904.fasta")
seq_8 <- readDNAStringSet("sequence_905.fasta")
seq_9 <- readDNAStringSet("sequence_906.fasta")
seq_10 <- readDNAStringSet("sequence_907.fasta")
seq_11 <- readDNAStringSet("sequence_908.fasta")
seq_12 <- readDNAStringSet("sequence_909.fasta")
seq_13 <- readDNAStringSet("sequence_910.fasta")
seq_14 <- readDNAStringSet("sequence_911.fasta")
seq_15 <- readDNAStringSet("sequence_918.fasta")
outgroup <-readDNAStringSet("sequence_outgroup.fasta") 
# ----

# combine samples into a single variable using the combine ('c') function
seqs <- c(seq_1, seq_2, seq_3, seq_4, seq_5, seq_6, seq_7, seq_8, seq_9, seq_10, seq_11, seq_12, seq_13, seq_14, seq_15, outgroup)
#----

# rename the samples to something shorter and more legible
# we do this by assigning a list of characters (using the same 'c' function)
# to the 'names' of the combined seqs variable
# check what these names are by first running just the names() function
names(seqs) <- c("Col_5_539", "Col_5_540", "Col_5_544",
                 "Col_5_901", "Col_5_902", "Col_5_903",
                 "Col_5_904", "Col_5_905", "Col_5_906",
                 "Col_5_907", "Col_5_908", "Col_5_909", 
                 "Col_5_910", "Col_5_911", "Col_5_918",
                 "Outgroup")
#----

# run the MSA! Assign it to a new variable
MyotisAln <- msa(seqs)
#----

# check the alignment length, two different ways
nchar(MyotisAln)
print(MyotisAln, show="complete") # here, you can also calculate the number of gaps by hand, or use the next step
#----

# Calculate the GC content. First, calculate the frequency of each nucleotide
alFreq <- alphabetFrequency(MyotisAln)
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
#----

# calculate the identity matrix
# first, convert the alignment to the seqinr format using msaConvert
# because the dist.alignment() function is part of the seqinr package
MyotisAln2 <- msaConvert(MyotisAln, type="seqinr::alignment")
d <- dist.alignment(MyotisAln2, "identity")
d
#----

# translate one sample to an amino acid sequence
# we again need to specify which package to use because the 'translate()' function 
# exists in both the Biostrings and seqinr packages
seq_1_AA <- Biostrings::translate(seq_10)
print(seq_1_AA)
#----

# write the alignment to a fasta file
# there is a write function in the phangorn package, but not one that I could find in seqinr or Biostrings
# Biostrings has a write function, but not for fasta-formatted files
MyotisAln_phyDat <- msaConvert(MyotisAln, type="phangorn::phyDat")
write.phyDat(MyotisAln_phyDat, "MyotisAln.fasta", format = "fasta")
#----

#Neighbor Joining Phylogenetic Tree
MyotisAlnTreeNJ <- nj(d)
plot(MyotisAlnTreeNJ, main="Phylogenetic Tree of Myotis myotis Gene Sequences")

#Maximum Likelihood Phylogenetic Tree
library(beast)
?beast
run_mcmc <- beast(myDataList = d, subsetIndex = seqs, 
                  zeroNormalization = TRUE)

print(run_mcmc)
plot(run_mcmc, fileName = "beast_plot.pdf", timeScale=1/6, xlab = "hours", ylab = "growth")
