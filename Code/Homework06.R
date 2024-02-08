if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("GenomicAlignments")
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
install.packages("seqinr")
library(seqinr)
library(ape)
install.packages(UniprotR)
install.packages(protti)
install.packages(r3dmol)


setwd("Data/Homework06/")
#Set Working directory for Data and Homework folder ####

getwd()
#Checked that directory was in place ####

mySequences01 <- readDNAStringSet("sequence01.fasta")

