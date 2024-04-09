# load in all of the libraries that you might need
# this should always be at the start of your script
library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)
library(ape)

# set the working directory to the folder containing all of your scripts and data
# filepaths and files should always be in quotes. Variables in R should not.
setwd("/Users/18638/OneDrive/Desktop/GitHub/Bioinformatics/Data/Final_Project/")
#----

# read in Myotis myotis D-loop sequences ####
# assign each one to a variable
# Note that these fasta files are contained in a folder called "Data/Final_Project"

seq_1 <- readDNAStringSet("sequence_col4_501.fasta")
seq_2 <- readDNAStringSet("sequence_col4_502.fasta")
seq_3 <- readDNAStringSet("sequence_col4_503.fasta")
seq_4 <- readDNAStringSet("sequence_col4_504.fasta")
seq_5 <- readDNAStringSet("sequence_col4_506.fasta")
seq_6 <- readDNAStringSet("sequence_col4_507.fasta")
seq_7 <- readDNAStringSet("sequence_col4_508.fasta")
seq_8 <- readDNAStringSet("sequence_col4_509.fasta")
seq_9 <- readDNAStringSet("sequence_col4_510.fasta")
seq_10 <- readDNAStringSet("sequence_col4_511.fasta")
seq_11 <- readDNAStringSet("sequence_col4_512.fasta")
seq_12 <- readDNAStringSet("sequence_col4_513.fasta")
seq_13 <- readDNAStringSet("sequence_col4_514.fasta")
seq_14 <- readDNAStringSet("sequence_col4_515.fasta")
seq_15 <- readDNAStringSet("sequence_col4_516.fasta")
seq_16 <- readDNAStringSet("sequence_col4_518.fasta")
seq_17 <- readDNAStringSet("sequence_col4_519.fasta") 
seq_18 <- readDNAStringSet("sequence_col4_520.fasta") 
seq_19 <- readDNAStringSet("sequence_col4_521.fasta") 
seq_20 <- readDNAStringSet("sequence_col4_522.fasta") 
seq_21 <- readDNAStringSet("sequence_col4_523.fasta") 
seq_22 <- readDNAStringSet("sequence_col4_524.fasta") 
seq_23 <- readDNAStringSet("sequence_col4_525.fasta") 
seq_24 <- readDNAStringSet("sequence_col4_526.fasta") 
seq_25 <- readDNAStringSet("sequence_col4_527.fasta") 
seq_26 <- readDNAStringSet("sequence_col4_529.fasta") 
seq_27 <- readDNAStringSet("sequence_col4_530.fasta")
seq_28 <- readDNAStringSet("sequence_col4_531.fasta")
seq_29 <- readDNAStringSet("sequence_col4_532.fasta")
seq_30 <- readDNAStringSet("sequence_col4_901.fasta")
seq_31 <- readDNAStringSet("sequence_col4_905.fasta")
seq_32 <- readDNAStringSet("sequence_col4_906.fasta")
seq_33 <- readDNAStringSet("sequence_col4_908.fasta")
seq_34 <- readDNAStringSet("sequence_col4_917.fasta")
seq_35 <- readDNAStringSet("sequence_col4_918.fasta")
seq_36 <- readDNAStringSet("sequence_col4_922.fasta")
seq_37 <- readDNAStringSet("sequence_col4_923.fasta")
seq_38 <- readDNAStringSet("sequence_col4_924.fasta")
seq_39 <- readDNAStringSet("sequence_col4_926.fasta")
seq_40 <- readDNAStringSet("sequence_col4_927.fasta")
seq_41 <- readDNAStringSet("sequence_col4_929.fasta")
seq_42 <- readDNAStringSet("sequence_col5_501.fasta")
seq_43 <- readDNAStringSet("sequence_col5_502.fasta")
seq_44 <- readDNAStringSet("sequence_col5_503.fasta")
seq_45 <- readDNAStringSet("sequence_col5_504.fasta")
seq_46 <- readDNAStringSet("sequence_col5_506.fasta")
seq_47 <- readDNAStringSet("sequence_col5_507.fasta")
seq_48 <- readDNAStringSet("sequence_col5_508.fasta")
seq_49 <- readDNAStringSet("sequence_col5_509.fasta")
seq_50 <- readDNAStringSet("sequence_col5_510.fasta")
seq_51 <- readDNAStringSet("sequence_col5_512.fasta")
seq_52 <- readDNAStringSet("sequence_col5_519.fasta")
seq_53 <- readDNAStringSet("sequence_col5_523.fasta")
seq_54 <- readDNAStringSet("sequence_col5_524.fasta")
seq_55 <- readDNAStringSet("sequence_col5_525.fasta")
seq_56 <- readDNAStringSet("sequence_col5_526.fasta")
seq_57 <- readDNAStringSet("sequence_col5_527.fasta")
seq_58 <- readDNAStringSet("sequence_col5_528.fasta")
seq_59 <- readDNAStringSet("sequence_col5_529.fasta")
seq_60 <- readDNAStringSet("sequence_col5_530.fasta")
seq_61 <- readDNAStringSet("sequence_col5_531.fasta")
seq_62 <- readDNAStringSet("sequence_col5_532.fasta")
seq_63 <- readDNAStringSet("sequence_col5_533.fasta")
seq_64 <- readDNAStringSet("sequence_col5_534.fasta")
seq_65 <- readDNAStringSet("sequence_col5_535.fasta")
seq_66 <- readDNAStringSet("sequence_col5_536.fasta")
seq_67 <- readDNAStringSet("sequence_col5_537.fasta")
seq_68 <- readDNAStringSet("sequence_col5_538.fasta")
seq_69 <- readDNAStringSet("sequence_col5_539.fasta")
seq_70 <- readDNAStringSet("sequence_col5_540.fasta")
seq_71 <- readDNAStringSet("sequence_col5_544.fasta")
seq_72 <- readDNAStringSet("sequence_col5_901.fasta")
seq_73 <- readDNAStringSet("sequence_col5_902.fasta")
seq_74 <- readDNAStringSet("sequence_col5_903.fasta")
seq_75 <- readDNAStringSet("sequence_col5_904.fasta")
seq_76 <- readDNAStringSet("sequence_col5_905.fasta")
seq_77 <- readDNAStringSet("sequence_col5_906.fasta")
seq_78 <- readDNAStringSet("sequence_col5_907.fasta")
seq_79 <- readDNAStringSet("sequence_col5_908.fasta")
seq_80 <- readDNAStringSet("sequence_col5_909.fasta")
seq_81 <- readDNAStringSet("sequence_col5_910.fasta")
seq_82 <- readDNAStringSet("sequence_col5_918.fasta")
seq_83 <- readDNAStringSet("sequence_col3_505.fasta")
seq_84 <- readDNAStringSet("sequence_col3_506.fasta")
seq_85 <- readDNAStringSet("sequence_col3_507.fasta")
seq_86 <- readDNAStringSet("sequence_col3_508.fasta")
seq_87 <- readDNAStringSet("sequence_col3_509.fasta")
seq_88 <- readDNAStringSet("sequence_col3_903.fasta")
seq_89 <- readDNAStringSet("sequence_col3_905.fasta")
seq_90 <- readDNAStringSet("sequence_col3_907.fasta")
seq_91 <- readDNAStringSet("sequence_col3_501.fasta")
seq_92 <- readDNAStringSet("sequence_col3_502.fasta")
seq_93 <- readDNAStringSet("sequence_col3_503.fasta")
seq_94 <- readDNAStringSet("sequence_col3_504.fasta")
seq_95 <- readDNAStringSet("sequence_col4_517.fasta")
seq_96 <- readDNAStringSet("sequence_col5_911.fasta")
seq_outgroup01 <- readDNAStringSet("sequence_J60.fasta")
seq_outgroup02 <- readDNAStringSet("sequence_J61.fasta")
seq_outgroup03 <- readDNAStringSet("sequence_J62.fasta")
seq_outgroup04 <- readDNAStringSet("sequence_J63.fasta")
seq_outgroup05 <- readDNAStringSet("sequence_J64.fasta")
seq_outgroup06 <- readDNAStringSet("sequence_J65.fasta")
# ----

# combine samples into a single variable using the combine ('c') function
seqs <- c(seq_1, seq_2, seq_3, seq_4, seq_5, seq_6, seq_7, seq_8, seq_9, seq_10, seq_11, seq_12, seq_13, seq_14, seq_15, seq_16, seq_17, seq_18, seq_19, seq_20, seq_21, seq_22, seq_23, seq_24, seq_25, seq_26, seq_27, seq_28, seq_29, seq_30, seq_31, seq_32, seq_33, seq_34, seq_35, seq_36, seq_37, seq_38, seq_39, seq_40, seq_41, seq_42, seq_43, seq_44, seq_45, seq_46, seq_47, seq_48, seq_49, seq_50, seq_51, seq_52, seq_53, seq_54, seq_55, seq_56, seq_57, seq_58, seq_59, seq_60, seq_61, seq_62, seq_63, seq_64, seq_65, seq_66, seq_67, seq_68, seq_69, seq_70, seq_71, seq_72, seq_73, seq_74, seq_75, seq_76, seq_77, seq_78, seq_79, seq_80, seq_81, seq_82, seq_83, seq_84, seq_85, seq_86, seq_87, seq_88, seq_89, seq_90, seq_91, seq_92, seq_93, seq_94, seq_95, seq_96, seq_outgroup01, seq_outgroup02, seq_outgroup03, seq_outgroup04, seq_outgroup05, seq_outgroup06)
#----

# rename the samples ####
#to something shorter and more legible
# we do this by assigning a list of characters (using the same 'c' function)
# to the 'names' of the combined seqs variable
# check what these names are by first running just the names() function
names(seqs) <- c("col4_501", "col4_502", "col4_503",
                 "col4_504", "col4_506", "col4_507",
                 "col4_508", "col4_509", "col4_510",
                 "col4_511", "col4_512", "col4_513", 
                 "col4_514", "col4_515", "col4_516",
                 "col4_518", "col4_519", "col4_520",
                 "col4_521", "col4_522", "col4_523",
                 "col4_524", "col4_525", "col4_526",
                 "col4_527", "col4_529", "col4_530", 
                 "col4_531", "col4_533", "col4_901",
                 "col4_905", "col4_906", "col4_908", 
                 "col4_917", "col4_918", "col4_922", 
                 "col4_923", "col4_924", "col4_926",
                 "col4_927", "col4_929", "col5_501",
                 "col5_502", "col5_503", "col5_504",
                 "col5_506", "col5_507", "col5_508",
                 "col5_509", "col5_510", "col5_512",
                 "col5_519", "col5_523", "col5_524",
                 "col5_525", "col5_526", "col5_527",
                 "col5_528", "col5_529", "col5_530",
                 "col5_531", "col5_532", "col5_533",
                 "col5_534", "col5_535", "col5_536",
                 "col5_537", "col5_538", "col5_539",
                 "col5_540", "col5_544", "col5_901",
                 "col5_902", "col5_903", "col5_904",
                 "col5_905", "col5_906", "col5_907",
                 "col5_908", "col5_909", "col5_910",
                 "col5_918", "col3_505", "col3_506",
                 "col3_507", "col3_508", "col3_509",
                 "col3_903", "col3_905", "col3_907",
                 "col3_501", "col3_502", "col3_503",
                 "col3_504", "col4_517", "col5_911", 
                 "J60", "J61", "J62", 
                 "J63", "J64", "J65")
#----
#MSA ####
# run the MSA! Assign it to a new variable
MyotisAln <- msa(seqs)
#----

# check the alignment length, two different ways
nchar(MyotisAln)
print(MyotisAln, show="complete") # here, you can also calculate the number of gaps by hand, or use the next step

#---

# calculate the identity matrix
# first, convert the alignment to the seqinr format using msaConvert
# because the dist.alignment() function is part of the seqinr package
MyotisAln2 <- msaConvert(MyotisAln, type="seqinr::alignment")
d <- dist.alignment(MyotisAln2, "identity")
d
#----

#Neighbor Joining Phylogenetic Tree
MyotisAlnTreeNJ <- nj(d)
plot(MyotisAlnTreeNJ, main="Phylogenetic Tree of Myotis myotis Gene Sequences")
#---
#Reroot the NJ tree ####
RerootNJ <- root(MyotisAlnTreeNJ, MyotisAlnTreeNJ$tip.label[99])
plot(RerootNJ)

#---
#Notes
#Look for conflicting topology
#Find package to line up the tips of the beast tree and the nj tree
#Robinson folds distance and KA/KP distance
#Measure if trees have the same topology
#needs to be .tre file
#Read up on Ape package
#Reroot the tree to the same outgroup and there is a command in Ape
#Command in R to drop the tip, give it a vector

#---
#Reading in BEAST tree ####
BEASTTree <- read.nexus("/Users/18638/OneDrive/Desktop/GitHub/Bioinformatics/Data/Final_Project/TreeAnnotator/MyotisTree09042024")

#---

#Reroot BEAST Tree ####
RerootBEAST <- root(BEASTTree, BEASTTree$tip.label[1])
plot(RerootBEAST)

#---

#Dropping tips ####
common_taxa <- intersect(RerootNJ$tip.label, RerootBEAST$tip.label)
common_taxa
#Prune trees to include only common taxa
RerootNJ_pruned <- drop.tip(RerootNJ, setdiff(RerootNJ$tip.label, common_taxa))
RerootBEAST_pruned <- drop.tip(RerootBEAST, setdiff(RerootBEAST$tip.label, common_taxa))


#---
#Robinson-Foulds distance ####
rf_distance <- dist.topo(RerootNJ_pruned, RerootBEAST_pruned)
rf_distance
#rf_distance = 4

#---

#Finding number of bipartitions ####
# Get the number of tips in the tree
num_tips <- length(RerootNJ_pruned$tip.label)

# Calculate the total number of bipartitions
total_bipartitions <- 2^num_tips - 3

total_bipartitions
#61 for NJ Tree

#---
num_tipsbeast <- length(RerootBEAST_pruned$tip.label)

# Calculate the total number of bipartitions
total_bipartitionsbeast <- 2^num_tipsbeast - 3

total_bipartitionsbeast
#61 for BEAST Tree

#---
#Both have the same number of bipartitions 











