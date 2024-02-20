
library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(UniprotR)
library(protti)

# install.packages("UniprotR")
# install.packages("protti")
# install.packages("DelayedArray")


mySequenceFile <- system.file("Horse/sequence-1.fasta", package="msa")

seq_1 <- readDNAStringSet("Horse/sequence-1.fasta")
seq_2 <- readDNAStringSet("Horse/sequence-2.fasta")
seq_3 <- readDNAStringSet("Horse/sequence-3.fasta")
seq_4 <- readDNAStringSet("Horse/sequence-4.fasta")
seq_5 <- readDNAStringSet("Horse/sequence-5.fasta")

MySequences <- readDNAStringSet("Horse/sequence-1.fasta")

# ClustalW

NewSequence <- readDNAStringSet("FullSequence.fasta")

myFirstAlignment <- msa(NewSequence)

myFirstAlignment

print(myFirstAlignment, show="complete")

print(NewSequence, show="complete")

DNA_Sequence <- DNAString(NewSequence)

Alignment <- DNAStringSet(c(NewSequence))

alignment_freq <- letterFrequency(Alignment)

?letteralignment
?letterFrequency

GC_content <- alphabetFrequency(myFirstAlignment)

hemoseq <- readDNAStringSet("FullSequence.fasta")
names(hemoseq) <- c("SF23-1", "SF124-2", "SF120-3", "SF115-4", "SF55-5")

hemoAln <- msa(hemoseq)

hemoAln

hemoAln2 <- msaConvert(hemoAln, type="seqinr::alignment")
hemoAln2 <- msaConvert(hemoAln, type = "seqinr::alignment")

d <- dist.alignment(hemoAln2, "identity")

d

min(d)
max(d)

?translate

translate(seq_1)

DNA_string <- DNAString(seq_2)

seq_2

AA_string

AA_string <- Biostrings::translate(seq_2)

# install.packages("phangorn")
library(phangorn)

Alignment_phyDat <- msaConvert(myFirstAlignment, type="phangorn::phyDat")

write.phyDat(Alignment_phyDat, "FUllSequence.fasta", format = "fasta")

# BiocManager::install("GenomicAlignments")

?writeXStringSet
Biostrings::writeXStringSet("Horse/sequence_2.fasta")
writeXStringSet(AA_string, "Horse/sequence_2_AminoAcid.fasta", append=FALSE)

# A0A5C6NVL2
# A0A7K5Q6Y4
# A0A9N7VP42
# A0A099YVD3
# A0A811Z0W7

ReadAccessionNumbers <- read.csv("UniprotRAccessionNumbers/AccessionNumbers.txt")
ReadAccessionNumbers
# change value after the dollar sign to the column name:
ReadAccessionNumbers$AccessionNumbers

# this would read all AccessionNumbers
GeneOntologyTerms <- GetProteinGOInfo(ReadAccessionNumbers$AccessionNumbers)

# adding the [1] reads just the first item, just testing to see if it works, and it does
GeneOntologyTerms <- GetProteinGOInfo(ReadAccessionNumbers$AccessionNumbers[1])
# GeneOntologyTerms <- GetProteinGOInfo(read.csv("UniprotRAccessionNumbers/AccessionNumbers.txt"))
# GeneOntologyTerms <- GetProteinGOInfo(Test$AccessionNumbers) # not sure what the Test variable is

?GetProteinGOInfo

# "P0A799"

# GeneOntologyTerms <- GetProteinGOInfo("P0A799")
# one of the three columns returned no results
GeneOntologyTerms
# so we if replace it with a character value of "NA", it will plot
GeneOntologyTerms[1,4] <- "NA"
# we can also add a directory path to save the file to your computer:
PlotGoInfo(GeneOntologyTerms, directorypath = "GoInfoPlot")

class(Test)

# test reads the information vertically

GeneOntologyTerms <- GetProteinGOInfo("A0A5FC6NVL2")

PlotGoInfo(AC1)

#Tried up to here with GO information with no success, moving on

Fetch <- fetch_uniprot("P0A799")

View(Fetch)

#1ZMR is our code
# this accession number had issues for some people. not sure why. The code looks correct
fetch_alphafold_prediction("1ZMR")