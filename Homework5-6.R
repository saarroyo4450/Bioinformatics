
library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(UniprotR)


mySequenceFile <- system.file("Horse/sequence-1.fasta", package="msa")

seq_1 <- readDNAStringSet("Horse/sequence-1.fasta")
seq_2 <- readDNAStringSet("Horse/sequence-2.fasta")
seq_3 <- readDNAStringSet("Horse/sequence-3.fasta")
seq_4 <- readDNAStringSet("Horse/sequence-4.fasta")
seq_5 <- readDNAStringSet("Horse/sequence-5.fasta")

MySequences <- readDNAStringSet("Horse/sequence-1.fasta")

ClustalW

NewSequence <- readDNAStringSet("FullSequence.fasta")

myFirstAlignment <- msa(NewSequence)

myFirstAlignment

print(myFirstAlignment, show="complete")

print(NewSequence, show="complete")

DNA_Sequence <- DNAString(NewSequence)

Aignment <- DNAStringSet(c(NewSequence))

alignment_freq <- letterFrequency(alignment)

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

DNA_string <- DNAString(seq_1)

AA_string <- Biostrings::translate(seq_2)

install.packages("phangorn")

Alignment_phyDat <- msaConvert(myFirstAlignment, type="phangorn::phyDat")

write.phyDat(Alignment_phyDat, "FUllSequence.fasta", format = "fasta")

BiocManager::install("GenomicAlignments")

install.packages("UniprotR")
install.packages("protti")
install.packages("DelayedArray")

?writeXStringSet
Biostrings::writeXStringSet("Horse/sequence_2.fasta")
writeXStringSet(AA_string, "Horse/sequence_2_AminoAcid.fasta", append=FALSE)

A0A5C6NVL2
A0A7K5Q6Y4
A0A9N7VP42
A0A099YVD3
A0A811Z0W7

ReadAccessionNumbers <- read.csv("UniprotRAccessionNumbers/AccessionNumbers.txt")
ReadAccessionNumbers
ReadAccessionNumbers$A0A5C6NVL2

GeneOntologyTerms <- GetProteinGOInfo(ReadAccessionNumbers)