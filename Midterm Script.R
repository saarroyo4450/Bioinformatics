# I will first begin by loading the packages necessary for this research lab.

library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(GenomicAlignments)
library(UniprotR)
library(protti)
library(dplyr)

# Now that I have installed the necessary packages I can first begin by 
# creating a variable that will read the 20 human sequences that were 
# provided for me.

# duplicate line. Used same readDNAStringSet function to read in sequences below
# My20Sequences <- readDNAStringSet("sequences.fasta")

# Now we can use this variable to align our 20 sequences in R.

# The20SequencesAligned <- msa(My20Sequences)

# We have successfully aligned our 20 sequences. lets see what we have done 
# with the print command. This will allow me to check for anything unusual 
# between the sequences.

# Before I do this, I am going to compare the sequences by converting them
# to seqinr format and computing a distance matrix. 

PatientsInSeqinr <- readDNAStringSet("sequences.fasta")
# generally not a good idea to name things as a number, because R considers it 
# to be a number instead of a character. Something like "H1", "H2" would be a good alternative
# names(PatientsInSeqinr) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
PatientsInSeqinrAligned <- msa(PatientsInSeqinr)
PatientsInSeqinrAligned
PatientsInSeqinr2 <- msaConvert(PatientsInSeqinrAligned, type = "seqinr::alignment")
PatientsInSeqinr2
DistanceMatrix <- dist.alignment(PatientsInSeqinr2, "identity")
DistanceMatrix

# Based on the distance matrix, the sequences for patient 4, patient 6,
# and patient 10, have differences. 

# Now I will use the print command in order to look for these differences

print(The20SequencesAligned, show="complete")

# Patient 4:At column 39, the patient had an "A" instead of a "C" as the 
# consensus sequence shows suggesting a point mutation.

# Patient 10:At column 39, the patient had a "G" instead of a "C" as the 
# consensus sequence shows suggesting a point mutation.

# Patient 6:there is a gap in alignment 1 on our 20th patient. For the
# first 4 pairs the consensus sequence shows, "AACT" but patient 6 
# shows "-AAT" they are missing C signaling a deletion. At column 47, the 
# patient had an "A" instead of a "G" as the consensus sequence shows 
# suggesting a point mutation.

# I exported my sequence and ran it through BLAST in order to look for 
# similar sequences and discovered that the most identical gene was called,
# "hbb gene for beta globin, partial cds" with the accession number, 
# "LC121775.1" 

# I will now translate the most different patient sequence to a protein.
# Based on my descriptions of the differences between the 3 unique patient
# sequences, it is obviously patient 6 with the most different sequence.
# For argument's sake of course, I will prove this with 2 simple commands

max(DistanceMatrix)
min(DistanceMatrix)

# a more programmatic way to do this: 
# convert from distance to matrix so that we can access row/column names
# then to a dataframe to use the dplyr package
df.dist <- as.data.frame(as.matrix(DistanceMatrix, labels=TRUE))
# filter the row containing the maximum value
dplyr::filter(df.dist, Homo_sapiens_6 == max(df.dist))

# or we could pull out which cell in the matrix has the maximum value:
which(DistanceMatrix == max(DistanceMatrix))

# The value "0.1184929" matches that of patient 6 based on the Distance Matrix

# Now I shall translate patient 6's sequence to a protein by first
# using the translate command on our MSA. 

Translated20Sequences <- Biostrings::translate(My20Sequences)
Translated20Sequences

# But we need our most different sequence, which is patient 6 so lets single that out

# align <- readDNAStringSet("/Users/sergi/Documents/GitHub/Bioinformatics/sequences.fasta")
# class(align)
Sequence6DNA <- PatientsInSeqinr$Homo_sapiens_6 # no need to re-read data, already exists

# Now we need to translate this sequence.

Sequence6AA <- Biostrings::translate(Sequence6DNA)
Sequence6AA

# or:
Sequence6AA <- Translated20Sequences$Homo_sapiens_6

# Finally we must put it into a fasta file

# write.fasta(Biostrings::translate(Sequence6DNA),"Homo_sapiens_6","/Users/sergi/Documents/GitHub/Bioinformatics/TranslatePatient6.fasta",open = "w",nbchar = 60,as.string = FALSE)

# a bit shorter way to write this. Your sample was already translated, to that's repetitive
# and the working directory is already the Bioinformatics folder, so you can just 
# use the file name
# if a line of code is very long, you can wrap it onto the next line right after any comma
write.fasta(Sequence6AA, "Homo_sapiens_6", "TranslatePatient6.fasta", 
            open = "w", nbchar = 60, as.string = FALSE)

# I ran this through Uniprot and discovered the accession number closest to 
# this sequence is A0A0J9YWK4, associated with the HBB gene. 

# Associated diseases found were blood related diseases such as malaria, 
# sickle cell disease, beta-thalassemia, and heinz body anemia. There was
# no evidence to suggest the patient had any of these associated diseases.

# I was able to obtain a picture of the structure of the gene through 
# Alpha Fold, it is labeled Patient6Protein.
