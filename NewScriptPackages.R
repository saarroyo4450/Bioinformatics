if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
n

BiocManager::install("Biostrings")
install.packages("seqinr") 

library(Biostrings)
