# Script to identify potential G-quadruplex in an input fasta file.
#
# It uses R/BioCondutor and the pqsfinder package with default
# parameters.
#
# PQS means: Potential Quadruplex Sequences...
#
# The package can be found here:
# https://bioconductor.org/packages/release/bioc/html/pqsfinder.html
# Here is the userguide:
# https://bioconductor.org/packages/release/bioc/vignettes/pqsfinder/inst/doc/pqsfinder.html
# And here is the article to cite:
# https://doi.org/10.1093/bioinformatics/btx413

#-----------------------------------
# TO RUN THIIS SCRIPT:
#-----------------------------------
#
# - Be sure to have the needed R/BioConductor LIBRARIES installed
#   (Biostrings, pqsfinder, GenomicRanges).
#
# - After having changed the PATHS, FOLDER AND FILENAMES according
#   to the analysis to be run you can launch this script using
#   the following command in the bash:
#
#   R CMD BATCH --no-save identify_pqs.R
#
#-----------------------------------

#-----------------------------------
# SET THE PATHS, FOLDER AND FILENAMES
#-----------------------------------
# Path to working directory, where to save the output...
wd = '/home/remo/Dropbox/WORK/PROJECTS/federico/damiano/gquad'
# Path and name of the fasta input to be analyzed.
input = '/home/remo/data/genomes/mmusculus/Mus_musculus.GRCm38.dna.primary_assembly.fa'
# Output table to write with the results.
output = 'mmusculus_grcm38_pqs.txt'
#-----------------------------------

#-----------------------------------
# LOAD THE LIBRARIES
#-----------------------------------
# Library to read the fasta.
library(Biostrings)
# Library to search for G-quadruplex.
library(pqsfinder)
# Library to work with ranges.
library(GenomicRanges)
#-----------------------------------

#-----------------------------------
# LET'S THE CODE BEGIN!
#-----------------------------------
# Set and go in the working directory.
setwd(wd)

# Load the fasta.
fasta = readDNAStringSet(input)

# Workaround needed to extract only the id of the sequences from the fasta 
# in case also a description is present....
# Extract the id from the fasta object.
ids = names(fasta)
# Split according to the space.
ids = strsplit(ids,' ')
# Extract only the first element from the split.
ids = sapply(ids,function(x)x[1])

# Get the number of sequences in the fasta.
nseq = length(fasta)

# Create the variable in which to write the output.
t = c()

# Main cycle to iterate the PQS search on each sequence and organize
# the results in a proper structure to be exported in a txt table.
for(i in 1:nseq) {
  seq = fasta[[i]]
  id = ids[i]
  print(paste0('Analyzing sequence ',id))
  pqs = pqsfinder(seq)
  if(length(pqs) > 0) {
    res = as(pqs,"GRanges")
    res = as.data.frame(res)
    res$seqnames = id
    res = res[,1:9]
    res$seq = as.data.frame(pqs)$seq
    t = rbind(t,res)
  }
}

# Dimension of the generated table with the results.
dim(t)

# Write the output file with the results.
write.table(t,file=output,sep='\t',row.names=F,quote=F)

# This is the End...
print("Done!")
