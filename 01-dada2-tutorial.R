library(dada2); packageVersion('dada2')

# THIS PIPELINE BEGINS WITH DEMULTIPLEXED FASTQ FILES
# so prior to this step you need to run Ramon's demultiplexer for dada2 
#-------------------------------------------------------------------------
#demultiplexer_for_DADA2
#If you want to use DADA2 with libraries created by adding adapters through ligation and with two levels of adapters - 
#see banzai (github.com/jimmyodonnell/banzai) - you need to demultiplex samples before pairing both ends.
#This script will look for your adapters and pcr primers on your reads, 
#and return 4 fastq files per unique sample: Fwd.1, Fwd.2, Rev.1 and Rev.2.
#---------------------------------------------------------------------------

#BEGIN dada2 tutorisl
#---------------------------------------------------------------------------

# path to the Miseq data files
path <- "./MiSeq_SOP"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#In gray-scale is a heat map of the frequency of each quality score at each base position.
#The mean quality score at each position is shown by the green line'''

# Place filtered files in filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

#learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

