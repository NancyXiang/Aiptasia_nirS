#########################
####   DADA2_nirS    ####
#########################


# run in R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"

setwd("/home/xiangn/R/nirS_full")

save.image("/home/xiangn/R/nirS_full/dada2_nirS_full.Rdata")

load("/home/xiangn/R/nirS_full/dada2_nirS_full.Rdata")


# getting ready
library(dada2)
packageVersion("dada2")
# 1.22.0
library(ShortRead)
packageVersion("ShortRead")
# 1.52.0
library(Biostrings)
packageVersion("Biostrings")
# 2.62.0


path <- "/home/xiangn/DADA2/nirS_full"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)

fnFs <- sort(list.files(path, pattern = "_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fq.gz", full.names = TRUE))

### ### ### ### ### ### ###
###   identify primers  ###
### ### ### ### ### ### ###

FWD <- "CCTAYTGGCCGCCRCART"  ## forward primer sequence nirS
REV <- "TCCMAGCCRCCRTCRTGCAG"  ## reverse primer sequence nirS
## nirS-1F, qR;
## REF: Lee, J. A., & Francis, C. A. (2017). DeepnirSamplicon sequencing of San Francisco Bay sediments enables prediction of geography and environmental conditions from denitrifying community composition. Environmental microbiology, 19(12), 4897-4912. doi:10.1111/1462-2920.13920


allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}


FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN sub-directory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))


filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Remove abiguous bases (Ns) in the sequencing reads


primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


### ### ### ### ### ### ###
###   remove primers    ###
### ### ### ### ### ### ###

# first-step cutadapt

cutadapt <- "/home/xiangn/miniconda3/envs/DADA2/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")
# cutadapt 3.5 with Python 3.9.9


path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)


# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)


# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


# sanity check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
# NOTE: still have primers in unexpected locations due to orientation mixups


# second-step cutadapt

# for orientation mixups
FWD2 <- FWD.orients[["RevComp"]]
REV2 <- REV.orients[["RevComp"]]


path.cut2 <- file.path(path, "cutadapt2")
if(!dir.exists(path.cut2)) dir.create(path.cut2)
fnFs.cut2 <- file.path(path.cut2, basename(fnFs))
fnRs.cut2 <- file.path(path.cut2, basename(fnRs))


FWD.RC2 <- dada2:::rc(FWD2)
REV.RC2 <- dada2:::rc(REV2)


# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD2, "-a", REV.RC2)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV2, "-A", FWD.RC2)


# Run Cutadapt - !Second Run!
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut2[i], "-p", fnRs.cut2[i], # output files
                             fnFs.cut[i], fnRs.cut[i])) # input files
}


# sanity check - !Second Run!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut2[[1]]),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut2[[1]]),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut2[[1]]),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut2[[1]]))
# Success! Primers are no longer detected in cut2files.


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut2, pattern = "_1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut2, pattern = "_2.fq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


filtFs <- file.path(path.cut2, "filtered", basename(cutFs))
filtRs <- file.path(path.cut2, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
    truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# !! minLen is defined by the length distribution of seqs and the identity of seqs of 100bp, 90bp, 70bp.
head(out)

# quality plot file
pdf(file = "QualityProfileForward_nirS_full_1_2.pdf",
     width = 4, # The width of the plot in inches
    height = 4)
plotQualityProfile(filtFs[1:2])
dev.off()

pdf(file = "QualityProfileReverse_nirS_full_1_2.pdf",
     width = 4, # The width of the plot in inches
    height = 4)
plotQualityProfile(filtRs[1:2])
dev.off()


# learn errors
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)


# plot error
pdf(file = "plotErrors_nirS_full_fwd.pdf")   # The directory you want to save the file in
plotErrors(errF, nominalQ = TRUE)
dev.off()



pdf(file = "plotErrors_nirS_full_rev.pdf")   # The directory you want to save the file in
plotErrors(errR, nominalQ = TRUE)
dev.off()


# dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# name derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# [1]   48 2068


# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 964 bimeras out of 2068 input sequences.


# seq length distribution
table(nchar(getSequences(seqtab.nochim)))
# As expected, quite a bit of length variability in the amplified nirS region.


dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,
    getN), rowSums(seqtab.nochim))


# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
    "nonchim")
rownames(track) <- sample.names
head(track)
# kept the majority of our raw reads, run successfully


### ### ### ### ### ###
###  Fasta output   ###
### ### ### ### ### ###

uniquesToFasta(seqtab.nochim, "dada2_nirS_full.fasta")

### ### ### ### ### ###
### Export table   ####
### ### ### ### ### ###

samples.out <- rownames(seqtab.nochim)
asv.1=as.data.frame(t(seqtab.nochim))
asv.1$sum=rowSums(asv.1)
asv.final=asv.1[order(-asv.1$sum),]
asv.final$Sequence=rownames(asv.final)
rownames(asv.final) = sprintf("preASV%04d", 1:nrow(asv.final))
write.table(asv.final, "NirS_full_ASV_table.txt",  quote = FALSE)
write.table(track, "NirS_full_ASV_stats.txt", quote = FALSE)

## DADA2_nirS_full pipeline finished on 13/01/2022 :D ##
