# AMB DADA2 script
#
# this script expects as input a directory with illumina paired-end reads,
# with adapters removed. each sample should be in its own directory <SAMPLE NAME>,
# and files containing reads should be named <SAMPLE NAME>_<1/2>.fastq(.gz)
#
# the script will produce filtered reads in the sample directories named
# <SAMPLE NAME>.filtered_<1/2>.fastq.gz
#
# this script is based on the DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html

rm(list=ls())

####################################################
### load required packages (install if necessary) ##
####################################################

install_if_needed <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# ensure BiocManager is installed first
install_if_needed("BiocManager")

# bioconductor packages
install_if_needed("dada2", bioc=T)
install_if_needed("phyloseq", bioc=T)
install_if_needed("DECIPHER", bioc=T)
install_if_needed("phangorn", bioc=T)
install_if_needed("ggplot2", bioc=F)
install_if_needed("tictoc", bioc=F)

library(parallel)
library(dada2)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(tictoc)

##########################
### TUNABLE PARAMETERS ###
##########################

# 16S files location
# if not set, ask the user for input directory
#input_dir <- "../16S_sequences/"

# output directory
# if not set, use DADA2_output/ relative to the input directory
#output_dir <- "~/Documents/dada2"

# set number of computational cores / CPUs used
n_cores <- detectCores()

# filterAndTrim parameters
truncLen<-0           # Default 0 (no truncation). Truncate reads after truncLen
                      # bases. Reads shorter than this are discarded.
maxEE<-c(2)           # Default Inf (no EE filtering). After truncation, reads with 
                      # higher than maxEE "expected errors" will be discarded.
                      # Expected errors are calculated from the nominal definition 
                      # of the quality score: EE = sum(10^(-Q/10))
truncQ<-2             # Default 2. Truncate reads at the first instance of a quality 
                      # score less than or equal to truncQ.

# taxonomic inference parameters
# use the SILVA db from here: https://www2.decipher.codes/Downloads.html
tax_db <- '~/Projects/resources/SILVA_SSU_r138_2_2024.RData'

###########################
### ACTUAL SCRIPT BELOW ###
###########################

if (!exists('input_dir')) {
  if (.Platform$OS.type == 'unix') {
    input_dir <- readline(prompt = "Enter directory to process: ")
  } else {
    input_dir <- choose.dir(getwd(), "Choose folder to process")
  }
}
if (!exists('output_dir')) output_dir <- file.path(input_dir, "DADA2_output/")  

# use multithreading only if we aren't on windows
multithread <- if (.Platform$OS.type == "windows") FALSE else n_cores

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

# write some basic info to log file
logfile <- file.path(output_dir, 'dada2_output.txt')

pkg_versions <- c(
  dada2 = as.character(packageVersion("dada2")),
  phyloseq = as.character(packageVersion("phyloseq")),
  DECIPHER = as.character(packageVersion("DECIPHER")),
  phangorn = as.character(packageVersion("phangorn"))
)

log_text <- paste0(
  "Script started at: ", format(Sys.time()), "\n\n",
  "Loaded packages:\n",
  paste0("  ", names(pkg_versions), " ", pkg_versions, collapse = "\n"), "\n\n",
  "Running DADA2 with parameters:\n",
  "  truncLen = (", paste(truncLen, collapse = ","), ")\n",
  "  maxEE = (", paste(maxEE, collapse = ","), ")\n",
  "  truncQ = ", truncQ, "\n\n"
)
cat(log_text)
writeLines(log_text, logfile)

cat('Loading files...\n')
tic("DADA2 pipeline")

# find all files matching *_[12].fastq{,.gz} and assign them as forward and reverse
fwd_files <- list.files(path=input_dir, pattern="_1\\.fastq(\\.gz)?$", full.names=T, recursive=T)
rev_files <- list.files(path=input_dir, pattern="_2\\.fastq(\\.gz)?$", full.names=T, recursive=T)

# only keep the trimmed ones
fwd_files <- fwd_files[!grepl("\\.raw_1\\.fastq(\\.gz)?$", fwd_files)]
rev_files <- rev_files[!grepl("\\.raw_2\\.fastq(\\.gz)?$", rev_files)]

# also remove any previously filtered ones
fwd_files <- fwd_files[!grepl("\\.filtered_1\\.fastq(\\.gz)?$", fwd_files)]
rev_files <- rev_files[!grepl("\\.filtered_2\\.fastq(\\.gz)?$", rev_files)]

# simple check to see we have same amount of fwd/rev files
stopifnot(length(fwd_files) == length(rev_files))

# plot aggregate quality data (add more detailed x axis)
cat('Plotting aggregate quality stats...\n')
pf <- plotQualityProfile(fwd_files, aggregate=T) + scale_x_continuous(breaks=seq(0,250,10))
pr <- plotQualityProfile(rev_files, aggregate=T) + scale_x_continuous(breaks=seq(0,250,10))
ggsave(plot=pf, filename = file.path(output_dir, "aggregate_quality_fwd.pdf"))
ggsave(plot=pr, filename = file.path(output_dir, "aggregate_quality_rev.pdf"))

# we'll create the filtered files in the same dir as the raw ones
fwd_filt <- sub("(_1\\.fastq(\\.gz)?)$", ".filtered\\1", fwd_files)
rev_filt <- sub("(_2\\.fastq(\\.gz)?)$", ".filtered\\1", rev_files)

# assign sample names
sample.names <- sapply(strsplit(basename(fwd_files), "_"), `[`, 1)
names(fwd_filt) <- sample.names
names(rev_filt) <- sample.names

# do the trimming and filtering
out <- filterAndTrim(fwd_files, fwd_filt, rev_files, rev_filt, 
                     truncLen=truncLen, maxEE=maxEE, truncQ=truncQ, maxN=0, 
                     rm.phix=TRUE, compress=TRUE, multithread=multithread)

# learn error rates
errF <- learnErrors(fwd_filt, multithread=multithread)
errR <- learnErrors(rev_filt, multithread=multithread)

# inspect results
svg(filename=file.path(output_dir, "errors_fwd.svg") ); plotErrors(errF, nominalQ=TRUE); dev.off() 
svg(filename=file.path(output_dir, "errors_rev.svg") ); plotErrors(errR, nominalQ=TRUE); dev.off()

# let's find out who's there!
dadaFs <- dada(fwd_filt, err=errF, multithread=multithread)
dadaRs <- dada(rev_filt, err=errR, multithread=multithread)

# merge the paired ends
mergers <- mergePairs(dadaFs, fwd_filt, dadaRs, rev_filt, verbose=TRUE)
if (length(sample.names) > 1) names(mergers) <- sample.names

# make the asv table
seqtab <- makeSequenceTable(mergers)
if (length(sample.names) == 1) rownames(seqtab) <- sample.names

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=TRUE)

# check reads throughout the process
getN <- function(x) sum(getUniques(x))
if (length(sample.names) > 1) {
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
} else {
  track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
}
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

# add fraction of non-chimeras column
track <- cbind(track, nonchim.ratio = track[, "nonchim"] / track[, "merged"])
track[, "nonchim.ratio"] <- round(track[, "nonchim.ratio"], 3)
rownames(track) <- sample.names

# output the filter stats to the input dir
write.table(track, file.path(output_dir, "filter_stats.tsv"), row.names=T, col.names=T, sep='\t')

# use DECIPHER IDTAXA method for better inference
# Murali, A., Bhargava, A. & Wright, E.S. IDTAXA: a novel approach for accurate
# taxonomic classification of microbiome sequences. Microbiome 6, 140 (2018).
# https://doi.org/10.1186/s40168-018-0521-5
#
# define the tax_db to use in the tunables section
load(tax_db)

# create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))
ids <- IdTaxa(dna, trainingSet, strand="top", processors=n_cores, verbose=FALSE) # processors=NULL means use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus")

# convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)

write.table(taxid, file.path(output_dir, "taxid.tsv"), row.names=T, col.names=T, sep='\t')
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab.nochim.rds"))
saveRDS(seqtab, file.path(output_dir, "seqtab.rds"))
saveRDS(taxid, file.path(output_dir, "taxa.rds"))

# note the running time
toc()

cat("DADA2 pipeline finished at: ", format(Sys.time()), "\n")
cat("Files are written to: ", "results", "\n")
