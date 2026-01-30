### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# helper stuff
no.threads <- 20L
max_targets = 100000L
blast_bin <- NULL

exe <- function(cmd) {
  if (is.null(blast_bin)) return(cmd)
  file.path(blast_bin, cmd)
}

# refseq.genomes.filename

# Create temporary directory
blast.tmp.dir <- "/Users/rmostowy/Data/tmp-ecfs/genome-megablast"
dir.create(blast.tmp.dir, showWarnings = FALSE, recursive = TRUE)

# Define helper files
db.prefix <- file.path(blast.tmp.dir, "refseq_genomes_db")
out.tsv <- file.path(blast.tmp.dir, "megablast_allbyall.tsv")

# Copy input FASTA file into the tmp directory to avoid issues
stopifnot(file.exists(refseq.genomes.filename))
src_fasta <- normalizePath(refseq.genomes.filename, mustWork = TRUE)

local_fasta <- file.path(blast.tmp.dir, "refseq-genomes.fa")
copying.fasta <- file.copy(src_fasta, local_fasta, overwrite = TRUE)
stopifnot(file.exists(local_fasta))

# Build BLAST DB
cat("Building BLAST DB...", "\n")
status <- system2(
  command = exe("makeblastdb"),
  args    = c("-in", local_fasta,
              "-dbtype", "nucl",
              "-parse_seqids",
              "-out", db.prefix),
  stdout  = TRUE,
  stderr  = TRUE
)
if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
  stop("makeblastdb failed.\n", paste(status, collapse = "\n"))
}
cat("Done!", "\n\n")

# Run MEGABLAST
outfmt <- paste(
  "6",
  "qseqid sseqid pident length mismatch gapopen",
  "qstart qend sstart send",
  "qlen slen evalue bitscore",
  sep = " "
)

cat("Running MegaBLAST...", "\n")
status2 <- system2(
  command = exe("blastn"),
  args    = c("-task", "megablast",
              "-query", local_fasta,
              "-db", db.prefix,
              "-evalue", format(GENOME.MEGABLAST.EVALUE, scientific = TRUE),
              "-num_threads", as.integer(no.threads),
              "-max_target_seqs", as.integer(max_targets),
              "-outfmt", shQuote(outfmt),
              "-out", out.tsv),
  stdout  = TRUE,
  stderr  = TRUE
)
if (!is.null(attr(status2, "status")) && attr(status2, "status") != 0L) {
  stop("blastn failed.\n", paste(status2, collapse = "\n"))
}
cat("Done!", "\n\n")


# read output table with data.table 
cat("Processing output...", "\n")
hits <- fread(out.tsv, header = FALSE, sep = "\t")
setnames(hits, c(
  "qname","sname","pident","alnlen","mismatch","gapopen",
  "qstart","qend","sstart","send","qlength","slength","evalue","bitscore"
))

# optional: coerce numeric/integer columns
int_cols <- c("alnlen","mismatch","gapopen","qstart","qend","sstart","send","qlength","slength")
hits[, (int_cols) := lapply(.SD, as.integer), .SDcols = int_cols]
hits[, `:=`(pident = as.numeric(pident),
            evalue = as.numeric(evalue),
            bitscore = as.numeric(bitscore))]
hits[, sname := sub("^ref\\|([^|]+)\\|$", "\\1", sname)]

fwrite(hits, file = genomes.aba.table.output)
cat("Done!", "\n\n")

# clean temporary dir
cat("Cleaning tmp directory...", "\n")
delete.dir(blast.tmp.dir)
cat("Done!", "\n\n")

