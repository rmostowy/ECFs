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

# Load refseq table filtered by NCBI
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Load reprseq sequences
reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
mmseqs.tmp.dir <- "/Users/rmostowy/Data/tmp-ecfs/recent-rec/mmseqs_tmp"
mmseqs.db <- file.path(mmseqs.tmp.dir, "reprseqDB")
mmseqs.res <- file.path(mmseqs.tmp.dir, "reprseq_vs_reprseq")
mmseqs.res.m8 <- file.path(mmseqs.tmp.dir, "reprseq_vs_reprseq.m8")

# Ensure the tmp directory exists
dir.create(mmseqs.tmp.dir, showWarnings = FALSE, recursive = TRUE)

# Step 1: Create MMseqs2 Database
mmseqs.createdb.cmd <- paste("mmseqs createdb", reprseq.seqs.filename, mmseqs.db)
system(mmseqs.createdb.cmd)

# Step 2: Run MMseqs2 All-vs-All Search
mmseqs.search.cmd <- paste(
  "mmseqs search", mmseqs.db, mmseqs.db, mmseqs.res, mmseqs.tmp.dir,
  "--min-seq-id 0.3 -s 7.5 --threads 12"  # Adjust parameters as needed
)
system(mmseqs.search.cmd)

# Step 3: Convert MMseqs2 Alignments to M8 Format (BLAST-like)
mmseqs.convert.cmd <- paste(
  "mmseqs convertalis", mmseqs.db, mmseqs.db, mmseqs.res, mmseqs.res.m8,
  "--format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov'"
)
system(mmseqs.convert.cmd)

# Step 4: Load BLAST-like results into a data.table
reprseq.vs.reprseq <- fread(mmseqs.res.m8, header = FALSE, sep = "\t")

# Assign meaningful column names
setnames(reprseq.vs.reprseq, c(
  "qname", "sname", "pident", "alnlen", "mismatch", "gapopen",
  "qstart", "qend", "qlength", "sstart", "send", "slength",
  "evalue", "bitscore", "qcov", "scov"
))
fwrite(reprseq.vs.reprseq, file = reprseq.vs.reprseq.blast.output)

# Load reprseq all-by-all cov table
# !!! This needs to be calculated from scratch in the new data release
reprseq.cov.filename <- '/Users/rmostowy/MGG Dropbox/Projects/divRBP/phage-pp-workdir-refseq-hhblits/output/prot-families/families/dataset-full/repr-hits-pairwise-prob95.csv'
reprseq.cov.orig <- fread(reprseq.cov.filename)
reprseq.cov.orig[, max_cov := NULL]
reprseq.cov.orig[, min_cov := NULL]

# adjust pair names so that qname < sname
reprseq.cov.ab <- reprseq.cov.orig
reprseq.cov.ba <- reprseq.cov.orig %>% 
  select(qname = sname, sname = qname, prob,
         scov_min = qcov_min, scov_max = qcov_max,
         qcov_min = scov_min, qcov_max = scov_max, pident)
reprseq.cov.both <- rbind(reprseq.cov.ab, reprseq.cov.ba)
reprseq.cov <- reprseq.cov.both[qname < sname]
setkey(reprseq.cov, qname, sname)

# find hits with a single alignment per pair when filtered at high pident
reprseq.vs.reprseq.prefiltered <- reprseq.vs.reprseq[qname != sname]
setorder(reprseq.vs.reprseq.prefiltered, qname, sname, bitscore)
all.pairs <- reprseq.vs.reprseq.prefiltered %>% select(qname, sname)
some.pairs.are.duplicated <- any(duplicated(all.pairs))
if(some.pairs.are.duplicated) cat("Some pairs are duplicated!\n")

# assign a unique pair.index per pair
all.pairs[, reprseq.a := pmin(qname, sname)]
all.pairs[, reprseq.b := pmax(qname, sname)]
all.pairs[, reprseq.ab := sprintf("%s-%s", reprseq.a, reprseq.b)]
unique.pairs <- unique(all.pairs$reprseq.ab)
no.unique.pairs <- length(unique.pairs)
ab2index <- data.table(ab = unique.pairs, index = 1:no.unique.pairs)
setkey(ab2index, ab)
all.pairs$pair.index <- ab2index[.(all.pairs$reprseq.ab)]$index
reprseq.vs.reprseq.prefiltered$pair.index <- all.pairs$pair.index

# consider only pairs with bidirectional hits
all.ab <- reprseq.vs.reprseq.prefiltered[qname < sname]$pair.index
all.ba <- reprseq.vs.reprseq.prefiltered[qname > sname]$pair.index
bidirectional.pairs <- intersect(all.ab, all.ba)
reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq.prefiltered[pair.index %chin% bidirectional.pairs]
reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq.filtered[qname < sname]

# adjust qcov/scov table to match the high pident table
reprseq.vs.reprseq.filtered.qcov.scov <- reprseq.cov[.(reprseq.vs.reprseq.filtered[, .(qname, sname)])]
reprseq.vs.reprseq.filtered.qcov.scov[, qcov := pmax(qcov_min, qcov_max, na.rm = T)]
reprseq.vs.reprseq.filtered.qcov.scov[, scov := pmax(scov_min, scov_max, na.rm = T)]

reprseq.vs.reprseq.full <- cbind(reprseq.vs.reprseq.filtered, reprseq.vs.reprseq.filtered.qcov.scov %>% select(pident.full = pident, qcov.full = qcov, scov.full = scov))
reprseq.vs.reprseq.full[, pident.full := pident.full*100]
reprseq.vs.reprseq.full[, q.hit.length := qend - qstart + 1]
reprseq.vs.reprseq.full[, s.hit.length := send - sstart + 1]

# Filter to create a table of putative recent recombinations
recent.recombination.table <- reprseq.vs.reprseq.full %>% 
  filter(q.hit.length >= RECENT.REC.MIN.Q.HIT.LENGTH &
           s.hit.length >= RECENT.REC.MIN.S.HIT.LENGTH &
           pident >= RECENT.REC.MIN.PIDENT &
           qcov.full <= RECENT.REC.MAX.QCOV &
           scov.full <= RECENT.REC.MAX.SCOV)
recent.recombination.table.final <- recent.recombination.table %>% 
  select(qname, sname, qstart, qend, qlength, sstart, send, slength,
         pident, qcov, qcov.full, scov, scov.full, 
         hit.length.query = q.hit.length, hit.length.subject = s.hit.length,
         pair.index)
setorder(recent.recombination.table.final, qname, sname)

# write to a file
fwrite(recent.recombination.table.final, file = reprseq.recent.rec.table.output)

