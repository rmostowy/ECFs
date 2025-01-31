### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))

### locate the RStudio project root and load the configuration file
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

##########################################
# 1. Identify reprseqs to remove 
#    based on missing characters
##########################################

# Read all coding sequences from the specified path using Biostrings
all.cds <- readAAStringSet(all.cds.filename, format = "fasta")
# Convert the AAStringSet object to a character vector
all.cds.s <- as.character(all.cds)

# count the number of times unknown characters (X) occur
count.number.of.x <- function(sequence) length(gregexpr("X", sequence)[[1]])
num.of.x <- sapply(all.cds.s, count.number.of.x)
# Identify sequences with more than the allowed number of 'X' characters
reprseqs.with.unknown.chars <- names(which(num.of.x > ACCEPTED.NUM.UNKNOWN.AA))

# Load the name table that maps NCBI IDs to representative names
name.table <- fread(name.table.filename)
name.table[, cds := gsub("\\.(\\d+):", ":", ncbi.id)]
# Mark sequences that have missing characters in the name table
name.table$has.missing.char <- name.table$ncbi.id %in% reprseqs.with.unknown.chars

##########################################
# 2. Identify reprseqs to remove 
#    based on overlap with NCBI
##########################################

# Format CDS predictions from the reprseq dataset
reprseq.table <- name.table %>% select(ncbi.id)
reprseq.table[, cds := gsub("\\.(\\d+):", ":", ncbi.id)]
reprseq.table[, c("accession", "range", "strand") := tstrsplit(cds, "[:()]", keep = c(1, 2, 3))]
reprseq.table[, c("start", "end") := tstrsplit(range, "-", keep = c(1, 2))]
reprseq.table[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
reprseq.table[, width := end - start + 1]
reprseq.table.cds <- reprseq.table[, .(accession, start, end, width, strand, cds)]
reprseq.table.cds.short <- reprseq.table.cds %>% select(accession, start, end, strand, cds)
setkey(reprseq.table.cds.short, accession, start, end)

# Load NCBI CDS predictions
ncbi.cds <- fread(accessions.cds.output)
ncbi.cds[, cds := sprintf("%s:%d-%d(%s)", accession, start, end, strand)]
ncbi.cds.short <- ncbi.cds %>% select(accession, start, end, strand)
setkey(ncbi.cds.short, accession, start, end)

# Merge two tables
merged.own.ncbi <- merge(reprseq.table.cds.short, ncbi.cds.short, by = "accession", suffixes = c("", ".ncbi"), allow.cartesian = T)

# Calculate overlap
merged.own.ncbi[, overlap := calculate_overlap_jaccard(start, end, start.ncbi, end.ncbi)]

# Calculate maximum overlap per reprseq CDS
merged.own.ncbi[, max.overlap := max(overlap), by = cds]
setorder(merged.own.ncbi, accession, start, end, strand, cds, -max.overlap)
reprseq.table.cds.short.overlap <- merged.own.ncbi[, .SD[1], by = .(accession, start, end, strand, cds)] %>% select(accession, start, end, strand, cds, overlap = max.overlap)
reprseq.table.cds.short.overlap[, .(cds, overlap)]

##########################################
# 3. Merge results and filter
##########################################

name.table.full <- merge(name.table, reprseq.table.cds.short.overlap[, .(cds, overlap)], by = "cds")
name.table.full[, sufficient.overlap := overlap >= MIN.CDS.OVERLAP.WITH.NCBI]

# Filter protein sequences
name.table.cleaned <- name.table.full[has.missing.char == FALSE & sufficient.overlap == TRUE]
all.reprseqs <- sort(unique(name.table.cleaned$repr.name))

##########################################
# 4. Extract and add family information
##########################################

# Load raw family data
reprseq.families.raw <- readLines(reprseq.families.raw.filename)
no.fam <- length(reprseq.families.raw)

# Clean family lists by removing sequences to be removed
reprseq.families.list.cleaned <- lapply(1:no.fam, function(k){
  this.fam.reprseqs <- strsplit(reprseq.families.raw[[k]], "\t")[[1]]
  this.fam.reprseqs[this.fam.reprseqs %in% all.reprseqs]
})
# Get the lengths of cleaned family lists
reprseq.families.list.cleaned.l <- sapply(reprseq.families.list.cleaned, length)
# Filter out empty families
reprseq.families.list.cleaned.corrected <- reprseq.families.list.cleaned[reprseq.families.list.cleaned.l>0]
no.fam.corrected <- length(reprseq.families.list.cleaned.corrected)

# Create a table with family assignments for representative sequences
family.table.list <- lapply(1:no.fam.corrected, function(index){
  reprseq <- reprseq.families.list.cleaned.corrected[[index]]
  family <- sprintf("fam%.*d", nchar(no.fam.corrected), index)
  data.table(family = rep(family, length(reprseq)), reprseq)
})
family.table <- rbindlist(family.table.list)

##########################################
# 5. Finalise the dataset
##########################################

# Merge cleaned name table with family table
name.table.cleaned.fam <- merge(name.table.cleaned, family.table, by.x = "repr.name", by.y = "reprseq", all.x = TRUE)
name.table.final <- name.table.cleaned.fam %>% select(reprseq = repr.name, cluster.name, ncbi.id, family, ncbi.overlap = overlap)
name.table.final[, ncbi.overlap := round(ncbi.overlap, d = 4)]

# Write the final name table to a file
fwrite(name.table.final, file = reprseq.table.filename)

# Write removed reprseqs to a file
reprseqs.to.remove <- name.table$repr.name[!name.table$repr.name %in% all.reprseqs]
write(reprseqs.to.remove, ncolumns = 1, reprseq.removed.filename)

# Filter sequences based on the final name table
# IMPORTANT: this step assumes that if the representative sequence was removed
# (for example due to insufficient overlap with NCBI) but other sequences in
# the cluster were not, then the representative sequence stays
# Note that it's a temporary solution before the final dataset is ready
filtered_seqs <- all.cds[names(all.cds) %in% name.table.final$cluster.name]
new_names <- name.table.final$reprseq[match(names(filtered_seqs), name.table.final$cluster.name)]
names(filtered_seqs) <- new_names
sorted_seqs <- filtered_seqs[order(names(filtered_seqs))]

# Write the filtered and sorted sequences to a file
writeXStringSet(sorted_seqs, filepath = reprseq.seqs.filename)

