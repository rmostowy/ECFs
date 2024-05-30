### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))

### locate the RStudio project root and load the configuration file
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

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
# Mark sequences that have missing characters in the name table
name.table$has.missing.char <- name.table$ncbi.id %in% reprseqs.with.unknown.chars
# Get the representative names of sequences with missing characters
reprseqs.with.missing.chars <- unique(name.table$repr.name[name.table$has.missing.char])

# Load tables with hits to other representative sequences and PHROGs
reprseq.vs.reprseq <- fread(reprseq.vs.reprseq.filename)
reprseq.vs.phrogs <- fread(phrog.hhpred.hits.filename)
reprseq.vs.phrogs <- reprseq.vs.phrogs[reprseq.vs.phrogs$prob>=min(reprseq.vs.reprseq$prob)]

# Identify reprseqs that have hits to other reprseqs
reprseq.with.hits.to.other.reprseq <- sort(unique(c(reprseq.vs.reprseq$qname, reprseq.vs.reprseq$sname)))
# Identify reprseqs that do not have hits to other reprseqs
reprseq.without.hits.to.other.reprseq <- unique(name.table$repr.name[!name.table$repr.name %in% reprseq.with.hits.to.other.reprseq])
# Identify reprseqs that have hits to PHROGs
reprseq.with.hits.to.phrogs <- sort(unique(reprseq.vs.phrogs$qname))
# Identify reprseqs that do not have hits to PHROGs
reprseq.without.hits.to.phrogs <- unique(name.table$repr.name[!name.table$repr.name %in% reprseq.with.hits.to.phrogs])

# Identify reprseqs that do not have hits to PHROGs or other reprseqs
reprseq.weird.singletons <- intersect(reprseq.without.hits.to.phrogs, reprseq.without.hits.to.other.reprseq)

# Create a vector of all sequences to remove
reprseqs.to.remove <- sort(unique(c(reprseq.weird.singletons, reprseqs.with.missing.chars)))
# Identify sequences that should not be removed
reprseqs.not.to.remove <- unique(name.table$repr.name[!name.table$repr.name %in% reprseqs.to.remove])

# Load raw family data
reprseq.families.raw <- readLines(reprseq.families.raw.filename)
no.fam <- length(reprseq.families.raw)

# Clean family lists by removing sequences to be removed
reprseq.families.list.cleaned <- lapply(1:no.fam, function(k){
  this.fam.reprseqs <- strsplit(reprseq.families.raw[[k]], "\t")[[1]]
  this.fam.reprseqs[!this.fam.reprseqs %in% reprseqs.to.remove]
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


# merge all tables into one, removing the iffy reprseqs
name.table.cleaned <- name.table[!name.table$repr.name %in% reprseqs.to.remove]
name.table.cleaned[, has.missing.char := NULL]

# Merge cleaned name table with family table
name.table.final <- merge(name.table.cleaned, family.table, by.x = "repr.name", by.y = "reprseq", all.x = TRUE)
setnames(name.table.final, old = "repr.name", new = "reprseq")

# Write the final name table to a file
fwrite(name.table.final, file = reprseq.table.filename)

# Filter sequences based on the final name table
name.table.final.repr <- name.table.final[name.table.final$cluster.name == name.table.final$ncbi.id]
filtered_seqs <- all.cds[names(all.cds) %in% name.table.final.repr$cluster.name]
new_names <- name.table.final.repr$reprseq[match(names(filtered_seqs), name.table.final.repr$cluster.name)]
names(filtered_seqs) <- new_names
sorted_seqs <- filtered_seqs[order(names(filtered_seqs))]

# Write the filtered and sorted sequences to a file
writeXStringSet(sorted_seqs, filepath = reprseq.seqs.filename)
