### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# merge original PHROG tables into a single table
phrogs.functions <- fread(phrogs.annot.filename)
phrogs.table.full <- fread(phrogs.table.filename, header = F)

# Read the custom PHROG mapping table from an Excel file
phrogs.custom <- as.data.table(read_excel(manual.phrog.mapping.filename, 1))
# Merge the annotation and full PHROG tables on the PHROG ID (phrog and V1 columns)
phrogs.merged <- merge.data.table(phrogs.functions, phrogs.table.full, by.x = "phrog", by.y = "V1")

# Create a new data table with selected columns and additional annotation
phrogs.functions <- data.table(phrog = phrogs.merged$phrog,
                               funct = phrogs.merged$annot,
                               class = phrogs.merged$category,
                               annotation.orig = sprintf("%s|%s", phrogs.merged$annot, phrogs.merged$category),
                               no.seqs = phrogs.merged$V5)

# Filter out singletons, NA functions, and "unknown function"
phrogs.functions <- phrogs.functions[phrogs.functions$no.seqs>1]
phrogs.functions <- phrogs.functions[!is.na(phrogs.functions$funct)]
phrogs.functions <- phrogs.functions[phrogs.functions$funct != "unknown function"]

# Create original and new annotation columns in the custom PHROG mapping table
phrogs.custom$annot.orig <- sprintf("%s|%s", phrogs.custom$funct.orig, phrogs.custom$class.orig)
phrogs.custom$annot.new <- sprintf("%s|%s", phrogs.custom$funct.new, phrogs.custom$class.new)
# Set the key for fast lookup
setkey(phrogs.custom, annot.orig)

# Map old annotations to new annotations in the PHROG functions table
phrogs.functions$annotation.new <- phrogs.custom[.(phrogs.functions$annotation.orig)]$annot.new
phrogs.functions$funct.new <- phrogs.custom[.(phrogs.functions$annotation.orig)]$funct.new
phrogs.functions$class.new <- phrogs.custom[.(phrogs.functions$annotation.orig)]$class.new

# Create the final PHROG functions table with updated annotations
phrogs.functions.final <- data.table(phrog = phrogs.functions$phrog,
                                     funct = phrogs.functions$funct.new,
                                     class = phrogs.functions$class.new,
                                     annotation = phrogs.functions$annotation.new,
                                     no.seqs = phrogs.functions$no.seqs,
                                     annotation.orig = phrogs.functions$annotation.orig)

# Calculate the number of sequences for each function, class, and annotation
funct.count.seqs <- phrogs.functions.final[, sum(no.seqs), by = .(funct, class, annotation)]
# Calculate the number of PHROGs for each function, class, and annotation
funct.count.phrogs <- phrogs.functions.final[, .N, by = .(funct, class, annotation)]
# Merge the sequence and PHROG counts
funct.count <- merge.data.table(funct.count.seqs, funct.count.phrogs, by = c("funct", "class", "annotation"))

# Rename columns for clarity
setnames(funct.count, old = c("V1", "N"), new = c("no.seqs", "no.phrogs"))
# Order the function count table by class, number of sequences, and function
funct.count <- funct.count[order(funct.count$class, -funct.count$no.seqs, funct.count$funct)]
# Add an index column
funct.count <- data.table(annot.index=1:nrow(funct.count), funct.count)

# Write the final PHROG functions table and function count table to output files
fwrite(phrogs.functions.final, file = phrogs.functions.output)
fwrite(funct.count, file = phrogs.count.output)

# Get antidefence sequences
antidefence.seqs <- readAAStringSet(antidefence.seqs.filename)
antidefence.seqs.n <- names(antidefence.seqs)

# Get antidefence table
antidefence.table <- fread(antidefence.descr.filename)
end.row <- which(antidefence.table$Name == "*References:")-2
antidefence.table <- antidefence.table[1:end.row]

# Format the antidefence table
antidefence.table[, base.accession := sub("\\.[0-9]+$", "", Accession)]
antidefence.table$has.seq <- antidefence.table$Accession %in% antidefence.seqs.n
antidefence.table.processed <- antidefence.table %>% 
  select(accession = base.accession, accession.full = Accession, name = Name, funct = Function, annotation, has.seq)
antidefence.table.processed.with.seq <- antidefence.table.processed[has.seq == TRUE]
antidefence.table.processed.with.seq[, has.seq := NULL]

fwrite(antidefence.table.processed.with.seq, file = antidefence.table.output)

