### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load refseq table
reprseq.table <- fread(reprseq.table.filename)
reprseq.table[, cds := gsub("\\.(\\d+):", ":", ncbi.id)]
reprseq.table.cds <- reprseq.table %>% select(cds)
# reprseq.table[, cds := NULL]
all.reprseq <- unique(reprseq.table$reprseq)

# Load NCBI (metadata + CDS) table
metadata.table <- fread(metadata.table.output)
ncbi.cds <- fread(accessions.cds.output)
ncbi.cds[, cds := sprintf("%s:%d-%d(%s)", accession, start, end, strand)]
ncbi.cds.short <- ncbi.cds %>% select(accession, start, end, strand)
setkey(ncbi.cds.short, accession, start, end)

# Modify reprseq.cds
reprseq.table.cds[, c("accession", "range", "strand") := tstrsplit(cds, "[:()]", keep = c(1, 2, 3))]
reprseq.table.cds[, c("start", "end") := tstrsplit(range, "-", keep = c(1, 2))]
reprseq.table.cds[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
reprseq.table.cds[, width := end - start + 1]
reprseq.table.cds <- reprseq.table.cds[, .(accession, start, end, width, strand, cds)]
reprseq.table.cds.short <- reprseq.table.cds %>% select(accession, start, end, strand, cds)
setkey(reprseq.table.cds.short, accession, start, end)
# setorder(reprseq.table.cds, accession, start)

# Merge two tables
merged.own.ncbi <- merge(reprseq.table.cds.short, ncbi.cds.short, by = "accession", suffixes = c("", ".ncbi"), allow.cartesian = T)

# Function to calculate overlap between two intervals
calculate_overlap <- function(start1, end1, start2, end2) {
  # Calculate the intersection (overlap length)
  overlap <- pmax(0, pmin(end1, end2) - pmax(start1, start2) + 1)
  
  # Calculate the union (total combined region covered by both CDSs)
  total_region <- (end1 - start1 + 1) + (end2 - start2 + 1) - overlap
  
  # Compute Jaccard similarity (Intersection over Union)
  jaccard <- overlap / total_region
  
  # Handle cases where there's no valid overlap
  jaccard[is.na(jaccard)] <- 0  # If NA, set to 0
  
  return(jaccard)
}

merged.own.ncbi[, overlap := calculate_overlap(start, end, start.ncbi, end.ncbi)]

# Compute the max overlap per unique CDS
cds_max_overlap_table <- merged.own.ncbi[, .(
  accession = first(accession),  # Take the first occurrence per CDS
  start = first(start),
  end = first(end),
  strand = first(strand),
  max_overlap = max(overlap, na.rm = TRUE)  # Get the max overlap per CDS
), by = cds]

merged.own.ncbi[, max_overlap := max(overlap), by = cds]
setorder(merged.own.ncbi, accession, start, end, strand, cds, -max_overlap)
reprseq.table.ncbi.overlap <- merged.own.ncbi[, .SD[1], by = .(accession, start, end, strand, cds)] %>% select(accession, start, end, strand, cds, max.overlap = max_overlap)

# Consider only CDSs which overlap with NCBI predictions
reprseq.table.ncbi.overlap.pass <- reprseq.table.ncbi.overlap[max.overlap >= MIN.CDS.OVERLAP.WITH.NCBI]

# Remove reprseq which do not contain filtered CDSs
reprseq.table.cds.filtered <- reprseq.table[cds %chin% reprseq.table.ncbi.overlap.pass$cds]
fwrite(reprseq.table.cds.filtered, file = reprseq.table.ncbi.filtered.filename)
