### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Define output folders
table.dir <- "tables"
output.dir <- file.path(visual.dir, "analysis_02-21", "overlaps")
make.dir(output.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Load recent recombination table
recent.rec.table <- fread(reprseq.recent.rec.table.output)
recent.rec.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]
recent.rec.table[, row.id := .I]

# Define functions
rel.pos <- function(pos, length) (pos-1)/(length-1)
# Function to assign multiple bins for each protein region
assign_bins <- function(start, end, bins = 99) {
  breakpoints <- seq(0, 1, length.out = bins + 1)  # Define bin boundaries
  
  # Find all bins that the start and end fall into
  start_bin <- findInterval(start, breakpoints, rightmost.closed = TRUE)
  end_bin <- findInterval(end, breakpoints, rightmost.closed = TRUE)
  
  # Assign to all bins in the range
  return(sprintf("Bin_%.*d", 2, start_bin:end_bin))
  # return(paste0("Bin_", start_bin:end_bin))
}

###############################################
### A) Plot density map of recombination overlaps in 2D
###############################################

# Calculate relative positions per rec pairs
rr <- recent.rec.table
pairs.rel.pos.list <- lapply(1:nrow(rr), function(pair.index){
  # this.rr <- recent.rec.table[pair.index]
  this.rr <- rr[pair.index]
  protein1.name <- this.rr$qname
  protein2.name <- this.rr$sname
  
  protein1.length <- this.rr$qlength
  aln1.start <- this.rr$qstart
  aln1.end <- this.rr$qend
  aln1.length <- aln1.end - aln1.start + 1
  
  prot1.rel.start <- round(rel.pos(aln1.start, protein1.length), d = 2)
  prot1.rel.end <- round(rel.pos(aln1.end, protein1.length), d = 2)
  
  protein2.length <- this.rr$slength
  aln2.start <- this.rr$sstart
  aln2.end <- this.rr$send
  aln2.length <- aln2.end - aln2.start + 1
  
  prot2.rel.start <- round(rel.pos(aln2.start, protein2.length), d = 2)
  prot2.rel.end <- round(rel.pos(aln2.end, protein2.length), d = 2)
  
  data.table(protein1.name, prot1.rel.start, prot1.rel.end,
             protein2.name, prot2.rel.start, prot2.rel.end)
})
pairs.rel.pos <- rbindlist(pairs.rel.pos.list)

# Apply function row-wise and store as lists
pairs.rel.pos[, prot1_bins := lapply(seq_len(.N), function(i) 
  assign_bins(prot1.rel.start[i], prot1.rel.end[i]))]
pairs.rel.pos[, prot2_bins := lapply(seq_len(.N), function(i) 
  assign_bins(prot2.rel.start[i], prot2.rel.end[i]))]

pairs.rel.pos[, row_id := .I]

# Expand protein1 and protein2 bins using a cross join (Cartesian join)
pairs.rel.pos.expanded <- pairs.rel.pos[
  , {
    # Perform a cross join between all bins for `prot1_bins` and `prot2_bins`
    expanded_bins <- CJ(prot1_bin = prot1_bins[[1]], prot2_bin = prot2_bins[[1]])
    # Add row_id only once and return expanded table row by row
    expanded_bins
  }, by = .(protein1.name, protein2.name, row_id)  # Ensure row_id is preserved
]

# Count occurrences of each (prot1_bin, prot2_bin) combination
heatmap_data <- pairs.rel.pos.expanded[, .N, by = .(prot1_bin, prot2_bin)]
# Extract numeric values from bin names (assuming format "Bin_X")
heatmap_data[, prot1_rel := as.numeric(gsub("Bin_", "", prot1_bin)) / max(as.numeric(gsub("Bin_", "", prot1_bin)))]
heatmap_data[, prot2_rel := as.numeric(gsub("Bin_", "", prot2_bin)) / max(as.numeric(gsub("Bin_", "", prot2_bin)))]

# Plot with relative positions on both axes
p.overlap.heatmap <- ggplot(heatmap_data, aes(x = prot1_rel, y = prot2_rel, fill = N)) +
  geom_tile() +
  scale_fill_viridis_c() +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0,1)) +  # X-axis from 0 to 1
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0,1)) +  # Y-axis from 0 to 1
  theme_minimal() +
  labs(title = "Heatmap of Overlapping Regions",
       x = "Protein 1 Relative Position",
       y = "Protein 2 Relative Position",
       fill = "Count")

p.overlap.heatmap.outfile <- file.path(output.dir, "plot1_overlap-heatmap.png")
ggsave(p.overlap.heatmap.outfile, p.overlap.heatmap, width = 10, height = 10)

###############################################
### B) Create a table of N-o-C assignment
###############################################

MY.MAX.DIST <- 6

assign.class <- function(start, end, length, max.dist.from.edge = 3){
  # is.nterm <- is.cterm <- is.other <- FALSE
  category <- NA
  rel.start <- rel.pos(start, length)
  rel.end <- rel.pos(end, length)
  
  is.upstream <- rel.end<0.5
  is.downstream <- rel.start>0.5
  
  if(is.upstream){
    dist.from.edge <- abs(start)
    if(dist.from.edge <= max.dist.from.edge){
      category <- "N"
    } else{
      category <- "o"
    }
  } else if(is.downstream){
    dist.from.edge <- abs(length - end + 1)
    if(dist.from.edge <= max.dist.from.edge){
      category <- "C"
    } else{ category <- "o" }
  } else{
    category <- "o"
  }
  return(category)
}


rr <- recent.rec.table

# MAX.DIST=
relative.classes.list <- lapply(1:nrow(rr), function(pair.index){
  
  this.rr <- rr[pair.index]
  protein1.name <- this.rr$qname
  protein2.name <- this.rr$sname
  
  protein1.length <- this.rr$qlength
  aln1.start <- this.rr$qstart
  aln1.end <- this.rr$qend
  protein1.class <- assign.class(aln1.start, aln1.end, protein1.length, max.dist.from.edge = MY.MAX.DIST)
  
  protein2.length <- this.rr$slength
  aln2.start <- this.rr$sstart
  aln2.end <- this.rr$send
  protein2.class <- assign.class(aln2.start, aln2.end, protein2.length, max.dist.from.edge = MY.MAX.DIST)
  
  data.table(protein1.name, aln1.start, aln1.end, protein1.length, protein1.class,
             protein2.name, aln2.start, aln2.end, protein2.length, protein2.class)
})
relative.classes <- rbindlist(relative.classes.list)
relative.classes[, class12 := sprintf("%s-%s", protein1.class, protein2.class)]
relative.classes[, rowid := .I]

# Convert to long format for ggupset
long_format <- relative.classes[, .(rowid, class12)]
long_format[, class_list := strsplit(class12, "-")]
long_format[, class_list := lapply(class_list, function(x) factor(x, levels = c("N", "o", "C")))]

# Create UpSet plot using ggupset
p.relative.classes <- ggplot(long_format, aes(x = class_list)) +
  geom_bar() +  # Creates the intersection counts
  scale_x_upset(sets = c("N", "o", "C")) +  # Use UpSet-style X axis
  labs(title = "UpSet Plot of Recent recombinations",
       x = "Intersection of Classes",
       y = "Count of Overlaps") +
  theme_minimal()
p.relative.classes.outfile <- file.path(output.dir, "plot2_rel-classes-upset.png")
ggsave(p.relative.classes.outfile, p.relative.classes, width = 7, height = 5)

p.relative.classes.p95 <- ggplot(long_format[recent.rec.table[pident >= 95]$row.id], aes(x = class_list)) +
  geom_bar() +  # Creates the intersection counts
  scale_x_upset(sets = c("N", "o", "C")) +  # Use UpSet-style X axis
  labs(title = "UpSet Plot of Recent recombinations",
       x = "Intersection of Classes",
       y = "Count of Overlaps") +
  theme_minimal()
p.relative.classes.p95.outfile <- file.path(output.dir, "plot3_rel-classes-upset-p95.png")
ggsave(p.relative.classes.p95.outfile, p.relative.classes.p95, width = 7, height = 5)


relative.classes.outfile <- file.path(table.dir, "relative-classes.txt")
fwrite(relative.classes, file = relative.classes.outfile)
