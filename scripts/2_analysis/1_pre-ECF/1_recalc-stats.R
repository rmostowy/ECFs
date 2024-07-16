### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# COMMENT: This file should become obsolete once Krzysztof recalculates
# qECF start and stop coordinates correctly

######################################
##### 1. Define paths and parameters
######################################
aba <- fread(reprseq.vs.reprseq.clust.filename)
setkey(aba, qname)
stats <- fread(qecf.stats.filename)
stats[, V1 := NULL]
no.cores <- 12

######################################
##### 2. Define functions
######################################
get.overlaps.fast <- function(coordinate.df){
  # convert
  colnames(coordinate.df) <- c("start", "end")
  fragments <- as.data.table(coordinate.df)
  
  # Add a unique identifier for each row
  fragments[, id := .I]
  
  # Self-join to compare all pairs, excluding self-comparisons
  pairs <- fragments[fragments, on=.(id < id), allow.cartesian = TRUE]
  
  # Calculate intersection and union lengths for each pair
  pairs[, `:=` (intersection_length = pmax(0, pmin(i.end, end) - pmax(i.start, start) + 1),
                union_length = pmax(end, i.end) - pmin(start, i.start) + 1)]
  
  # Calculate overlap proportion for each pair
  pairs[, overlap_proportion := intersection_length / union_length]
  
  # Calculate the mean overlap, including zeros
  mean_overlap <- mean(pairs$overlap_proportion, na.rm = TRUE)
  
  # Return the mean overlap
  return(mean_overlap)
}

coverage_proportion <- function(start, end, span_start, span_end) {
  # Calculate the overlap between the fragment and the span
  overlap_start <- max(start, span_start)
  overlap_end <- min(end, span_end)
  overlap_length <- max(0, overlap_end - overlap_start + 1)
  
  # Calculate the total length of the span
  span_length <- span_end - span_start + 1
  
  # Return the proportion of the span covered by the fragment
  overlap_length / span_length
}


######################################
##### 3. Calculate stats
######################################
sh <- aba[aba$qname == aba$sname]
sh[, qcov := (qend - qstart + 1) / qlength]
sh[, scov := (send - sstart + 1) / slength]
# setkey(sh, qname)

# remove self-hits
aba <- aba[aba$qname != aba$sname]
all.reprseqs <- unique(aba$qname)

# recalculate stats
my.stats.list <- mclapply(all.reprseqs, function(query.name){
  query.aba <- aba[.(query.name)]
  no.clusters <- max(query.aba$cluster)+1
  query.stats.list <- lapply(0:max(query.aba$cluster), function(cluster.index){
    this.cluster.frags <- query.aba[query.aba$cluster==cluster.index]
    qlength <- this.cluster.frags$qlength[1]
    fragments <- this.cluster.frags[,c("qstart", "qend")]
    no.fragments <- nrow(fragments)
    # fragment.overlaps <- get.overlaps(fragments)
    average.frag.overlap <- get.overlaps.fast(as.data.frame(fragments))
    fragments.start <- fragments$qstart
    fragments.end <- fragments$qend
    qecf.start <- floor(as.numeric(quantile(fragments.start)["25%"]))
    qecf.end <- ceiling(as.numeric(quantile(fragments.end)["75%"]))
    qecf.span <- qecf.end - qecf.start + 1
    frag.coverage <- sapply(1:no.fragments, function(i) coverage_proportion(fragments.start[i], fragments.end[i], qecf.start, qecf.end))
    average.coverage <- mean(frag.coverage)
    max.span <- max(fragments$qend) - min(fragments$qstart) + 1
    frag.lengths <- fragments$qend-fragments$qstart+1
    al2ms <- mean(frag.lengths)/max.span
    data.table(qname = query.name, qlength,
               cluster=cluster.index, nb.hits = no.fragments, average.frag.overlap, qecf.span, average.coverage, span2ql = qecf.span/qlength, qecf.start, qecf.end, al2ms)  
  })
  query.stats <- rbindlist(query.stats.list)
}, mc.cores = no.cores, mc.silent = F)
my.stats <- rbindlist(my.stats.list)

######################################
##### 4. Process and save data
######################################

# Round all float columns to 4 significant digits
my.stats[, (names(my.stats)) := lapply(.SD, function(x) {
  if (is.numeric(x)) return(round(x, 4))
  else return(x)
})]

fwrite(my.stats, file = "reprseq-all-stat-mine.txt")

