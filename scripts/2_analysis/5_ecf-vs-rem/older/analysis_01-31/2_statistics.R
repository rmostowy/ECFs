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

# Load refseq table
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Calculate NCBI overlaps per reprseq
max.overlap.per.reprseq <- reprseq.table[, .(reprseq, ncbi.overlap)]
setorder(max.overlap.per.reprseq, reprseq, -ncbi.overlap)
max.overlap.per.reprseq <- max.overlap.per.reprseq[, .SD[1], by = reprseq]
setkey(max.overlap.per.reprseq, reprseq)

# Load recent recombination table
recent.rec.table <- fread(reprseq.recent.rec.table.output)
recent.rec.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load ECF families
ecf.family.table <- fread(ecf.family.table.output)
setkey(ecf.family.table, ecf)

# Load ECF table
ecf.class.table <- fread(ecf.classification.table.output)
setkey(ecf.class.table, ecf)

ecf.class.table$family <- ecf.family.table[.(ecf.class.table$ecf)]$family
ecf2ecod <- ecf.class.table %>% select(ecf, family, t.id)
ecf2ecod$t.name <- sapply(1:nrow(ecf2ecod), function(k){
  out <- "unknown"
  this.t.id <- ecf2ecod$t.id[k]
  if(this.t.id != "unknown") out <- ecod.domains[t_id == this.t.id]$t_name[1]  
  return(out)
})
setorder(ecf2ecod, family)

# Load functional annotations
reprseq.functions <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions, reprseq)

# Load reprseq vs ECFs
reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecfs$EF <- ecf.family.table[.(reprseq.ecfs$sname)]$family
setkey(reprseq.ecfs, qname)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)


overlap.scores.list <- lapply(1:nrow(recent.rec.table), function(pair.index){
  this.rr <- recent.rec.table[pair.index]
  protein1.name <- this.rr$qname
  protein2.name <- this.rr$sname
  
  protein1.length <- this.rr$qlength
  aln1.start <- this.rr$qstart
  aln1.end <- this.rr$qend
  aln1.length <- aln1.end - aln1.start + 1
  
  protein1.ecod.hits <- reprseq.ecod[.(protein1.name)]
  protein1.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
  setorder(protein1.ecod.hits, -rr.overlap.jaccard)
  protein1.ecod.best.overlap <- protein1.ecod.hits[1]
  ecod.overlap.parms <- calculate_overlap_parms(aln1.start, aln1.end, protein1.ecod.best.overlap$qstart, protein1.ecod.best.overlap$qend)
  ecod.aln.length <- protein1.ecod.best.overlap$qend - protein1.ecod.best.overlap$qstart + 1
  prop.of.ecod.hit.that.overlaps.with.rec <- ecod.overlap.parms$overlap.length/ecod.aln.length
  prop.of.full.ecod.covered.by.overlap <- prop.of.ecod.hit.that.overlaps.with.rec*protein1.ecod.best.overlap$scov
  protein1.ecod.score <- prop.of.full.ecod.covered.by.overlap
  
  protein1.ecf.hits <- reprseq.ecfs[.(protein1.name)]
  protein1.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
  setorder(protein1.ecf.hits, -rr.overlap.jaccard)
  protein1.ecf.best.overlap <- protein1.ecf.hits[1]
  ecf.overlap.parms <- calculate_overlap_parms(aln1.start, aln1.end, protein1.ecf.best.overlap$qstart, protein1.ecf.best.overlap$qend)
  ecf.aln.length <- protein1.ecf.best.overlap$qend - protein1.ecf.best.overlap$qstart + 1
  prop.of.ecf.hit.that.overlaps.with.rec <- ecf.overlap.parms$overlap.length/ecf.aln.length
  prop.of.full.ecf.covered.by.overlap <- prop.of.ecf.hit.that.overlaps.with.rec*protein1.ecf.best.overlap$scov
  protein1.ecf.score <- prop.of.full.ecf.covered.by.overlap
  
  protein2.length <- this.rr$slength
  aln2.start <- this.rr$sstart
  aln2.end <- this.rr$send
  aln2.length <- aln2.end - aln2.start + 1
  
  protein2.ecod.hits <- reprseq.ecod[.(protein2.name)]
  protein2.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln2.start, aln2.end, qstart, qend)]
  setorder(protein2.ecod.hits, -rr.overlap.jaccard)
  protein2.ecod.best.overlap <- protein2.ecod.hits[1]
  ecod.overlap.parms <- calculate_overlap_parms(aln2.start, aln2.end, protein2.ecod.best.overlap$qstart, protein2.ecod.best.overlap$qend)
  ecod.aln.length <- protein2.ecod.best.overlap$qend - protein2.ecod.best.overlap$qstart + 1
  prop.of.ecod.hit.that.overlaps.with.rec <- ecod.overlap.parms$overlap.length/ecod.aln.length
  prop.of.full.ecod.covered.by.overlap <- prop.of.ecod.hit.that.overlaps.with.rec*protein2.ecod.best.overlap$scov
  protein2.ecod.score <- prop.of.full.ecod.covered.by.overlap
  
  protein2.ecf.hits <- reprseq.ecfs[.(protein2.name)]
  protein2.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln2.start, aln2.end, qstart, qend)]
  setorder(protein2.ecf.hits, -rr.overlap.jaccard)
  protein2.ecf.best.overlap <- protein2.ecf.hits[1]
  ecf.overlap.parms <- calculate_overlap_parms(aln2.start, aln2.end, protein2.ecf.best.overlap$qstart, protein2.ecf.best.overlap$qend)
  ecf.aln.length <- protein2.ecf.best.overlap$qend - protein2.ecf.best.overlap$qstart + 1
  prop.of.ecf.hit.that.overlaps.with.rec <- ecf.overlap.parms$overlap.length/ecf.aln.length
  prop.of.full.ecf.covered.by.overlap <- prop.of.ecf.hit.that.overlaps.with.rec*protein2.ecf.best.overlap$scov
  protein2.ecf.score <- prop.of.full.ecf.covered.by.overlap
  
  if(is.na(protein1.ecod.score)) protein1.ecod.score <- 0
  if(is.na(protein1.ecf.score)) protein1.ecf.score <- 0
  if(is.na(protein2.ecod.score)) protein2.ecod.score <- 0
  if(is.na(protein2.ecf.score)) protein2.ecf.score <- 0
  
  rbind(data.table(pair.index, reprseq = protein1.name, ecod.score = protein1.ecod.score, ecf.score = protein1.ecf.score), 
        data.table(pair.index, reprseq = protein2.name, ecod.score = protein2.ecod.score, ecf.score = protein2.ecf.score))
  
})
overlap.scores <- rbindlist(overlap.scores.list)

# Reshape data to long format for faceting
overlap.scores.long <- melt(overlap.scores, 
                            id.vars = c("pair.index", "reprseq"), 
                            measure.vars = c("ecod.score", "ecf.score"), 
                            variable.name = "Score_Type", 
                            value.name = "Score_Value")

# Plot histograms with facet_grid
ggplot(overlap.scores.long, aes(x = Score_Value, fill = Score_Type)) +
  geom_histogram(binwidth = 0.05, color = "black", alpha = 0.7) +
  facet_grid(rows = vars(Score_Type), scales = "free_y") +  # Free y-scales for better visualization
  labs(title = "Distribution of ECOD and ECF Scores",
       x = "Score Value", y = "Frequency") +
  theme_minimal()



#############################
### Relative overlap distribution
#############################

rel.pos <- function(pos, length) (pos-1)/(length-1)
rr <- recent.rec.table[pident>=95]
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

# Apply function row-wise and store as lists
pairs.rel.pos[, prot1_bins := lapply(seq_len(.N), function(i) 
  assign_bins(prot1.rel.start[i], prot1.rel.end[i]))]

pairs.rel.pos[, prot2_bins := lapply(seq_len(.N), function(i) 
  assign_bins(prot2.rel.start[i], prot2.rel.end[i]))]

# Expand `prot1_bins` and `prot2_bins` while ensuring equal row repetition
pairs.rel.pos.expanded <- pairs.rel.pos[
  , .(protein1.name = rep(protein1.name, times = lengths(prot1_bins)),
      protein2.name = rep(protein2.name, times = lengths(prot2_bins)),
      prot1_bin = unlist(prot1_bins),
      prot2_bin = unlist(prot2_bins))
]

# Count occurrences of each (prot1_bin, prot2_bin) combination
heatmap_data <- pairs.rel.pos.expanded[, .N, by = .(prot1_bin, prot2_bin)]


# heatmap_data$prot1_bin <- factor(heatmap_data$prot1_bin, levels = )
ggplot(heatmap_data, aes(x = prot1_bin, y = prot2_bin, fill = N)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Heatmap of Overlapping Regions",
       x = "Protein 1 Bins",
       y = "Protein 2 Bins",
       fill = "Count")


#############################
### N/C/other assignment
#############################

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


# rel.pos <- function(pos, length) (pos-1)/(length-1)
rr <- recent.rec.table[pident >= 95]
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

#############################
### ECFs associated with overlaps
#############################

# rr <- recent.rec.table[pident >= 95]
rr <- recent.rec.table
overlapping.ecf.pairs.list <- lapply(1:nrow(rr), function(pair.index){
  this.rr <- rr[pair.index]
  protein1.name <- this.rr$qname
  protein2.name <- this.rr$sname
  
  protein1.length <- this.rr$qlength
  aln1.start <- this.rr$qstart
  aln1.end <- this.rr$qend
  aln1.length <- aln1.end - aln1.start + 1
  
  protein1.ecf.hits <- reprseq.ecfs[.(protein1.name)]
  protein1.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
  setorder(protein1.ecf.hits, -rr.overlap.jaccard)
  protein1.ecf.best.overlap <- protein1.ecf.hits[1]
  ecf.overlap.parms <- calculate_overlap_parms(aln1.start, aln1.end, protein1.ecf.best.overlap$qstart, protein1.ecf.best.overlap$qend)
  ecf.aln.length <- protein1.ecf.best.overlap$qend - protein1.ecf.best.overlap$qstart + 1
  prop.of.ecf.hit.that.overlaps.with.rec <- ecf.overlap.parms$overlap.length/ecf.aln.length
  prop.of.full.ecf.covered.by.overlap <- prop.of.ecf.hit.that.overlaps.with.rec*protein1.ecf.best.overlap$scov
  protein1.ecf.score <- prop.of.full.ecf.covered.by.overlap
  protein1.best.overlapping.ecf <- protein1.ecf.best.overlap$sname
  protein1.best.overlapping.ef <- ecf.family.table[.(protein1.best.overlapping.ecf)]$family
  
  protein2.length <- this.rr$slength
  aln2.start <- this.rr$sstart
  aln2.end <- this.rr$send
  aln2.length <- aln2.end - aln2.start + 1
  
  protein2.ecf.hits <- reprseq.ecfs[.(protein2.name)]
  protein2.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln2.start, aln2.end, qstart, qend)]
  setorder(protein2.ecf.hits, -rr.overlap.jaccard)
  protein2.ecf.best.overlap <- protein2.ecf.hits[1]
  ecf.overlap.parms <- calculate_overlap_parms(aln2.start, aln2.end, protein2.ecf.best.overlap$qstart, protein2.ecf.best.overlap$qend)
  ecf.aln.length <- protein2.ecf.best.overlap$qend - protein2.ecf.best.overlap$qstart + 1
  prop.of.ecf.hit.that.overlaps.with.rec <- ecf.overlap.parms$overlap.length/ecf.aln.length
  prop.of.full.ecf.covered.by.overlap <- prop.of.ecf.hit.that.overlaps.with.rec*protein2.ecf.best.overlap$scov
  protein2.ecf.score <- prop.of.full.ecf.covered.by.overlap
  protein2.best.overlapping.ecf <- protein2.ecf.best.overlap$sname
  protein2.best.overlapping.ef <- ecf.family.table[.(protein2.best.overlapping.ecf)]$family
  
  data.table(protein1.name, protein2.name,
             protein1.ecf.score, protein1.best.overlapping.ecf, protein1.best.overlapping.ef,
             protein2.ecf.score, protein2.best.overlapping.ecf, protein2.best.overlapping.ef)
})
overlapping.ecf.pairs <- rbindlist(overlapping.ecf.pairs.list)

recombination.ecfs <- overlapping.ecf.pairs[protein1.ecf.score >= 0.8 | protein2.ecf.score > 0.8]
rec.ecfs.prot1 <- overlapping.ecf.pairs[protein1.ecf.score >= 0.9]$protein1.best.overlapping.ecf

rec.ecfs.prot1.metadata <- ecf.class.table[.(rec.ecfs.prot1)]








# ####################################################################################
# ####################################################################################
# ####################################################################################
# ####################################################################################
# 
# 
# 
# 
# rrr <- recent.rec.table[pident >= 95]
# no.pairs <- nrow(rrr)
# 
# # 1: visualse direct comparison of the two proteins
# for(pair.index in 1:no.pairs){
#   this.pair.output.dir <- file.path(output.dir, sprintf("pair%.*d", nchar(no.pairs), pair.index))
#   make.dir(this.pair.output.dir)
#   protein1_name <- rrr$qname[pair.index]
#   protein2_name <- rrr$sname[pair.index]
#   protein1_overlap <- max.overlap.per.reprseq[.(protein1_name)]$ncbi.overlap
#   protein2_overlap <- max.overlap.per.reprseq[.(protein2_name)]$ncbi.overlap
#   protein1_length <- rrr$qlength[pair.index]
#   protein2_length <- rrr$slength[pair.index]
#   aln1_start <- rrr$qstart[pair.index]
#   aln1_end <- rrr$qend[pair.index]
#   aln2_start <- rrr$sstart[pair.index]
#   aln2_end <- rrr$send[pair.index]
#   aln_identity <- rrr$pident[pair.index]
#   protein1_description <- sprintf("F: %s, O: %.2f", reprseq.functions[.(protein1_name)]$funct, protein1_overlap)
#   protein2_description <- sprintf("F: %s, O: %.2f", reprseq.functions[.(protein2_name)]$funct, protein2_overlap)
#   
#   p.align <- plot_protein_alignment(protein1_name, protein2_name,
#                                     protein1_length, protein2_length,
#                                     aln1_start, aln1_end, aln2_start, aln2_end,
#                                     ecf1_start, ecf1_end, ecf2_start, ecf2_end,
#                                     aln_identity, protein1_description, protein2_description)
#   
#   
#   # 2: visualise domains within each protein
#   ### get ECOD hits
#   protein1.ecod.hits <- reprseq.ecod[.(protein1_name)]
#   protein2.ecod.hits <- reprseq.ecod[.(protein2_name)]
#   # 
#   protein1.ecod.hits.best <- select.best.ecod.hits(protein1.ecod.hits)
#   protein2.ecod.hits.best <- select.best.ecod.hits(protein2.ecod.hits)
#   
#   ### get ECF hits
#   protein1.ecf.hits <- reprseq.ecfs[.(protein1_name)]
#   protein2.ecf.hits <- reprseq.ecfs[.(protein2_name)]
#   # 
#   protein1.ecf.hits.all <- select.all.ecf.hits(protein1.ecf.hits)
#   protein2.ecf.hits.all <- select.all.ecf.hits(protein2.ecf.hits)
#   # 
#   protein1.ecf.hits.best <- select.best.ecf.hits(protein1.ecf.hits)
#   protein2.ecf.hits.best <- select.best.ecf.hits(protein2.ecf.hits)
#   
#   
#   # Combine tables
#   protein1.domain.hits <- combine_domain_hits(protein1.ecod.hits.best, protein1.ecf.hits.best, types = c("ECOD", "ECF"))
#   protein2.domain.hits <- combine_domain_hits(protein2.ecod.hits.best, protein2.ecf.hits.best, types = c("ECOD", "ECF"))
#   
#   p.ecfs1 <- plot_ecf_domains_with_alignment(protein1_name, protein1_length, aln1_start, aln1_end, protein1.ecf.hits.all)
#   p.ecfs2 <- plot_ecf_domains_with_alignment(protein2_name, protein2_length, aln2_start, aln2_end, protein2.ecf.hits.all)
#   
#   p.domains1 <- plot_protein_domains_compact(protein1_length, aln1_start, aln1_end, protein1.domain.hits, protein1_name)
#   p.domains2 <- plot_protein_domains_compact(protein2_length, aln2_start, aln2_end, protein2.domain.hits, protein2_name)
#   
#   p.align.fn <- file.path(this.pair.output.dir, "alignment.png")
#   p.ecfs1.fn <- file.path(this.pair.output.dir, "ecfs1.png")
#   p.ecfs2.fn <- file.path(this.pair.output.dir, "ecfs2.png")
#   p.domains1.fn <- file.path(this.pair.output.dir, "domains1.png")
#   p.domains2.fn <- file.path(this.pair.output.dir, "domains2.png")
#   ggsave(p.align.fn, p.align, width = 6, height = 4)
#   if(!is.null(p.ecfs1)) ggsave(p.ecfs1.fn, p.ecfs1, width = 6, height = 4)
#   if(!is.null(p.ecfs2)) ggsave(p.ecfs2.fn, p.ecfs2, width = 6, height = 4)
#   if(!is.null(p.domains1)) ggsave(p.domains1.fn, p.domains1, width = 6, height = 4)
#   if(!is.null(p.domains2)) ggsave(p.domains2.fn, p.domains2, width = 6, height = 4)
# }
# 
