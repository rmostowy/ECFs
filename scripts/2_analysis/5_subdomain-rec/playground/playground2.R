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

output.dir <- "vis-rec-pairs"
make.dir(output.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.ncbi.filtered.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load ECF families
ecf.family.table <- fread(ecf.family.table.output)
setkey(ecf.family.table, ecf)

# Load ECF table
ecf.class.table <- fread(ecf.classification.table.output)
ecf.class.table$family <- ecf.family.table[.(ecf.class.table$ecf)]$family
ecf2ecod <- ecf.class.table %>% select(ecf, family, t.id)
ecf2ecod$t.name <- sapply(1:nrow(ecf2ecod), function(k){
  out <- "unknown"
  this.t.id <- ecf2ecod$t.id[k]
  if(this.t.id != "unknown") out <- ecod.domains[t_id == this.t.id]$t_name[1]  
  return(out)
})
setorder(ecf2ecod, family)


# Load all-by-all reprseq
reprseq.vs.reprseq <- fread(reprseq.vs.reprseq.filename)
reprseq.vs.reprseq[, qcov := round((qend - qstart + 1)/qlength, d = 4)]
reprseq.vs.reprseq[, scov := round((send - sstart + 1)/slength, d = 4)]
reprseq.vs.reprseq <- reprseq.vs.reprseq[qname %chin% all.reprseq & sname %chin% all.reprseq]

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions.unique, reprseq)

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
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)

# Load reprseq qcov/scov
# IMPORTANT: loads a dataset outside of config;
# if used in real analysis, it needs to be solved

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


# find hits with the highest pident for all reprseq pairs
reprseq.vs.reprseq.prefiltered <- reprseq.vs.reprseq[prob >= 95]
setorder(reprseq.vs.reprseq.prefiltered, qname, sname, -pident)
reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq.prefiltered[, .SD[1], by = .(qname, sname)]
reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq.filtered[qname < sname]
reprseq.vs.reprseq.filtered.qcov.scov <- reprseq.cov[.(reprseq.vs.reprseq.filtered[, .(qname, sname)])]
reprseq.vs.reprseq.filtered.qcov.scov[, qcov := pmin(qcov_min, qcov_max, na.rm = T)]
reprseq.vs.reprseq.filtered.qcov.scov[, scov := pmin(scov_min, scov_max, na.rm = T)]
rr <- cbind(reprseq.vs.reprseq.filtered, reprseq.vs.reprseq.filtered.qcov.scov %>% select(pident.full = pident, qcov.full = qcov, scov.full = scov))
rr[, pident.full := pident.full*100]
rr[, hit.length := qend - qstart + 1]

rr[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

ggplot(rr, aes(y = pident, x = cov.max.full)) + geom_bin2d(bins = 100)

rrr <- rr[pident > 90 & cov.max.full < 0.5 & hit.length >= 30]
# rrr <- rr[pident > 60 & cov.max.full < 0.8 & hit.length >= 10]
no.pairs <- nrow(rrr)
setkey(rrr, qname, sname)

# all.shared <- unique(c(rrr$qname, rrr$sname))
# my.table(reprseq.functions[.(all.shared)]$funct)[1:30]


# VISUALISATION

# get annotations
reprseq.functions <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions, reprseq)


# 1: visualse direct comparison of the two proteins
for(pair.index in 1:no.pairs){
  this.pair.output.dir <- file.path(output.dir, sprintf("pair%.*d", nchar(no.pairs), pair.index))
  make.dir(this.pair.output.dir)
  protein1_name <- rrr$qname[pair.index]
  protein2_name <- rrr$sname[pair.index]
  protein1_length <- rrr$qlength[pair.index]
  protein2_length <- rrr$slength[pair.index]
  aln1_start <- rrr$qstart[pair.index]
  aln1_end <- rrr$qend[pair.index]
  aln2_start <- rrr$sstart[pair.index]
  aln2_end <- rrr$send[pair.index]
  aln_identity <- rrr$pident[pair.index]
  protein1_description <- sprintf("F: %s", reprseq.functions[.(protein1_name)]$funct)
  protein2_description <- sprintf("F: %s", reprseq.functions[.(protein2_name)]$funct)
  
  p.align <- plot_protein_alignment(protein1_name, protein2_name,
                                    protein1_length, protein2_length,
                                    aln1_start, aln1_end, aln2_start, aln2_end,
                                    ecf1_start, ecf1_end, ecf2_start, ecf2_end,
                                    aln_identity, protein1_description, protein2_description)
  
  
  # 2: visualise domains within each protein
  ### get ECOD hits
  protein1.ecod.hits <- reprseq.ecod[.(protein1_name)]
  protein2.ecod.hits <- reprseq.ecod[.(protein2_name)]
  # 
  protein1.ecod.hits.best <- select.best.ecod.hits(protein1.ecod.hits)
  protein2.ecod.hits.best <- select.best.ecod.hits(protein2.ecod.hits)
  
  ### get ECF hits
  protein1.ecf.hits <- reprseq.ecfs[.(protein1_name)]
  protein2.ecf.hits <- reprseq.ecfs[.(protein2_name)]
  # 
  protein1.ecf.hits.all <- select.all.ecf.hits(protein1.ecf.hits)
  protein2.ecf.hits.all <- select.all.ecf.hits(protein2.ecf.hits)
  # 
  protein1.ecf.hits.best <- select.best.ecf.hits(protein1.ecf.hits)
  protein2.ecf.hits.best <- select.best.ecf.hits(protein2.ecf.hits)
  
  
  # Combine tables
  protein1.domain.hits <- combine_domain_hits(protein1.ecod.hits.best, protein1.ecf.hits.best, types = c("ECOD", "ECF"))
  protein2.domain.hits <- combine_domain_hits(protein2.ecod.hits.best, protein2.ecf.hits.best, types = c("ECOD", "ECF"))
  
  p.ecfs1 <- plot_ecf_domains_with_alignment(protein1_name, protein1_length, aln1_start, aln1_end, protein1.ecf.hits.all)
  p.ecfs2 <- plot_ecf_domains_with_alignment(protein2_name, protein2_length, aln2_start, aln2_end, protein2.ecf.hits.all)
  
  p.domains1 <- plot_protein_domains_compact(protein1_length, protein1.domain.hits, protein1_name)
  p.domains2 <- plot_protein_domains_compact(protein2_length, protein2.domain.hits, protein2_name)
  
  p.align.fn <- file.path(this.pair.output.dir, "alignment.png")
  p.ecfs1.fn <- file.path(this.pair.output.dir, "ecfs1.png")
  p.ecfs2.fn <- file.path(this.pair.output.dir, "ecfs2.png")
  p.domains1.fn <- file.path(this.pair.output.dir, "domains1.png")
  p.domains2.fn <- file.path(this.pair.output.dir, "domains2.png")
  ggsave(p.align.fn, p.align, width = 6, height = 4)
  if(!is.null(p.ecfs1)) ggsave(p.ecfs1.fn, p.ecfs1, width = 6, height = 4)
  if(!is.null(p.ecfs2)) ggsave(p.ecfs2.fn, p.ecfs2, width = 6, height = 4)
  if(!is.null(p.domains1)) ggsave(p.domains1.fn, p.domains1, width = 6, height = 4)
  if(!is.null(p.domains2)) ggsave(p.domains2.fn, p.domains2, width = 6, height = 4)
}

###############################################
################## FUNCTIONS ##################
###############################################


# Funciton to visulaise the protein alignment
plot_protein_alignment <- function(protein1_name, protein2_name,
                                   protein1_length, protein2_length, 
                                   aln1_start, aln1_end, aln2_start, aln2_end, 
                                   ecf1_start, ecf1_end, ecf2_start, ecf2_end, 
                                   aln_identity, protein1_description, protein2_description) {
  
  # Customization options for the alignment polygon) {
  
  # Customization options for the alignment polygon
  aln_colour <- "lightblue"  # Color of the alignment polygon
  domain_colour <- "red"     # Color of the domain regions
  
  # Create a base plot with horizontal lines for proteins
  p <- ggplot() +
    # Draw Protein 1 and Protein 2 as black horizontal lines
    geom_segment(aes(x = 0, xend = protein1_length, y = 2, yend = 2), linewidth = 2, color = "black") +
    geom_segment(aes(x = 0, xend = protein2_length, y = 1, yend = 1), linewidth = 2, color = "black")
  
  # Add the connecting rectangle for the alignment
  p <- p + geom_polygon(data = data.frame(
    x = c(aln1_start, aln1_end, aln2_end, aln2_start),
    y = c(2, 2, 1, 1)
  ), aes(x = x, y = y), fill = aln_colour, alpha = aln_identity) +
    # Add text annotations for protein descriptions
    annotate("text", x = 0, y = 1.03, label = protein2_description, hjust = 0, size = 5, vjust = 0) +
    annotate("text", x = 0, y = 1.97, label = protein1_description, hjust = 0, size = 5, vjust = 1)
  
  # Customize plot
  p <- p +
    labs(x = "Position", y = "") +
    scale_y_continuous(breaks = c(2, 1), labels = c(protein1_name, protein2_name)) +
    theme(axis.text.y = element_text(size = 12),
          plot.margin = margin(20, 30, 20, 30))  # Adjust margins to fit text
  
  return(p)
}

select.best.ecod.hits <- function(ecod.hits){
  ecod.hits.out <- NULL
  all.ecod.names <- ecod.hits$sname
  all.ecod.names <- all.ecod.names[!is.na(all.ecod.names)]
  n <- length(all.ecod.names)
  if(n > 0){
    setorder(ecod.hits, qname, t.id, -bitscore)
    ecod.hits <- ecod.hits[, .SD[1], by = .(qname, t.id)]  
    ecod.hits.out <- ecod.hits %>% select(qstart, qend, name = t.id)
  }
  return(ecod.hits.out)
}

select.all.ecf.hits <- function(ecf.hits){
  ecf.hits.out <- NULL
  all.ecf.names <- ecf.hits$sname
  all.ecf.names <- all.ecf.names[!is.na(all.ecf.names)]
  n <- length(all.ecf.names)
  if(n > 0) ecf.hits.out <- ecf.hits %>% select(qstart, qend, name = sname, family = EF)
  return(ecf.hits.out)
}

select.best.ecf.hits <- function(ecf.hits){
  ecf.hits.out <- NULL
  all.ecf.names <- ecf.hits$sname
  all.ecf.names <- all.ecf.names[!is.na(all.ecf.names)]
  n <- length(all.ecf.names)
  if(n > 0){
    setorder(ecf.hits, qname, EF, -bitscore)
    ecf.hits <- ecf.hits[, .SD[1], by = .(qname, EF)]  
    ecf.hits[, ecf.name := sprintf("%s:%s", EF, sname)]
    ecf.hits.out <- ecf.hits %>% select(qstart, qend, name = ecf.name)
  }
  return(ecf.hits.out)
}

# Function to combine domain hits
combine_domain_hits <- function(..., types) {
  # List of input tables
  input_tables <- list(...)
  
  # Check for NULL input
  if (all(sapply(input_tables, is.null))) {
    return(NULL)  # Return NULL if all input tables are NULL
  }
  
  # Combine non-NULL tables with their types
  rbindlist(
    lapply(seq_along(input_tables), function(i) {
      if (!is.null(input_tables[[i]])) {
        data.table(type = types[i], input_tables[[i]])
      } else {
        NULL
      }
    }),
    use.names = TRUE, fill = TRUE
  )
}

# Function to visualise ECF domains and alignment on the query
plot_ecf_domains_with_alignment <- function(protein_name, protein_length, aln1_start, aln1_end, protein_ecf_hits_all) {
  
  # Check if domain_data is NULL or empty
  if (is.null(protein_ecf_hits_all) || nrow(protein_ecf_hits_all) == 0) {
    message("No domains detected for this protein.")
    return(NULL)
  }
  
  # Add a y-axis level for each ECF domain, reducing vertical spread for clarity
  protein_ecf_hits_all$y_level <- seq_len(nrow(protein_ecf_hits_all)) * 0.5
  
  # Define colors for families
  family_colors <- unique(protein_ecf_hits_all$family)
  color_palette <- setNames(RColorBrewer::brewer.pal(length(family_colors), "Set3"), family_colors)
  
  # Create the plot
  p <- ggplot() +
    # Plot the query protein as a black line with the aligned segment in red
    geom_segment(aes(x = 0, xend = protein_length, y = 0, yend = 0), linewidth = 2, color = "black") +
    geom_segment(aes(x = aln1_start, xend = aln1_end, y = 0, yend = 0), linewidth = 2, color = "red") +
    # Plot each ECF domain as a coloured segment, positioned on its y_level
    geom_segment(data = protein_ecf_hits_all, aes(x = qstart, xend = qend, y = y_level, yend = y_level, color = family),
                 linewidth = 4) +
    # Add labels for each ECF domain
    geom_text(data = protein_ecf_hits_all, aes(x = (qstart + qend) / 2, y = y_level + 0.1, label = name),
              size = 3, hjust = 0.5, color = "black") +
    # Define custom colours for families
    scale_color_manual(values = color_palette) +
    # Customise the plot
    labs(
      x = sprintf("Position along protein: %s", protein_name),
      y = "",
      title = "ECF Domains and Query Alignment",
      color = "Family"
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(20, 30, 20, 30)
    ) +
    coord_cartesian(xlim = c(0, protein_length))
  
  return(p)
}

# Function to visualise protein domains with protein name in x-axis label
plot_protein_domains_compact <- function(protein_length, domain_data, protein_name) {
  
  # Check if domain_data is NULL or empty
  if (is.null(domain_data) || nrow(domain_data) == 0) {
    message("No domains detected for this protein.")
    return(NULL)
  }
  
  # Add a y-axis level for each domain, with reduced spacing for compactness
  domain_data$y_level <- seq_len(nrow(domain_data)) * 0.5  # Adjust spacing
  
  # Define colours for domain types
  domain_colours <- c("ECOD" = "red", "ECF" = "green")
  
  # Create the plot
  p <- ggplot(data = domain_data) +
    # Plot the full protein as a horizontal black line at the bottom
    geom_segment(aes(x = 0, xend = protein_length, y = 0, yend = 0),
                 linewidth = 2, color = "black") +
    # Plot each domain as a segment coloured by its type
    geom_segment(aes(x = qstart, xend = qend, y = y_level, yend = y_level, color = type),
                 linewidth = 4) +
    # Add text labels for each domain's name above the segment
    geom_text(aes(x = qstart, y = y_level + 0.1, label = name),
              size = 4, hjust = 0, vjust = 0, color = "blue") +
    # Define custom colours for domain types
    scale_color_manual(values = domain_colours) +
    # Customise the plot
    labs(
      x = sprintf("%s", protein_name),  # Add protein name to x-axis
      y = "",
      color = "Domain Type"
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(20, 30, 20, 30)
    ) +
    coord_cartesian(xlim = c(0, protein_length))  # Ensure x-axis matches protein length
  
  return(p)
}


# a <- "reprseq133624"
# b <- "reprseq120726"
# reprseq.vs.reprseq.filtered[qname == a & sname == b]
# reprseq.vs.reprseq.filtered[qname == b & sname == a]


# # Load taxonomy data
# reprseq2taxonomy <- unique(reprseq.table[, c("reprseq", "cluster.name")])
# reprseq2taxonomy[, accession := sub("\\..*", "", cluster.name)]
# reprseq2taxonomy[, cluster.name := NULL]
# metadata.table <- fread(metadata.table.output)
# setkey(metadata.table, accession)
# 
# reprseq2taxonomy$family.ncbi <- metadata.table[.(reprseq2taxonomy$accession)]$family.ncbi
# reprseq2taxonomy$genus.ncbi <- metadata.table[.(reprseq2taxonomy$accession)]$genus.ncbi
# setkey(reprseq2taxonomy, reprseq)


# # ############
# # 1. Find all overlapping ECFs and ECODs in reprseqs and calculate overlaps
# # ############
# 
# # Load reprseq vs ECFs
# reprseq.ecfs.all <- fread(reprseq.ecfs.output)
# reprseq.ecfs <- reprseq.ecfs.all %>%
#   filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
# 
# # Load reprseq vs ECOD
# reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
# reprseq.ecod <- reprseq.ecod.all %>%
#   filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE)
# 
# # Calculate overlap between all ECFs and ECODs in all reprseq
# reprseq.ecfs.overlap <- reprseq.ecfs %>% select(qname, qlength, qstart.ecf = qstart , qend.ecf = qend,
#                                                 sname.ecf = sname, sstart.ecf = sstart, send.ecf = send, slength.ecf = slength,
#                                                 pident.ecf = pident, hit.length.ecf = hit.length)
# reprseq.ecods.overlap <- reprseq.ecod %>% select(qname, qstart.ecod = qstart , qend.ecod = qend,
#                                                  sname.ecod = sname, sstart.ecod = sstart, send.ecod = send, slength.ecod = slength,
#                                                  pident.ecod = pident, hit.length.ecod = hit.length, bitscore.ecod = bitscore)
# combined.ecfs.ecods <- reprseq.ecfs.overlap[reprseq.ecods.overlap, on = .(qname), allow.cartesian = TRUE]
# combined.ecfs.ecods[, c("overlap.start", "overlap.end", "overlap.length") := {
#   res <- calculate_overlap_parms(qstart.ecf, qend.ecf, qstart.ecod, qend.ecod)
#   list(res$overlap.start, res$overlap.end, res$overlap.length)
# }]
# combined.ecfs.ecods <- na.omit(combined.ecfs.ecods, cols = "overlap.length")[overlap.length > 0]
# 
# # Filter by minimum length
# reprseq.ecf.ecod.overlaps <- combined.ecfs.ecods %>% filter(overlap.length >= MIN.OVERLAP.LENGTH)
# setkey(reprseq.ecf.ecod.overlaps, sname.ecf)
# 
# # Calculate overlaps between ECFs and ECODs:
# 
# # a) 'ecf.cov.overlap' is the proportion of the ECF hit (from the query's point of view) 
# # that is covered by the overlap (if high, then most of the ECF is covered by the ECF-ECOD)
# reprseq.ecf.ecod.overlaps[, ecf.cov.overlap := round(overlap.length/hit.length.ecf, d = 4)]
# 
# # b) 'ecod.cov.overlap' is the proportion of the ECOD hit (not entire ECOD) covered by the overlap
# # (again, if ~1, then the ECOD hit is equivalent to the overlap)
# reprseq.ecf.ecod.overlaps[, ecod.cov.overlap := round(overlap.length/hit.length.ecod, d = 4)]
# 
# # c) 'scov.ecod' is the proportion of the ECOD length that hits the query
# reprseq.ecf.ecod.overlaps[, scov.ecod := round((send.ecod - sstart.ecod + 1)/slength.ecod, d = 4)]
# 
# # d) 'ecod.cov.by.ecf' is the proportion of the entire ECOD covered by the overlap,
# # and it is calculated by weighing the proportion of full ECOD's length that hits the query
# # by the proportion of that hit length that overlaps with the ECF.
# reprseq.ecf.ecod.overlaps[, ecod.cov.by.ecf := round(scov.ecod * ecod.cov.overlap, d = 4)]
# 
# # Consider only ECFs which are nearly covered by ECODs
# reprseq.ecf.ecod.overlaps.shortlist <- reprseq.ecf.ecod.overlaps[reprseq.ecf.ecod.overlaps$ecf.cov.overlap > MIN.QCOV.FOR.DOMAIN.EQUIVALENCE]
# 
# # ############
# # 2. Identify all pairs of ECFs which connect ECOF
# # ############
# 
# # Identify all ECFs which overlap with at least 2 different H level ECODs
# ecf.ecod.overlapping <- reprseq.ecf.ecod.overlaps.shortlist %>% select(ecf = sname.ecf, ecod = sname.ecod)
# ecf.ecod.overlapping$x_id <- ecod.domains[.(ecf.ecod.overlapping$ecod)]$x_id
# ecf.ecod.overlapping$h_id <- ecod.domains[.(ecf.ecod.overlapping$ecod)]$h_id
# ecf.ecod.overlapping$t_id <- ecod.domains[.(ecf.ecod.overlapping$ecod)]$t_id
# no.t.per.ecf <- ecf.ecod.overlapping[, .(no.t.levels = uniqueN(t_id)), by = .(ecf)]
# ecfs.with.multiple.t.levels <- sort(no.t.per.ecf[no.t.levels>1]$ecf)
# 
# # Of all reprseq-ecf pairs, consider only a single ECOD hit amongst those with the same H level (with the highest bitscore)
# reprseq.ecf.ecod.overlaps.filtered1 <- reprseq.ecf.ecod.overlaps.shortlist[sname.ecf %in% ecfs.with.multiple.t.levels]
# reprseq.ecf.ecod.overlaps.filtered1$x.id <- ecod.domains[.(reprseq.ecf.ecod.overlaps.filtered1$sname.ecod)]$x_id
# reprseq.ecf.ecod.overlaps.filtered1$h.id <- ecod.domains[.(reprseq.ecf.ecod.overlaps.filtered1$sname.ecod)]$h_id
# reprseq.ecf.ecod.overlaps.filtered1$t.id <- ecod.domains[.(reprseq.ecf.ecod.overlaps.filtered1$sname.ecod)]$t_id
# setorder(reprseq.ecf.ecod.overlaps.filtered1, qname, sname.ecf, t.id, -bitscore.ecod)
# reprseq.ecf.ecod.overlaps.filtered2 <- reprseq.ecf.ecod.overlaps.filtered1[, .SD[1], by = .(qname, sname.ecf, t.id)]
# 
# # Inner-join the table to list of pairs of reprseqs
# left.join.table <- reprseq.ecf.ecod.overlaps.filtered2 %>% select(query.a = qname, qlength.a = qlength, shared.ecf = sname.ecf, ecf.length = slength.ecf,
#                                          ecf.start.a = qstart.ecf, ecf.end.a = qend.ecf, 
#                                          ecod.a = sname.ecod, x.id.a = x.id, h.id.a = h.id, t.id.a = t.id, ecod.start.a = qstart.ecod, ecod.end.a = qend.ecod, ecod.cov.by.ecf.a = ecod.cov.by.ecf)
# 
# right.join.table <- reprseq.ecf.ecod.overlaps.filtered2 %>% select(query.b = qname, qlength.b = qlength, shared.ecf = sname.ecf,
#                                           ecf.start.b = qstart.ecf, ecf.end.b = qend.ecf, 
#                                           ecod.b = sname.ecod, x.id.b = x.id, h.id.b = h.id, t.id.b = t.id, ecod.start.b = qstart.ecod, ecod.end.b = qend.ecod, ecod.cov.by.ecf.b = ecod.cov.by.ecf)
# 
# putative.recomb.pairs.temp <- left.join.table[right.join.table, on = .(shared.ecf), allow.cartesian = TRUE]
# putative.recomb.pairs.temp <- putative.recomb.pairs.temp[query.a < query.b]
# # putative.recomb.pairs.temp <- putative.recomb.pairs.temp[query.a < query.b & t.id.a != t.id.b]
# 
# putative.recomb.pairs <- putative.recomb.pairs.temp %>% 
#   select(query.a, query.b, shared.ecf, qlength.a, qlength.b,
#          ecf.start.a, ecf.end.a, ecf.start.b, ecf.end.b,
#          ecod.a, ecod.b, x.id.a, h.id.a, x.id.b, h.id.b, t.id.a, t.id.b,
#          ecod.start.a, ecod.end.a, ecod.start.b, ecod.end.b,
#          ecod.cov.by.ecf.a, ecod.cov.by.ecf.b)
# # putative.recomb.pairs[ , different.x := x.id.a != x.id.b]
# 
# 
# # ############
# # 3. Merge-in the information about reprseq alignment
# # ############
# 
# # Filter 'reprseq.vs.reprseq' to leave only the strongest hit between each pair
# # setorder(reprseq.vs.reprseq, qname, sname, -bitscore)
# setorder(reprseq.vs.reprseq, qname, sname, -pident)
# reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq[, .SD[1], by = .(qname, sname)] %>%
#   select(qname, sname, aln.start.a = qstart, aln.end.a = qend, aln.start.b = sstart, aln.end.b = send,
#          cov.a = qcov, cov.b = scov, aln.prob = prob, aln.eval = eval, aln.pident = pident)
# 
# # Merge full reprseq alignments to add to the 'putative.recomb.pairs' table
# putative.recomb.pairs.aln <- merge(putative.recomb.pairs, reprseq.vs.reprseq.filtered,
#                                       by.x = c("query.a", "query.b"),
#                                       by.y = c("qname", "sname"),
#                                       all.x = TRUE, all.y = FALSE)
# putative.recomb.pairs.aln <- putative.recomb.pairs.aln[!is.na(aln.pident)]
# putative.recomb.pairs.aln[, max.cov := pmax(cov.a, cov.b, na.rm = T)]
# 
# ggplot(putative.recomb.pairs.aln[aln.pident>20], aes(x = max.cov, y = aln.pident)) + geom_bin2d()
# 
# hist(putative.recomb.pairs.aln$aln.pident)
# sh <- putative.recomb.pairs.aln[max.cov < 0.5 & aln.pident > 70]
# 
# 
# 
# setorder(reprseq.vs.reprseq, qname, sname, -pident)
# rr <- reprseq.vs.reprseq[, .SD[1], by = .(qname, sname)]
# rrr <- rr[prob>95 & pident > 95 & max.cov < 0.3]
# 
# 
