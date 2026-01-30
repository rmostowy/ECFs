### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(gggenomes)))


### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

NO.GENES.TO.EXTEND <- 4L
PAD.PARM <- 10L

ecf.subseqs.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs", "ecf-subseqs")

# Load data
metadata.table <- fread(metadata.table.output)
refseq.protein.stats.table <- fread(refseq.protein.stats.table.filename)
name2row <- setNames(seq_len(nrow(refseq.protein.stats.table)), refseq.protein.stats.table$name)
# setkey(refseq.protein.stats.table, name)

reprseq.table <- fread(reprseq.table.output)
reprseq.table[, accession := tstrsplit(ncbi.id, ":", fixed = TRUE)[[1]]]
ncbi2repr <- setNames(reprseq.table$reprseq, reprseq.table$ncbi.id)
setkey(reprseq.table, reprseq)

genomes.aba.table <- fread(genomes.aba.table.output)
setkey(genomes.aba.table, qname, sname)

reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq2funct <- setNames(reprseq.functions.unique$funct, reprseq.functions.unique$reprseq)
reprseq2funct[reprseq2funct==""] <- "unknown"

# input.reprseq.names <- c("reprseq141902","reprseq098632","reprseq000317",
#                          "reprseq031519","reprseq095134","reprseq091785")
# input.ecf <- "ecf_00381"

input.reprseq.names <- c("reprseq109697","reprseq058101","reprseq094488",
                         "reprseq081403","reprseq032135","reprseq103340")
input.ecf <- "ecf_00076"

res <- run_ecf_synteny_plot(
  ecf_id = input.ecf,
  reprseq_names = input.reprseq.names,
  n_extend = NO.GENES.TO.EXTEND,
  pad = PAD.PARM,
  ecf_subseqs_dir = ecf.subseqs.dir,
  metadata_table = metadata.table,
  refseq_protein_stats = refseq.protein.stats.table,
  name2row = name2row,
  reprseq_table = reprseq.table,
  ncbi2repr = ncbi2repr,
  genomes_aba_table = genomes.aba.table,
  reprseq2funct = reprseq2funct
)

####################################################################
############################## HELPER ##############################
####################################################################


get_gene_neighborhood <- function(reprseq_names,
                                  reprseq_table,
                                  refseq_protein_stats,
                                  name2row,
                                  ncbi2repr,
                                  n_extend = 4L) {
  # Returns gene.table: seq_id, start, end, conf, feat_id
  gene.table.list <- lapply(seq_along(reprseq_names), function(i) {
    repr <- reprseq_names[i]
    
    # reprseq_table keyed by reprseq
    this.link <- reprseq_table[repr]
    this.ncbi.id <- this.link$cluster.name[1]
    
    row_i <- name2row[[this.ncbi.id]]
    row.range <- (row_i - n_extend):(row_i + n_extend)
    row.range <- row.range[row.range >= 1L & row.range <= nrow(refseq_protein_stats)]
    
    dt <- refseq_protein_stats[row.range]
    
    # restrict to the accession/genome of the focal protein
    accession <- strsplit(this.ncbi.id, ":")[[1]][1]
    dt <- dt[seqid %in% accession]
    
    dt[, reprseq := ncbi2repr[name]]
    dt_out <- dt %>% dplyr::select(seq_id = seqid, start, end, conf, feat_id = reprseq)
    data.table::as.data.table(dt_out)
  })
  
  data.table::rbindlist(gene.table.list, use.names = TRUE, fill = TRUE)
}

make_seq_table <- function(all_genomes, metadata_table) {
  seq.table <- data.table::data.table(seq_id = all_genomes)
  seq.metadata <- metadata_table[seq.table, on = .(accession = seq_id)]
  cbind(
    seq.table,
    seq.metadata[, .(length = genome.length, host, genus.ncbi)]
  )
}

make_link_table <- function(all_genomes, genomes_aba_table) {
  input.aba <- genomes_aba_table[qname %chin% all_genomes & sname %chin% all_genomes]
  input.aba <- input.aba[qname < sname]
  input.aba %>%
    dplyr::select(
      seq_id = qname, start = qstart, end = qend,
      seq_id2 = sname, start2 = sstart, end2 = send
    ) %>%
    data.table::as.data.table()
}

make_funct_table <- function(gene_table, reprseq2funct) {
  funct.table <- data.table::data.table(feat_id = sort(unique(gene_table$feat_id)))
  funct.table[, cluster_id := reprseq2funct[feat_id]]
  funct.table
}

load_ecf_coords <- function(ecf_id, reprseq_names, ecf_subseqs_dir) {
  fn <- file.path(ecf_subseqs_dir, sprintf("%s-frag.txt", ecf_id))
  dt <- data.table::fread(fn)
  dt[name %in% reprseq_names]
}


build_ecf_synteny_tables <- function(ecf_id,
                                     reprseq_names,
                                     n_extend = 4L,
                                     ecf_subseqs_dir,
                                     # data objects
                                     metadata_table,
                                     refseq_protein_stats,
                                     name2row,
                                     reprseq_table,
                                     ncbi2repr,
                                     genomes_aba_table,
                                     reprseq2funct) {
  
  # (A) genes around focal protein in each genome
  gene.table <- get_gene_neighborhood(
    reprseq_names = reprseq_names,
    reprseq_table = reprseq_table,
    refseq_protein_stats = refseq_protein_stats,
    name2row = name2row,
    ncbi2repr = ncbi2repr,
    n_extend = n_extend
  )
  
  all.genomes <- unique(gene.table$seq_id)
  
  # (B) genomes table
  seq.table <- make_seq_table(all.genomes, metadata_table)
  
  # (C) links table
  link.table <- make_link_table(all.genomes, genomes_aba_table)
  
  # (D) functions / clusters table
  funct.table <- make_funct_table(gene.table, reprseq2funct)
  
  # (E) ECF feature table (map AA coords -> genome coords)
  ecf.coords <- load_ecf_coords(ecf_id, reprseq_names, ecf_subseqs_dir)
  
  ecf.genome.table <- map_ecf_to_genome(ecf.coords, gene.table,
                                        aa_start_col = "start",
                                        aa_end_col = "end",
                                        ecf_id_col = "name",
                                        feat_id_col = "feat_id",
                                        seq_id_col = "seq_id",
                                        gene_start_col = "start",
                                        gene_end_col = "end",
                                        include_stop = FALSE)
  
  ecf.feature.table <- ecf.genome.table %>%
    dplyr::select(seq_id, start = ecf_genome_start, end = ecf_genome_end) %>%
    data.table::as.data.table()
  
  list(
    ecf_id = ecf_id,
    reprseq_names = reprseq_names,
    gene = gene.table,
    seq = seq.table,
    link = link.table,
    funct = funct.table,
    feat = ecf.feature.table
  )
}

zoom_ecf_tables <- function(tables, pad = 10L, clip_links = TRUE) {
  stopifnot(all(c("gene","seq","link","feat") %in% names(tables)))
  
  win <- define_genome_windows(tables$gene, tables$seq, pad = pad)
  
  list(
    win = win,
    seq  = make_zoomed_seqs(tables$seq, win),
    gene = rebase_genes_to_windows(tables$gene, win),
    link = rebase_links_to_windows(tables$link, win, clip_links = clip_links, require_both = TRUE),
    feat = rebase_features_to_windows(tables$feat, win),
    funct = tables$funct,          # unchanged
    ecf_id = tables$ecf_id,
    reprseq_names = tables$reprseq_names
  )
}

plot_ecf_gggenomes <- function(tables,
                               feat_col = "red",
                               feat_size = 8,
                               fill_palette = "Spectral") {
  
  # Minimal required set
  stopifnot(all(c("seq","gene","link","feat","funct") %in% names(tables)))
  
  gggenomes(
    seqs  = tibble::as_tibble(tables$seq),
    genes = tibble::as_tibble(tables$gene),
    links = tibble::as_tibble(tables$link),
    feats = tibble::as_tibble(tables$feat)
  ) %>%
    add_clusters(tibble::as_tibble(tables$funct)) +
    geom_seq() + geom_bin_label() +
    geom_gene(aes(fill = cluster_id)) +
    geom_link() +
    geom_feat(size = feat_size, col = feat_col) +
    scale_fill_brewer(palette = fill_palette)
}

run_ecf_synteny_plot <- function(ecf_id,
                                 reprseq_names,
                                 n_extend = 4L,
                                 pad = 10L,
                                 ecf_subseqs_dir,
                                 metadata_table,
                                 refseq_protein_stats,
                                 name2row,
                                 reprseq_table,
                                 ncbi2repr,
                                 genomes_aba_table,
                                 reprseq2funct) {
  
  base <- build_ecf_synteny_tables(
    ecf_id = ecf_id,
    reprseq_names = reprseq_names,
    n_extend = n_extend,
    ecf_subseqs_dir = ecf_subseqs_dir,
    metadata_table = metadata_table,
    refseq_protein_stats = refseq_protein_stats,
    name2row = name2row,
    reprseq_table = reprseq_table,
    ncbi2repr = ncbi2repr,
    genomes_aba_table = genomes_aba_table,
    reprseq2funct = reprseq2funct
  )
  
  zoom <- zoom_ecf_tables(base, pad = pad, clip_links = TRUE)
  
  return.plot <- plot_ecf_gggenomes(zoom, fill_palette = "Set1")
  return(return.plot)
  
  # list(
  #   base_tables = base,
  #   zoom_tables = zoom,
  #   base_plot = plot_ecf_gggenomes(base, fill_palette = "Spectral"),
  #   zoom_plot = plot_ecf_gggenomes(zoom, fill_palette = "Set1")
  # )
}






# # Define parms
# NO.GENES.TO.EXTEND <- 4
# 
# ecf.subseqs.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs", "ecf-subseqs")
# 
# # Load data
# metadata.table <- fread(metadata.table.output)
# refseq.protein.stats.table <- fread(refseq.protein.stats.table.filename)
# name2row <- setNames(seq_len(nrow(refseq.protein.stats.table)), refseq.protein.stats.table$name)
# # setkey(refseq.protein.stats.table, name)
# 
# reprseq.table <- fread(reprseq.table.output)
# reprseq.table[, accession := tstrsplit(ncbi.id, ":", fixed = TRUE)[[1]]]
# ncbi2repr <- setNames(reprseq.table$reprseq, reprseq.table$ncbi.id)
# setkey(reprseq.table, reprseq)
# 
# genomes.aba.table <- fread(genomes.aba.table.output)
# setkey(genomes.aba.table, qname, sname)
# 
# reprseq.functions.unique <- fread(reprseq.functions.unique.output)
# reprseq2funct <- setNames(reprseq.functions.unique$funct, reprseq.functions.unique$reprseq)
# reprseq2funct[reprseq2funct==""] <- "unknown"
# 
# # Define input pamrs
# input.reprseq.names <- c("reprseq141902", "reprseq098632", "reprseq000317",
#                          "reprseq031519", "reprseq095134", "reprseq091785")
# input.ecf <- "ecf_00381"
# ecf.coords.table.fn <- file.path(ecf.subseqs.dir, sprintf("%s-frag.txt", input.ecf))
# ecf.coords.table <- fread(ecf.coords.table.fn)
# ecf.coords.table <- ecf.coords.table[name %in% input.reprseq.names]
# 
# ecf.genome.table <- map_ecf_to_genome(ecf.coords.table, gene.table)
# ecf.feature.table <- ecf.genome.table %>% select(seq_id, start = ecf_genome_start, end = ecf_genome_end)
# 
# 
# # 1. Make a table of gene coordinates for the focal region
# gene.table.list <- lapply(1:length(input.reprseq.names), function(index){
#   this.input.seq <- input.reprseq.names[index]
#   this.input.seq.genome.link <- reprseq.table[this.input.seq]
#   this.input.seq.ncbi.id <- this.input.seq.genome.link$cluster.name[1]
#   this.input.seq.accession <- strsplit(this.input.seq.ncbi.id, ":")[[1]][1]
#   # this.input.seq.ncbi.id.row.index <- refseq.protein.stats.table[this.input.seq.ncbi.id, which = T]
#   this.input.seq.ncbi.id.row.index <- name2row[[this.input.seq.ncbi.id]]
#   row.range <- (this.input.seq.ncbi.id.row.index-NO.GENES.TO.EXTEND):(this.input.seq.ncbi.id.row.index+NO.GENES.TO.EXTEND)
#   row.range <- row.range[row.range >= 1L & row.range <= nrow(refseq.protein.stats.table)]
#   this.input.seq.gene.table <- refseq.protein.stats.table[row.range]
#   this.input.seq.gene.table <- this.input.seq.gene.table[seqid %in% this.input.seq.accession]
#   this.input.seq.gene.table[, reprseq := ncbi2repr[name]]
#   this.input.seq.gene.table.out <- this.input.seq.gene.table %>% select(seq_id = seqid, start, end, conf, feat_id = reprseq)
#   return(this.input.seq.gene.table.out)
# })
# gene.table <- rbindlist(gene.table.list)
# 
# all.genomes <- unique(gene.table$seq_id)
# 
# # 2. Make a table of genomes
# seq.table <- data.table(seq_id = all.genomes)
# seq.metadata <- metadata.table[seq.table, on = .(accession = seq_id)]
# seq.table <- cbind(seq.table,
#                    seq.metadata[, .(length = genome.length, host, genus.ncbi)])
# 
# # 3. Make a table of links
# input.aba.table <- genomes.aba.table[qname %chin% all.genomes & sname %chin% all.genomes]
# input.aba.table <- input.aba.table[qname < sname]
# link.table <- input.aba.table %>% select(seq_id = qname, start = qstart, end = qend,
#                                          seq_id2 = sname, start2 = sstart, end2 = send)
# 
# # 4. Make functions table
# funct.table <- data.table(feat_id = sort(unique(gene.table$feat_id)))
# funct.table[, cluster_id := reprseq2funct[feat_id]]
# 
# # TRY PLOTTING 
# 
# # Prepare tracks
# ggg.seq.track <- tibble::tibble(seq.table)
# ggg.gene.track <- tibble::tibble(gene.table)
# ggg.link.track <- tibble::tibble(link.table)
# ggg.funct.table <- tibble::tibble(funct.table)
# ggg.feat.table <- tibble::tibble(ecf.feature.table)
# 
# gggenomes(seqs = ggg.seq.track, 
#           genes = ggg.gene.track, 
#           links = ggg.link.track,
#           feats = ggg.feat.table) %>%
#   add_clusters(ggg.funct.table) + 
#   geom_seq() + geom_bin_label() +
#   geom_gene(aes(fill = cluster_id)) +
#   geom_link() + geom_feat(size = 8, col = "red") +
#   scale_fill_brewer(palette="Spectral")
# 
# 
# # ZOOM IN ON THE SELECTED REGION
# 
# PAD.PARM=10L
# 
# win <- define_genome_windows(gene.table, seq.table, pad = PAD.PARM)
# genes_zoom <- rebase_genes_to_windows(gene.table, win)
# links_zoom <- rebase_links_to_windows(link.table, win, clip_links = TRUE)
# seqs_zoom <- make_zoomed_seqs(seq.table, win)
# feats_zoom <- rebase_features_to_windows(ecf.feature.table, win)
# 
# 
# gggenomes(
#   seqs  = as_tibble(seqs_zoom),
#   genes = as_tibble(genes_zoom),
#   links = as_tibble(links_zoom),
#   feats = as_tibble(feats_zoom)
# ) %>%
#   add_clusters(ggg.funct.table) +
#   geom_seq() + geom_bin_label() +
#   geom_feat(size=8, col = "red") +
#   geom_gene(aes(fill = cluster_id)) +
#   geom_link() + 
#   scale_fill_brewer(palette="Set1")
# 
# 
# 
# ######################################################################################################
# ######################################################################################################
# ######################################################################################################
# ######################################################################################################
# 
map_ecf_to_genome <- function(ecf_dt, gene_dt,
                              aa_start_col = "start",
                              aa_end_col   = "end",
                              ecf_id_col   = "name",
                              feat_id_col  = "feat_id",
                              seq_id_col   = "seq_id",
                              gene_start_col = "start",
                              gene_end_col   = "end",
                              include_stop = FALSE) {

  ecf_dt  <- copy(as.data.table(ecf_dt))
  gene_dt <- copy(as.data.table(gene_dt))

  # Basic checks
  stopifnot(all(c(ecf_id_col, aa_start_col, aa_end_col) %in% names(ecf_dt)))
  stopifnot(all(c(seq_id_col, feat_id_col, gene_start_col, gene_end_col) %in% names(gene_dt)))

  # Rename to standard internal names
  ecf <- ecf_dt[, .(
    feat_id  = get(ecf_id_col),
    aa_start = as.integer(get(aa_start_col)),
    aa_end   = as.integer(get(aa_end_col))
  )]

  genes <- gene_dt[, .(
    seq_id = get(seq_id_col),
    feat_id = get(feat_id_col),
    g_start = as.integer(get(gene_start_col)),
    g_end   = as.integer(get(gene_end_col))
  )]

  # Join ECF definitions onto each genome occurrence of that RepSeq
  out <- merge(genes, ecf, by = "feat_id", allow.cartesian = TRUE)

  # Strand inference
  out[, strand := fifelse(g_start <= g_end, "+", "-")]

  # Gene bounds and length (nt)
  out[, g_min := pmin(g_start, g_end)]
  out[, g_max := pmax(g_start, g_end)]
  out[, gene_len_nt := g_max - g_min + 1L]

  # Convert AA positions to nucleotide offsets within CDS (1-based nt within CDS)
  # AA i spans nt: (i-1)*3+1 ... (i-1)*3+3
  out[, nt_start_in_cds := (aa_start - 1L) * 3L + 1L]
  out[, nt_end_in_cds   := (aa_end   - 1L) * 3L + 3L]

  # Optionally include stop codon after aa_end (rarely desired for domains)
  if (include_stop) {
    out[, nt_end_in_cds := nt_end_in_cds + 3L]
  }

  # Sanity checks: warn/clip if ECF exceeds CDS length
  out[, exceeds := nt_end_in_cds > gene_len_nt]
  if (any(out$exceeds, na.rm = TRUE)) {
    warning("Some ECF coordinates extend beyond CDS length in at least one genome occurrence. See column `exceeds`.")
  }

  # Map into genome coordinates depending on strand
  # Forward: genome_nt = g_min + offset - 1
  out[strand == "+", `:=`(
    ecf_start_nt = g_min + nt_start_in_cds - 1L,
    ecf_end_nt   = g_min + nt_end_in_cds   - 1L
  )]

  # Reverse: CDS runs from g_max downwards
  # nt_start_in_cds counted from 5' end of CDS (which is at g_max on genome if reverse)
  # so genome coordinate decreases with offset
  out[strand == "-", `:=`(
    ecf_start_nt = g_max - (nt_start_in_cds - 1L),
    ecf_end_nt   = g_max - (nt_end_in_cds   - 1L)
  )]

  # For consistency, store genomic interval as start<=end plus strand
  out[, ecf_min_nt := pmin(ecf_start_nt, ecf_end_nt)]
  out[, ecf_max_nt := pmax(ecf_start_nt, ecf_end_nt)]

  # Return a clean table
  out[
    , .(
      seq_id,
      feat_id,
      strand,
      gene_start = g_start,
      gene_end   = g_end,
      aa_start,
      aa_end,
      ecf_genome_start = ecf_min_nt,
      ecf_genome_end   = ecf_max_nt,
      gene_len_nt,
      nt_start_in_cds,
      nt_end_in_cds,
      exceeds
    )
  ][order(seq_id, feat_id, ecf_genome_start)]
}

define_genome_windows <- function(gene_dt, seq_dt, pad = 0L) {
  gene_dt <- copy(as.data.table(gene_dt))
  seq_dt  <- copy(as.data.table(seq_dt))

  stopifnot(all(c("seq_id", "start", "end") %in% names(gene_dt)))
  stopifnot(all(c("seq_id", "length") %in% names(seq_dt)))

  win <- gene_dt[
    , .(
      win_start = as.integer(min(pmin(start, end)) - pad),
      win_end   = as.integer(max(pmax(start, end)) + pad)
    ),
    by = seq_id
  ]

  # add genome length and clamp window to [1, length]
  win <- merge(
    win,
    seq_dt[, .(seq_id, genome_length = as.integer(length))],
    by = "seq_id",
    all.x = TRUE
  )

  win[is.na(genome_length), genome_length := .Machine$integer.max]
  win[win_start < 1L, win_start := 1L]
  win[win_end   > genome_length, win_end := genome_length]

  win[, offset := win_start - 1L]
  win[, zoom_length := as.integer(win_end - win_start + 1L)]

  setkey(win, seq_id)
  win[]
}

rebase_genes_to_windows <- function(gene_dt, win_dt, keep_only_overlapping = TRUE) {
  gene_dt <- copy(as.data.table(gene_dt))
  win_dt  <- copy(as.data.table(win_dt))

  stopifnot(all(c("seq_id", "start", "end") %in% names(gene_dt)))
  stopifnot(all(c("seq_id", "win_start", "win_end", "offset") %in% names(win_dt)))

  setkey(win_dt, seq_id)

  genes <- merge(gene_dt, win_dt[, .(seq_id, win_start, win_end, offset)], by = "seq_id", all.x = FALSE)

  if (keep_only_overlapping) {
    genes <- genes[
      pmax(start, end) >= win_start &
        pmin(start, end) <= win_end
    ]
  }

  genes[, `:=`(
    start = start - offset,
    end   = end   - offset
  )]

  genes[]
}

rebase_features_to_windows <- function(feat_dt, win_dt, keep_only_overlapping = TRUE) {
  feat_dt <- copy(as.data.table(feat_dt))
  win_dt  <- copy(as.data.table(win_dt))

  stopifnot(all(c("seq_id", "start", "end") %in% names(feat_dt)))
  stopifnot(all(c("seq_id", "win_start", "win_end", "offset") %in% names(win_dt)))

  setkey(win_dt, seq_id)

  feats <- merge(feat_dt, win_dt[, .(seq_id, win_start, win_end, offset)], by = "seq_id", all.x = FALSE)

  if (keep_only_overlapping) {
    feats <- feats[
      pmax(start, end) >= win_start &
        pmin(start, end) <= win_end
    ]
  }

  feats[, `:=`(
    start = start - offset,
    end   = end   - offset
  )]

  feats[]
}

rebase_links_to_windows <- function(link_dt, win_dt, clip_links = TRUE, require_both = TRUE) {
  link_dt <- copy(as.data.table(link_dt))
  win_dt  <- copy(as.data.table(win_dt))

  stopifnot(all(c("seq_id", "start", "end", "seq_id2", "start2", "end2") %in% names(link_dt)))
  stopifnot(all(c("seq_id", "win_start", "win_end", "offset") %in% names(win_dt)))

  w1 <- win_dt[, .(seq_id, win_start1 = win_start, win_end1 = win_end, offset1 = offset)]
  w2 <- win_dt[, .(seq_id2 = seq_id, win_start2 = win_start, win_end2 = win_end, offset2 = offset)]

  links <- merge(link_dt, w1, by = "seq_id",  all.x = TRUE)
  links <- merge(links,  w2, by = "seq_id2", all.x = TRUE)

  # Keep only links that have windows defined
  if (require_both) {
    links <- links[!is.na(offset1) & !is.na(offset2)]
  } else {
    links <- links[!is.na(offset1) | !is.na(offset2)]
  }

  if (clip_links) {
    # overlap tests
    ok1 <- links[, end  >= win_start1 & start  <= win_end1]
    ok2 <- links[, end2 >= win_start2 & start2 <= win_end2]

    links <- if (require_both) links[ok1 & ok2] else links[ok1 | ok2]

    # clip (only where window exists)
    links[!is.na(win_start1), `:=`(
      start = pmax(start, win_start1),
      end   = pmin(end,   win_end1)
    )]
    links[!is.na(win_start2), `:=`(
      start2 = pmax(start2, win_start2),
      end2   = pmin(end2,   win_end2)
    )]
  }

  # rebase (only where window exists)
  links[!is.na(offset1), `:=`(
    start = start - offset1,
    end   = end   - offset1
  )]
  links[!is.na(offset2), `:=`(
    start2 = start2 - offset2,
    end2   = end2   - offset2
  )]

  links[]
}

make_zoomed_seqs <- function(seq_dt, win_dt, keep_only_windowed = TRUE) {
  seq_dt <- copy(as.data.table(seq_dt))
  win_dt <- copy(as.data.table(win_dt))

  stopifnot(all(c("seq_id", "length") %in% names(seq_dt)))
  stopifnot(all(c("seq_id", "win_start", "win_end") %in% names(win_dt)))

  seqs <- merge(seq_dt, win_dt[, .(seq_id, win_start, win_end)], by = "seq_id", all.x = TRUE)

  if (keep_only_windowed) seqs <- seqs[!is.na(win_start)]

  seqs[, orig_length := as.integer(length)]
  seqs[, length := as.integer(win_end - win_start + 1L)]

  seqs[]
}


