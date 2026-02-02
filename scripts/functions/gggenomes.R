plot.ecf.synteny <- function(x,
                             palette = "Set1",
                             feat_col = "red",
                             feat_size = 8) {
  
  stopifnot(is.list(x))
  stopifnot(all(c("seq.table", "gene.table", "funct.table", "feat.table") %in% names(x)))
  
  seqs  <- tibble::as_tibble(x$seq.table)
  genes <- tibble::as_tibble(x$gene.table)
  feats <- tibble::as_tibble(x$feat.table)
  funct <- tibble::as_tibble(x$funct.table)
  
  has_links <- ("link.table" %in% names(x)) &&
    !is.null(x$link.table) &&
    nrow(x$link.table) > 0
  
  if (has_links) {
    links <- tibble::as_tibble(x$link.table)
    
    p <- gggenomes::gggenomes(
      seqs  = seqs,
      genes = genes,
      links = links,
      feats = feats
    ) %>%
      gggenomes::add_clusters(funct) +
      gggenomes::geom_seq() + gggenomes::geom_bin_label() +
      gggenomes::geom_feat(size = feat_size, col = feat_col) +
      gggenomes::geom_gene(ggplot2::aes(fill = cluster_id)) +
      gggenomes::geom_link() +
      ggplot2::scale_fill_brewer(palette = palette)
    
  } else {
    p <- gggenomes::gggenomes(
      seqs  = seqs,
      genes = genes,
      feats = feats
    ) %>%
      gggenomes::add_clusters(funct) +
      gggenomes::geom_seq() + gggenomes::geom_bin_label() +
      gggenomes::geom_feat(size = feat_size, col = feat_col) +
      gggenomes::geom_gene(ggplot2::aes(fill = cluster_id)) +
      ggplot2::scale_fill_brewer(palette = palette)
  }
  
  p
}

prepare.gggenomes.tables.zoomed <- function(input.list = gggenomes.inputs, pad.size = 0L){
  window.table <- define_genome_windows(input.list$gene.table, input.list$seq.table, pad.size)
  gene.table.zoom <- rebase_genes_to_windows(input.list$gene.table, window.table)
  link.table.zoom <- rebase_links_to_windows(input.list$link.table, window.table, clip_links = TRUE)
  seq.table.zoom <- make_zoomed_seqs(input.list$seq.table, window.table)
  feat.table.zoom <- rebase_features_to_windows(input.list$feat.table, window.table)  
  
  list(seq.table = seq.table.zoom,
       gene.table = gene.table.zoom, 
       link.table = link.table.zoom,
       funct.table = input.list$funct.table,
       feat.table = feat.table.zoom)  
  
}

# Redefine coordinates to zoom in on the selected region



prepare.gggenomes.tables <- function(reprseq.names,
                                     ecf.coords.dt = ecf.coords.table,
                                     metadata = metadata.table,
                                     prodigal.table = refseq.protein.stats.table,
                                     reprseq.dt = reprseq.table,
                                     genomes.aba.dt = genomes.aba.table,
                                     reprseq2functs = reprseq.functions.unique,
                                     no.genes.upstream.downstream = 2L){
  
  reprseq.dt[, accession := tstrsplit(ncbi.id, ":", fixed = TRUE)[[1]]]
  reprseq.table.ref <- reprseq.dt[cluster.name == ncbi.id]
  repr2accession <- setNames(reprseq.table.ref$accession, reprseq.table.ref$reprseq)
  input.reprseq.genomes <- repr2accession[reprseq.names]
  name2row <- setNames(seq_len(nrow(prodigal.table)), prodigal.table$name)
  ncbi2repr <- setNames(reprseq.dt$reprseq, reprseq.dt$ncbi.id)
  setkey(reprseq.dt, reprseq)
  setkey(genomes.aba.dt, qname, sname)
  reprseq2funct <- setNames(reprseq2functs$funct, reprseq2functs$reprseq)
  reprseq2funct[reprseq2funct==""] <- "unknown"
  
  
  # 1. Make a table of gene coordinates for the focal region
  gene.table.list <- lapply(1:length(reprseq.names), function(index){
    this.input.seq <- reprseq.names[index]
    this.input.seq.genome.link <- reprseq.dt[this.input.seq]
    this.input.seq.ncbi.id <- this.input.seq.genome.link$cluster.name[1]
    this.input.seq.accession <- strsplit(this.input.seq.ncbi.id, ":")[[1]][1]
    this.input.seq.ncbi.id.row.index <- name2row[[this.input.seq.ncbi.id]]
    row.range <- (this.input.seq.ncbi.id.row.index-no.genes.upstream.downstream):(this.input.seq.ncbi.id.row.index+no.genes.upstream.downstream)
    row.range <- row.range[row.range >= 1L & row.range <= nrow(prodigal.table)]
    this.input.seq.gene.table <- prodigal.table[row.range]
    this.input.seq.gene.table <- this.input.seq.gene.table[seqid %in% this.input.seq.accession]
    this.input.seq.gene.table[, reprseq := ncbi2repr[name]]
    this.input.seq.gene.table.out <- this.input.seq.gene.table %>% select(seq_id = seqid, start, end, conf, feat_id = reprseq, strand)
    return(this.input.seq.gene.table.out)
  })
  gene.table <- rbindlist(gene.table.list)
  # all.genomes <- unique(gene.table$seq_id)
  
  # 2. Map ECF to genomes
  
  ecf.feature.table.list <- lapply(1:length(reprseq.names), function(index){
    this.reprseq <- reprseq.names[index]
    this.reprseq.accession <- input.reprseq.genomes[index]
    this.reprseq.ecf.table <- ecf.coords.table[name == this.reprseq]
    this.reprseq.ecf.start <- this.reprseq.ecf.table$start
    this.reprseq.ecf.end <- this.reprseq.ecf.table$end
    this.reprseq.length <- this.reprseq.ecf.table$length
    this.reprseq.coords <- gene.table[seq_id == this.reprseq.accession & feat_id == this.reprseq]
    this.reprseq.start <- this.reprseq.coords$start
    this.reprseq.end <- this.reprseq.coords$end
    this.reprseq.strand <- this.reprseq.coords$strand
    this.ecf.coords <- map.ecf.coords.to.genome(this.reprseq.ecf.start, this.reprseq.ecf.end, this.reprseq.start, this.reprseq.end, this.reprseq.strand, this.reprseq.length)
    data.table(seq_id = this.reprseq.accession, feat_id = this.reprseq, start = this.ecf.coords["start"], end = this.ecf.coords["end"])
  })
  ecf.feature.table <- rbindlist(ecf.feature.table.list)
  
  # 2. Make a table of genomes
  seq.table <- data.table(seq_id = input.reprseq.genomes)
  seq.metadata <- metadata[seq.table, on = .(accession = seq_id)]
  seq.table <- cbind(seq.table,
                     seq.metadata[, .(length = genome.length, host, genus.ncbi)])
  
  # 3. Make a table of links
  input.aba.table <- genomes.aba.dt[qname %chin% input.reprseq.genomes & sname %chin% input.reprseq.genomes]
  input.aba.table <- input.aba.table[qname < sname]
  link.table <- input.aba.table %>% select(seq_id = qname, start = qstart, end = qend,
                                           seq_id2 = sname, start2 = sstart, end2 = send)
  
  # 4. Make functions table
  funct.table <- data.table(feat_id = sort(unique(gene.table$feat_id)))
  funct.table[, cluster_id := reprseq2funct[feat_id]]
  
  list(seq.table = seq.table,
       gene.table = gene.table, link.table = link.table,
       funct.table = funct.table,
       feat.table = ecf.feature.table)
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


