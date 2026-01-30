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


