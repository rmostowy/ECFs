### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(gggenomes)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

source("plot-ecf-genomes-functs.R")

# Define parms
NO.GENES.TO.EXTEND <- 4

output.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs")
ecf.subseqs.dir <- file.path(output.dir, "ecf-subseqs")
ecf.clusters.dir <- file.path(output.dir, "ecf-clusters")

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

reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
reprseq.seqs.n <- names(reprseq.seqs)

# Define input pamrs
# input.reprseq.names <- c("reprseq141902", "reprseq098632", "reprseq000317",
#                          "reprseq031519", "reprseq095134", "reprseq091785")
# input.ecf <- "ecf_00381"

input.reprseq.names <- c("reprseq109697","reprseq058101","reprseq094488",
                         "reprseq081403","reprseq037135","reprseq103340")
input.ecf <- "ecf_00076"
input.ecf.dir <- file.path(ecf.clusters.dir, input.ecf)
make.dir(input.ecf.dir)

input.reprseq.seqs <- reprseq.seqs[reprseq.seqs.n %in% input.reprseq.names]
input.reprseq.seqs.outfile <- file.path(input.ecf.dir, "custom-reprseqs.fa")
writeXStringSet(input.reprseq.seqs, filepath = input.reprseq.seqs.outfile)


ecf.coords.table.fn <- file.path(ecf.subseqs.dir, sprintf("%s-frag.txt", input.ecf))
ecf.coords.table <- fread(ecf.coords.table.fn)
ecf.coords.table <- ecf.coords.table[name %in% input.reprseq.names]


# 1. Make a table of gene coordinates for the focal region
gene.table.list <- lapply(1:length(input.reprseq.names), function(index){
  this.input.seq <- input.reprseq.names[index]
  this.input.seq.genome.link <- reprseq.table[this.input.seq]
  this.input.seq.ncbi.id <- this.input.seq.genome.link$cluster.name[1]
  this.input.seq.accession <- strsplit(this.input.seq.ncbi.id, ":")[[1]][1]
  # this.input.seq.ncbi.id.row.index <- refseq.protein.stats.table[this.input.seq.ncbi.id, which = T]
  this.input.seq.ncbi.id.row.index <- name2row[[this.input.seq.ncbi.id]]
  row.range <- (this.input.seq.ncbi.id.row.index-NO.GENES.TO.EXTEND):(this.input.seq.ncbi.id.row.index+NO.GENES.TO.EXTEND)
  row.range <- row.range[row.range >= 1L & row.range <= nrow(refseq.protein.stats.table)]
  this.input.seq.gene.table <- refseq.protein.stats.table[row.range]
  this.input.seq.gene.table <- this.input.seq.gene.table[seqid %in% this.input.seq.accession]
  this.input.seq.gene.table[, reprseq := ncbi2repr[name]]
  this.input.seq.gene.table.out <- this.input.seq.gene.table %>% select(seq_id = seqid, start, end, conf, feat_id = reprseq)
  return(this.input.seq.gene.table.out)
})
gene.table <- rbindlist(gene.table.list)
all.genomes <- unique(gene.table$seq_id)

# 2. Map ECF to genomes
ecf.genome.table <- map_ecf_to_genome(ecf.coords.table, gene.table)
ecf.feature.table <- ecf.genome.table %>% select(seq_id, start = ecf_genome_start, end = ecf_genome_end)

# 2. Make a table of genomes
seq.table <- data.table(seq_id = all.genomes)
seq.metadata <- metadata.table[seq.table, on = .(accession = seq_id)]
seq.table <- cbind(seq.table,
                   seq.metadata[, .(length = genome.length, host, genus.ncbi)])

# 3. Make a table of links
input.aba.table <- genomes.aba.table[qname %chin% all.genomes & sname %chin% all.genomes]
input.aba.table <- input.aba.table[qname < sname]
link.table <- input.aba.table %>% select(seq_id = qname, start = qstart, end = qend,
                                         seq_id2 = sname, start2 = sstart, end2 = send)

# 4. Make functions table
funct.table <- data.table(feat_id = sort(unique(gene.table$feat_id)))
funct.table[, cluster_id := reprseq2funct[feat_id]]

# TRY PLOTTING 

# Prepare tracks
ggg.seq.track <- tibble::tibble(seq.table)
ggg.gene.track <- tibble::tibble(gene.table)
ggg.link.track <- tibble::tibble(link.table)
ggg.funct.table <- tibble::tibble(funct.table)
ggg.feat.table <- tibble::tibble(ecf.feature.table)

gggenomes(seqs = ggg.seq.track, 
          genes = ggg.gene.track, 
          links = ggg.link.track,
          feats = ggg.feat.table) %>%
  add_clusters(ggg.funct.table) + 
  geom_seq() + geom_bin_label() +
  geom_gene(aes(fill = cluster_id)) +
  geom_link() + geom_feat(size = 8, col = "red") +
  scale_fill_brewer(palette="Spectral")


# ZOOM IN ON THE SELECTED REGION

PAD.PARM=10L

win <- define_genome_windows(gene.table, seq.table, pad = PAD.PARM)
genes_zoom <- rebase_genes_to_windows(gene.table, win)
links_zoom <- rebase_links_to_windows(link.table, win, clip_links = TRUE)
seqs_zoom <- make_zoomed_seqs(seq.table, win)
feats_zoom <- rebase_features_to_windows(ecf.feature.table, win)


gggenomes(
  seqs  = as_tibble(seqs_zoom),
  genes = as_tibble(genes_zoom),
  links = as_tibble(links_zoom),
  feats = as_tibble(feats_zoom)
) %>%
  add_clusters(ggg.funct.table) +
  geom_seq() + geom_bin_label() +
  geom_feat(size=8, col = "red") +
  geom_gene(aes(fill = cluster_id)) +
  geom_link() + 
  scale_fill_brewer(palette="Set1")

p <- gggenomes(
  seqs  = as_tibble(seqs_zoom),
  genes = as_tibble(genes_zoom),
  feats = as_tibble(feats_zoom)
) %>%
  add_clusters(ggg.funct.table) +
  geom_seq() + geom_bin_label() +
  geom_feat(size=8, col = "red") +
  geom_gene(aes(fill = cluster_id)) +
  scale_fill_brewer(palette="Set1")
p.outfile <- file.path(input.ecf.dir, "custom-plot.pdf")
ggsave(p.outfile, p, width = 8, height = 5)
