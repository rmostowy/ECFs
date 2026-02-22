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

# source("plot-ecf-genomes-functs.R")

#### Define parms
NO.GENES.TO.EXTEND <- 4L
PAD.PARM=10L

#### Define paths
output.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs")
ecf.subseqs.dir <- file.path(output.dir, "ecf-subseqs")
ecf.clusters.dir <- file.path(output.dir, "ecf-clusters")

#### Load data
metadata.table <- fread(metadata.table.output)
refseq.protein.stats.table <- fread(refseq.protein.stats.table.filename)
reprseq.table <- fread(reprseq.table.output)
genomes.aba.table <- fread(genomes.aba.table.output)
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
reprseq.seqs.n <- names(reprseq.seqs)

#### Define inpupt parms
# input.reprseq.names <- c("reprseq141902", "reprseq098632", "reprseq000317", "reprseq031519", "reprseq095134", "reprseq091785")
# input.ecf <- "ecf_00381"

# input.reprseq.names <- c("reprseq109697","reprseq058101","reprseq094488", "reprseq081403","reprseq037135","reprseq103340")
# input.ecf <- "ecf_00076"

# input.reprseq.names <- c("reprseq097811","reprseq016828", "reprseq041485", "reprseq052714", "reprseq138546", "reprseq088582", "reprseq071489", "reprseq090351")
# input.ecf <- "ecf_00470"

input.reprseq.names <- c("reprseq146207", "reprseq120314", "reprseq136748", "reprseq027134", "reprseq130930", "reprseq007270")
input.ecf <- "ecf_01247"

input.reprseq.names <- c("reprseq133525", "reprseq000993", "reprseq096185", "reprseq069508", "reprseq010606")
input.ecf <- "ecf_00014"

# input.reprseq.names <- c("reprseq", "reprseq", "reprseq", "reprseq", "reprseq", "reprseq")
# input.ecf <- "ecf_0"


#### Create output directory
input.ecf.dir <- file.path(ecf.clusters.dir, input.ecf)
make.dir(input.ecf.dir)


# Extract sequences of reprseqs
input.reprseq.seqs <- reprseq.seqs[reprseq.seqs.n %in% input.reprseq.names]
input.reprseq.seqs.outfile <- file.path(input.ecf.dir, sprintf("custom-%s.fa", input.ecf))
writeXStringSet(input.reprseq.seqs, filepath = input.reprseq.seqs.outfile)


ecf.coords.table.fn <- file.path(ecf.subseqs.dir, sprintf("%s-frag.txt", input.ecf))
ecf.coords.table <- fread(ecf.coords.table.fn)

# Prepare gggenomes inputs
gggenomes.inputs <- prepare.gggenomes.tables(input.reprseq.names,
                                             ecf.coords.table,
                                             metadata.table,
                                             refseq.protein.stats.table,
                                             reprseq.table,
                                             genomes.aba.table,
                                             reprseq.functions.unique,
                                             no.genes.upstream.downstream = NO.GENES.TO.EXTEND)


gggenomes.inputs.zoom <- prepare.gggenomes.tables.zoomed(gggenomes.inputs, pad = PAD.PARM)

p <- plot.ecf.synteny(gggenomes.inputs.zoom, palette = "Set3")
p.outfile <- file.path(input.ecf.dir, sprintf("custom-%s.pdf", input.ecf))
ggsave(p.outfile, p, width = 8, height = 5)
