### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(gggenomes)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

source("functions.R")

# Define output folder
output.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs")
ecf.clusters.dir <- file.path(output.dir, "ecf-clusters")
ecf.subseqs.dir <- file.path(output.dir, "ecf-subseqs")
make.dir(ecf.clusters.dir)

#### Load data
metadata.table <- fread(metadata.table.output)
refseq.protein.stats.table <- fread(refseq.protein.stats.table.filename)
reprseq.table <- fread(reprseq.table.output)
reprseq.table.ref <- reprseq.table[cluster.name == ncbi.id]
reprseq.table.ref[, accession := tstrsplit(ncbi.id, ":", fixed = TRUE)[[1]]]
reprseq2accession <- setNames(reprseq.table.ref$accession, reprseq.table.ref$reprseq)

genomes.aba.table <- fread(genomes.aba.table.output)
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
reprseq.seqs.n <- names(reprseq.seqs)

all.ecfs <- list.files(ecf.subseqs.dir)
all.ecfs <- unique(substr(all.ecfs, 1, 9))

for(index in 1:length(all.ecfs)){
# for(index in 1:10){
  this.ecf.name <- all.ecfs[index]
  this.ecf.dir <- file.path(ecf.clusters.dir, this.ecf.name)
  this.ecf.aln.filename <- file.path(ecf.subseqs.dir, sprintf("%s.aln", this.ecf.name))
  this.ecf.cl.filename <- file.path(ecf.subseqs.dir, sprintf("%s-cl.txt", this.ecf.name))
  
  cat(sprintf("Examining %s:\n", this.ecf.name))
  this.ecf.aln <- readAAStringSet(this.ecf.aln.filename)
  
  this.ecf.coords.table.fn <- file.path(ecf.subseqs.dir, sprintf("%s-frag.txt", this.ecf.name))
  this.ecf.coords.table <- fread(this.ecf.coords.table.fn)
  
  this.ecf.cl <- fread(this.ecf.cl.filename)
  this.ecf.cl.sizes <- table(this.ecf.cl$cluster)
  no.suffiently.large.clusters <- lw(this.ecf.cl.sizes >= RE.MIN.CLUSTER.SIZE)
  suffiently.large.clusters.exist <- no.suffiently.large.clusters>0
  # suffiently.large.clusters.exist
  
  if(suffiently.large.clusters.exist){
    cat(sprintf("Found %d large clusters! Making plots...\n", no.suffiently.large.clusters))
    make.dir(this.ecf.dir)
    cl.indices <- which(this.ecf.cl.sizes >= RE.MIN.CLUSTER.SIZE)
    
    for(ci in 1:no.suffiently.large.clusters){
      this.cluster.index <- cl.indices[ci]
      this.cluster.members <- this.ecf.cl[cluster == this.cluster.index]$reprseq
      
      # if members come from identical genomes, remove them to facilitate plotting
      this.cluster.members.accessions <- reprseq2accession[this.cluster.members]
      duplicated.accessions <- any(duplicated(this.cluster.members.accessions))
      if(duplicated.accessions) this.cluster.members <- this.cluster.members[!duplicated(this.cluster.members.accessions)]

      # save ECF sequences of all members
      this.cluster.aln <- extract.aln.subset(this.ecf.aln, this.cluster.members)
      this.cluster.aln.outfile <- file.path(this.ecf.dir, sprintf("cl%d-ecf.aln", ci))
      writeXStringSet(this.cluster.aln, filepath = this.cluster.aln.outfile)
      
      this.cluster.reprseqs.index <- match(this.cluster.members, reprseq.seqs.n)
      this.cluster.reprseqs.seqs <- reprseq.seqs[this.cluster.reprseqs.index]
      this.cluster.reprseqs.seqs.outfile <- file.path(this.ecf.dir, sprintf("cl%d-reprseqs.fa", ci))
      writeXStringSet(this.cluster.reprseqs.seqs, filepath = this.cluster.reprseqs.seqs.outfile)
      
      # Prepare gggenomes inputs
      gggenomes.inputs <- prepare.gggenomes.tables(reprseq.names = this.cluster.members,
                                                   ecf.coords.dt = this.ecf.coords.table,
                                                   metadata.table,
                                                   refseq.protein.stats.table,
                                                   reprseq.table,
                                                   genomes.aba.table,
                                                   reprseq.functions.unique,
                                                   no.genes.upstream.downstream = NO.GENES.TO.EXTEND)
      
      gggenomes.inputs.zoom <- prepare.gggenomes.tables.zoomed(gggenomes.inputs, pad = PAD.PARM)
      
      this.cluster.synteny.plot <- plot.ecf.synteny(gggenomes.inputs.zoom, palette = "Set3")
      this.cluster.synteny.plot.outfile <- file.path(this.ecf.dir, sprintf("cl%d-plot.pdf", ci))
      ggsave(this.cluster.synteny.plot.outfile, this.cluster.synteny.plot, width = 8, height = 5)
    }
    cat("Done!\n\n")
  } else{
    cat("No sufficiently large clusters found. Moving on!\n\n")
  }
}  


