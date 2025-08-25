### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### locate the RStudio project root and load the configuration file
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

##########################################
# 1. Create a reference table of all
#    proteins used in the study
##########################################

# Load tables
name.table.origin <- fread(name.table.origin.filename)
refseq.protein.clustering <- fread(refseq.protein.clustering.filename, header = F)
refseq.protein.stats.table <- fread(refseq.protein.stats.table.filename)
setkey(refseq.protein.stats.table, name)

# Merge tables to obtain the basic version of name.table
merged.protein.tables <- merge.data.table(name.table.origin, refseq.protein.clustering, by.x = "prot.name", by.y = "V1")
name.table <- merged.protein.tables %>% select(reprseq = repr.name, cluster.name = prot.name, ncbi.id = V2)
setorder(name.table, reprseq, ncbi.id)

# Add gene confidence from Prodigal
name.table$prodigal.conf <- refseq.protein.stats.table[.(name.table$ncbi.id)]$conf

##########################################
# 3. Save reprseq sequences
##########################################

all.reprseqs <- unique(name.table$reprseq)
N <- length(all.reprseqs)

all.cds <- readAAStringSet(all.cds.filename)
all.cds.n <- names(all.cds)

protname2seq <- data.table(prot.name = all.cds.n, index = 1:length(all.cds.n))
setkey(protname2seq, prot.name)

repr2clust <- unique(name.table %>% select(reprseq, cluster.name))
repr2clust$seq.index <- protname2seq[.(repr2clust$cluster.name)]$index

reprseq.seqs <- all.cds[repr2clust$seq.index]
names(reprseq.seqs) <- repr2clust$reprseq

writeXStringSet(reprseq.seqs, filepath = reprseq.seqs.filename)

##########################################
# 2. Calculate reprseq qcov/scov
##########################################

use.cores <- 8
no.iterations <- 20

reprseq.vs.reprseq <- fread(reprseq.vs.reprseq.filename)
reprseq.vs.reprseq.cov <- reprseq.vs.reprseq %>% filter(prob >= MIN.PROB.FOR.REPRSEQ.FAM)

pair.table <- get.all.pairs(reprseq.vs.reprseq.cov)
no.pairs <- nrow(pair.table)

pairs.per.iteration <- floor(no.pairs/no.iterations)
start.end.per.iteration <- lapply(1:no.iterations, function(index){
  end.pair <- pairs.per.iteration*index
  start.pair <- end.pair - pairs.per.iteration + 1
  if(index == no.iterations) end.pair <- no.pairs
  return(c(start.pair, end.pair))
})

pb <- txtProgressBar(min = 0, max = no.iterations, style = 3, initial = 0)
pair.table.hits.list <- lapply(1:no.iterations, function(iteration.index){
  this.iteration.pairs.index <- start.end.per.iteration[[iteration.index]][1]:start.end.per.iteration[[iteration.index]][2]
  this.iteration.pairs <- pair.table[this.iteration.pairs.index]
  this.iteration.hits <- get.pair.table.hits(this.iteration.pairs, reprseq.vs.reprseq.cov, no.cores = use.cores)
  this.iteration.output <- cbind(this.iteration.pairs, this.iteration.hits)
  setTxtProgressBar(pb, iteration.index)
  return(this.iteration.output)
})
close(pb)

reprseq.pair.table.hits <- rbindlist(pair.table.hits.list)
new.col.names <- c("query", "subject", "prob", "scov.min", "scov.max", "qcov.min", "qcov.max", "pident", "hit.length.min", "hit.length.max")
setnames(reprseq.pair.table.hits, new = new.col.names)
reprseq.pair.table.hits[, prob := prob * 100]
reprseq.pair.table.hits[, pident := pident * 100]
setorder(reprseq.pair.table.hits, query, subject)

reprseq.vs.reprseq.cov <- reprseq.pair.table.hits
fwrite(reprseq.vs.reprseq.cov, file = reprseq.vs.reprseq.cov.output)

##########################################
# 3. Run MCL to identify reprseq families
##########################################

# calculate weights for REPRSEQ families
reprseq.vs.reprseq.cov[, qcov := pmin(qcov.min, qcov.max)]
reprseq.vs.reprseq.cov[, scov := pmin(scov.min, scov.max)]
reprseq.vs.reprseq.cov[, weight := scov*qcov*prob/100]
reprseq.nonself.hits <- reprseq.vs.reprseq.cov[qcov >= REPRSEQ.FAMILY.QCOV.THR & scov >= REPRSEQ.FAMILY.SCOV.THR] %>%
  select(from = query, to = subject, weight)
reprseq.self.hits <- data.table(from = name.table$reprseq, to = name.table$reprseq, weight = 1)
reprseq.mcl.table <- rbind(reprseq.nonself.hits, reprseq.self.hits)

# -- run MCL --
# create a temp mcl directory
mcl.run.dir <- file.path(data.tmp.dir, "mcl-reprseq")
mcl.tmp.dir <- file.path(mcl.run.dir, "tmp")
make.dir(mcl.tmp.dir)

# run MCL
input.filename <- file.path(mcl.run.dir, "mcl-input.txt")
fwrite(reprseq.mcl.table, file = input.filename, sep = " ")
output.filename <- file.path(mcl.run.dir, "mcl-output.txt")
run.mcl(input.filename, output.filename, inflation.factor = REPRSEQ.FAMILY.INFLATION.FACTOR, tmp.dir = mcl.tmp.dir)

# process MCL output
mcl.output <- readLines(output.filename)
Nf <- length(mcl.output)

mcl.output.dt <- data.table(index = 1:Nf, mcl_output = mcl.output)
mcl.output.dt <- mcl.output.dt[, .(reprseq = unlist(tstrsplit(mcl_output, "\t", fixed = TRUE))), by = index]
mcl.output.dt[, family := sprintf("fam%.*d", nchar(Nf), index)]

# clear the directory
delete.dir(mcl.run.dir)

# save the output
reprseq.family.table <- mcl.output.dt %>% select(reprseq, family)
setkey(reprseq.family.table, reprseq)

##########################################
# 4. Merge families into name.table
##########################################

# Write the final name table to a file
name.table$family <- reprseq.family.table[.(name.table$reprseq)]$family
fwrite(name.table, file = reprseq.table.output)

