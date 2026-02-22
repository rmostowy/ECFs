### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

source("functions.R")

tables.dir <- "tables"
make.dir(tables.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$qname)

# Load REM table
rem.table <- fread(rem.table.output)
rem.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load PFAM domain table
pfam.domains <- fread(pfam.domains.output)
setkey(pfam.domains, accession)

# Load ECF families
ecf.family.table <- fread(ecf.family.table.output)
setkey(ecf.family.table, ecf)

# Load ECF table
ecf.class.table <- fread(ecf.classification.table.output)
setkey(ecf.class.table, ecf)

# Translate ECFs to ECODs
ecf.class.table$family <- ecf.family.table[.(ecf.class.table$ecf)]$family
ecf2ecod <- ecf.class.table %>% select(ecf, family, t.id)
ecod.tid2tname <- unique(ecod.domains %>% select(t.id = t_id, t.name = t_name))
setkey(ecod.tid2tname, t.id)
ecf2ecod$t.name <- ecod.tid2tname[.(ecf2ecod$t.id)]$t.name
ecf2ecod[is.na(t.name), t.name := "unknown"]
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
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)

# Load reprseq vs Pfam
reprseq.pfam.all <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.pfam.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.pfam, qname)


rec.domain.overlap.table.list <- mclapply(1:nrow(rem.table), function(rem.index){
# rec.domain.overlap.table.list <- mclapply(1:1000, function(rem.index){
  cat(rem.index,"\n")
  this.rr <- rem.table[rem.index]
  this.rr.pair.index <- this.rr$pair.index.rem

  # overlaps with protein1
  protein1.name <- this.rr$qname
  protein1.length <- this.rr$qlength
  aln.start <- this.rr$qstart
  aln.end <- this.rr$qend
  
  # ECOD
  protein1.ecod.hits <- reprseq.ecod[.(protein1.name)]
  protein1.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln.start, aln.end, qstart, qend)]
  setorder(protein1.ecod.hits, -rr.overlap.jaccard)
  protein1.ecod.best.overlap <- protein1.ecod.hits[1]
  protein1.ecod.overlap.score <- get.domain.overlap.score(rec.start = aln.start, rec.end = aln.end,
                                                          domain.start = protein1.ecod.best.overlap$qstart, domain.end = protein1.ecod.best.overlap$qend,
                                                          domain.cov = protein1.ecod.best.overlap$scov)
  protein1.ecod.best.hit <- protein1.ecod.best.overlap$sname
  
  # PFAM
  protein1.pfam.hits <- reprseq.pfam[.(protein1.name)]
  protein1.pfam.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln.start, aln.end, qstart, qend)]
  setorder(protein1.pfam.hits, -rr.overlap.jaccard)
  protein1.pfam.best.overlap <- protein1.pfam.hits[1]
  protein1.pfam.overlap.score <- get.domain.overlap.score(rec.start = aln.start, rec.end = aln.end,
                                                          domain.start = protein1.pfam.best.overlap$qstart, domain.end = protein1.pfam.best.overlap$qend,
                                                          domain.cov = protein1.pfam.best.overlap$scov)
  protein1.pfam.best.hit <- protein1.pfam.best.overlap$sname
  
  # ECF
  protein1.ecf.hits <- reprseq.ecfs[.(protein1.name)]
  protein1.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln.start, aln.end, qstart, qend)]
  setorder(protein1.ecf.hits, -rr.overlap.jaccard)
  protein1.ecf.best.overlap <- protein1.ecf.hits[1]
  protein1.ecf.overlap.score <- get.domain.overlap.score(rec.start = aln.start, rec.end = aln.end,
                                                         domain.start = protein1.ecf.best.overlap$qstart, domain.end = protein1.ecf.best.overlap$qend,
                                                         domain.cov = protein1.ecf.best.overlap$scov)
  protein1.ecf.best.hit <- protein1.ecf.best.overlap$sname
  
  
  # overlaps with protein2
  protein2.name <- this.rr$sname
  protein2.length <- this.rr$slength
  aln.start <- this.rr$sstart
  aln.end <- this.rr$send
  
  # ECOD
  protein2.ecod.hits <- reprseq.ecod[.(protein2.name)]
  protein2.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln.start, aln.end, qstart, qend)]
  setorder(protein2.ecod.hits, -rr.overlap.jaccard)
  protein2.ecod.best.overlap <- protein2.ecod.hits[1]
  protein2.ecod.overlap.score <- get.domain.overlap.score(rec.start = aln.start, rec.end = aln.end,
                                                          domain.start = protein2.ecod.best.overlap$qstart, domain.end = protein2.ecod.best.overlap$qend,
                                                          domain.cov = protein2.ecod.best.overlap$scov)
  protein2.ecod.best.hit <- protein2.ecod.best.overlap$sname
  
  # PFAM
  protein2.pfam.hits <- reprseq.pfam[.(protein2.name)]
  protein2.pfam.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln.start, aln.end, qstart, qend)]
  setorder(protein2.pfam.hits, -rr.overlap.jaccard)
  protein2.pfam.best.overlap <- protein2.pfam.hits[1]
  protein2.pfam.overlap.score <- get.domain.overlap.score(rec.start = aln.start, rec.end = aln.end,
                                                          domain.start = protein2.pfam.best.overlap$qstart, domain.end = protein2.pfam.best.overlap$qend,
                                                          domain.cov = protein2.pfam.best.overlap$scov)
  protein2.pfam.best.hit <- protein2.pfam.best.overlap$sname
  
  # ECF
  protein2.ecf.hits <- reprseq.ecfs[.(protein2.name)]
  protein2.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln.start, aln.end, qstart, qend)]
  setorder(protein2.ecf.hits, -rr.overlap.jaccard)
  protein2.ecf.best.overlap <- protein2.ecf.hits[1]
  protein2.ecf.overlap.score <- get.domain.overlap.score(rec.start = aln.start, rec.end = aln.end,
                                                         domain.start = protein2.ecf.best.overlap$qstart, domain.end = protein2.ecf.best.overlap$qend,
                                                         domain.cov = protein2.ecf.best.overlap$scov)
  protein2.ecf.best.hit <- protein2.ecf.best.overlap$sname
  
  this.rr.out <- data.table(rem.index, this.rr.pair.index, 
                            protein1.name, 
                            protein1.ecod.best = protein1.ecod.best.hit, protein1.ecod.score = round(protein1.ecod.overlap.score, d=2),
                            protein1.pfam.best = protein1.pfam.best.hit, protein1.pfam.score = round(protein1.pfam.overlap.score, d=2),
                            protein1.ecf.best = protein1.ecf.best.hit, protein1.ecf.score = round(protein1.ecf.overlap.score, d=2),
                            protein2.name, 
                            protein2.ecod.best = protein2.ecod.best.hit, protein2.ecod.score = round(protein2.ecod.overlap.score, d=2),
                            protein2.pfam.best = protein2.pfam.best.hit, protein2.pfam.score = round(protein2.pfam.overlap.score, d=2),
                            protein2.ecf.best = protein2.ecf.best.hit, protein2.ecf.score = round(protein2.ecf.overlap.score, d=2)
  )
  return(this.rr.out)
}, mc.cores = 18)

rec.domain.overlap.table <- rbindlist(rec.domain.overlap.table.list)
rec.domain.overlap.table.outfile <- file.path(tables.dir, "rem-domain-overlaps.txt")
fwrite(rec.domain.overlap.table, file = rec.domain.overlap.table.outfile)

