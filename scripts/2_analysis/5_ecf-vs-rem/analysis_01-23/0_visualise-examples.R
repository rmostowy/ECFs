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

# Load visualisation functions
source("functions.R")

output.dir <- file.path(visual.dir, "part5_rem", "all-pairs-sub")
make.dir(output.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$reprseq)
setkey(reprseq.table, reprseq)

# # Calculate NCBI overlaps per reprseq
# max.overlap.per.reprseq <- reprseq.table[, .(reprseq, ncbi.overlap)]
# setorder(max.overlap.per.reprseq, reprseq, -ncbi.overlap)
# max.overlap.per.reprseq <- max.overlap.per.reprseq[, .SD[1], by = reprseq]
# setkey(max.overlap.per.reprseq, reprseq)

# Load recent recombination table
rem.table <- fread(rem.table.output)
rem.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

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
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)

rem.table[, rem.index := .I]
rrr <- rem.table[pident >= 70]
no.pairs <- nrow(rrr)

# 1: visualse direct comparison of the two proteins
for(index in 1:no.pairs){
  cat(index,"\n")
  this.rem.index <- rrr$rem.index[index]
  this.pair.output.dir <- file.path(output.dir, sprintf("rem%.*d", nchar(no.pairs), this.rem.index))
  make.dir(this.pair.output.dir)
  protein1_name <- rrr$qname[index]
  protein2_name <- rrr$sname[index]
  protein1_conf <- reprseq.table[.(protein1_name)]$prodigal.conf/100
  protein2_conf <- reprseq.table[.(protein2_name)]$prodigal.conf/100
  protein1_length <- rrr$qlength[index]
  protein2_length <- rrr$slength[index]
  aln1_start <- rrr$qstart[index]
  aln1_end <- rrr$qend[index]
  aln2_start <- rrr$sstart[index]
  aln2_end <- rrr$send[index]
  aln_identity <- rrr$pident[index]
  protein1_description <- sprintf("F: %s, C: %.2f", reprseq.functions[.(protein1_name)]$funct, protein1_conf)
  protein2_description <- sprintf("F: %s, C: %.2f", reprseq.functions[.(protein2_name)]$funct, protein2_conf)
  
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
  
  p.domains1 <- plot_protein_domains_compact(protein1_length, aln1_start, aln1_end, protein1.domain.hits, protein1_name)
  p.domains2 <- plot_protein_domains_compact(protein2_length, aln2_start, aln2_end, protein2.domain.hits, protein2_name)
  
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

