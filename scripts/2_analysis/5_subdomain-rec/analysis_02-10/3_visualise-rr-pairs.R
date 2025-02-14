### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Define output folder
output.dir <- file.path(visual.dir, "analysis_02-10", "rr-network")
make.dir(output.dir)

# Load reprseq table
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Load recent recombination table
recent.rec.table <- fread(reprseq.recent.rec.table.output)
recent.rec.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

recent.rec.table.sub <- fread("rececent-rec-table-sub.txt")
recent.rec.table.sub[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

recent.rec.table.sub20 <- fread("rececent-rec-table-sub-20p.txt")
recent.rec.table.sub20[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

recent.rec.table.sub30 <- fread("rececent-rec-table-sub-30p.txt")
recent.rec.table.sub30[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

recent.rec.table.sub40 <- fread("rececent-rec-table-sub-40p.txt")
recent.rec.table.sub40[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

recent.rec.table.sub70 <- fread("rececent-rec-table-sub-70p.txt")
recent.rec.table.sub70[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load functional annotations
reprseq.functions <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions, reprseq)

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.ecod, qname)

reprseq2ecod <- reprseq.ecod %>% select(qname, sname, bitscore)
setorder(reprseq2ecod, qname, -bitscore)
reprseq2ecod <- reprseq2ecod[, .SD[1], by = .(qname)]
reprseq2ecod$t.id <- ecod.domains[.(reprseq2ecod$sname)]$t_id
reprseq2ecod$t.name <- ecod.domains[.(reprseq2ecod$sname)]$t_name
reprseq2ecod.final <- reprseq2ecod %>% select(reprseq = qname, t.id, t.name)


reprseq.metadata <- reprseq.functions %>% select(reprseq, phrogs.score = score, funct, class, annotation)
reprseq.metadata <- merge(reprseq.metadata, reprseq2ecod.final, by = 'reprseq', all.x = T) 
reprseq.metadata.outfile <- file.path(output.dir, "reprseq-table.txt")
fwrite(reprseq.metadata, file = reprseq.metadata.outfile, sep = "\t")
# setkey(reprseq.metadata, reprseq)

rr.network <- recent.rec.table %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.outfile <- file.path(output.dir, "rr-network.txt")
fwrite(rr.network, file = rr.network.outfile, sep = "\t")

rr.network100 <- recent.rec.table[qlength > 100 & slength > 100] %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.outfile <- file.path(output.dir, "rr-network-long.txt")
fwrite(rr.network100, file = rr.network.outfile, sep = "\t")

rr.network.sub <- recent.rec.table.sub %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.sub.outfile <- file.path(output.dir, "rr-network-sub.txt")
fwrite(rr.network.sub, file = rr.network.sub.outfile, sep = "\t")

rr.network.sub20 <- recent.rec.table.sub20 %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.sub20.outfile <- file.path(output.dir, "rr-network-sub-20p.txt")
fwrite(rr.network.sub20, file = rr.network.sub20.outfile, sep = "\t")

rr.network.sub30 <- recent.rec.table.sub30 %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.sub30.outfile <- file.path(output.dir, "rr-network-sub-30p.txt")
fwrite(rr.network.sub30, file = rr.network.sub30.outfile, sep = "\t")

rr.network.sub40 <- recent.rec.table.sub40 %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.sub40.outfile <- file.path(output.dir, "rr-network-sub-40p.txt")
fwrite(rr.network.sub40, file = rr.network.sub40.outfile, sep = "\t")

rr.network.sub50 <- recent.rec.table.sub50 %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.sub50.outfile <- file.path(output.dir, "rr-network-sub-50p.txt")
fwrite(rr.network.sub50, file = rr.network.sub50.outfile, sep = "\t")

rr.network.sub70 <- recent.rec.table.sub70 %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.sub70.outfile <- file.path(output.dir, "rr-network-sub-70p.txt")
fwrite(rr.network.sub70, file = rr.network.sub70.outfile, sep = "\t")
 

# # #############
# # rr.network$q.funct <- reprseq.metadata[.(rr.network$qname)]$funct
# # rr.network$s.funct <- reprseq.metadata[.(rr.network$sname)]$funct
# # 
# # ab <- grepl("fiber", rr.network$q.funct) & grepl("methyltransf", rr.network$s.funct)
# # ba <- grepl("fiber", rr.network$s.funct) & grepl("methyltransf", rr.network$q.funct)
# # rr.network[ab | ba]
# 
# 
# 
# rr.network.q.funct <- reprseq.functions[.(rr.network$qname)]$funct
# rr.network$q.cat <- "other"
# rr.network$q.cat[rr.network.q.funct %in% c("tail fiber", "tail spike")] <- "RBP"
# rr.network$q.cat[grepl("methyltransferase", rr.network.q.funct) | grepl("kinase", rr.network.q.funct)] <- "MK"
# 
# rr.network.s.funct <- reprseq.functions[.(rr.network$sname)]$funct
# rr.network$s.cat <- "other"
# rr.network$s.cat[rr.network.s.funct %in% c("tail fiber", "tail spike")] <- "RBP"
# rr.network$s.cat[grepl("methyltransferase", rr.network.s.funct) | grepl("kinase", rr.network.s.funct)] <- "MK"
# 
# 
# # Convert rr.network to graph format
# rr.network.graph <- graph_from_data_frame(rr.network[, .(qname, sname)], directed = FALSE)
# 
# # Extract nodes with category RBP and MK
# rbp_nodes <- unique(rr.network[q.cat == "RBP", qname])
# mk_nodes <- unique(rr.network[s.cat == "MK", sname])
# 
# # Find the shortest path between any RBP and any MK
# shortest_paths_matrix <- shortest.paths(rr.network.graph, v = rbp_nodes, to = mk_nodes)
# 
# # Get the minimum shortest path between any RBP and MK
# min_distance <- min(shortest_paths_matrix, na.rm = TRUE)
# 
# # Find which pairs have this shortest path
# shortest_path_pairs <- which(shortest_paths_matrix == min_distance, arr.ind = TRUE)
# 
# # Extract actual protein names
# rbp_proteins <- rownames(shortest_paths_matrix)[shortest_path_pairs[, 1]]
# mk_proteins <- colnames(shortest_paths_matrix)[shortest_path_pairs[, 2]]
# 
# # Create a data.table of shortest paths
# shortest_paths_dt <- data.table(
#   rbp_protein = rbp_proteins,
#   mk_protein = mk_proteins,
#   shortest_distance = min_distance
# )
