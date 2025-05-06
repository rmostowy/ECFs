### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load ECF tables
ecf.table <- fread(ecf.table.output)
setkey(ecf.table, ecf.name)
all.ecfs <- unique(ecf.table$ecf.name)
# ecf.lengths <- ecf.table$length

MAX.ECF.LENGTH <- 800
ECF.MIN.COV.THR <- 0.1
ECF.MAX.COV.THR <- 0.9
ecf4network <- ecf.table[length <= MAX.ECF.LENGTH]$ecf.name

# Calculate coverage
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]

ecfs.vs.ecfs.cov.filtered <- ecfs.vs.ecfs.cov[query %chin% ecf4network & subject %chin% ecf4network]
ecfs.vs.ecfs.cov.filtered <- ecfs.vs.ecfs.cov.filtered[cov.min >= ECF.MIN.COV.THR & cov.max >= ECF.MAX.COV.THR]
# ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[!(ecfs.vs.ecfs.cov$qlength*ECF.LENGTH.DIFF.THR < ecfs.vs.ecfs.cov$slength | ecfs.vs.ecfs.cov$qlength/ECF.LENGTH.DIFF.THR > ecfs.vs.ecfs.cov$slength)]
# ecf.network.full.self <- data.table(query = all.ecfs, subject = all.ecfs, prob = 100, pident = 100, hit.length = ecf.lengths$length, qcov = 1, scov = 1, cov.min = 1, cov.max = 1)
# ecf.network.full <- rbind(ecf.network.full.nonself, ecf.network.full.self)

# # filter ECF pairs with an extreme length difference
# filter by min/max coverage

# # create an ECF network
# ecf.network.full.nonself <- ecfs.vs.ecfs.cov %>%
#   select(query, subject, prob, pident, hit.length, qcov, scov, cov.min, cov.max)

ecf.network <- ecfs.vs.ecfs.cov.filtered %>%
  select(query, subject, weight = cov.max)
setkey(ecf.network, query, subject)

ecf.network.outfile <- file.path(ecf.network.dir, "ecf-vis-network.txt")
fwrite(ecf.network, file = ecf.network.outfile, sep = "\t")

# containment.threshold <- 0.9
# edges.q2s <- ecfs.vs.ecfs.cov %>% select(query, subject, qcov, scov)
# edges.s2q <- ecfs.vs.ecfs.cov %>% select(query = subject, subject = query, qcov = scov, scov = qcov)
# all.edges <- rbind(edges.q2s, edges.s2q)
# all.edges[, query.contained := ifelse(qcov >= containment.threshold, TRUE, FALSE)]
# 
# containment.network <- graph_from_data_frame(all.edges[query.contained==TRUE], directed = TRUE)
# full.network <- graph_from_data_frame(ecf.network, directed = FALSE)
# no.nodes.contain <- length(V(containment.network)$name)
# 
# full.network.tr <- transitivity(full.network, type = "local")
# 
# 
# 
# ###########################################
# 
# # 1. Take a query ECF from the containment network
# contain.neighbour.max.cov.list <- lapply(1:no.nodes.contain, function(node.index){
#   # cat(node.index,"\n")
#   out <- NULL
#   this.node <- V(containment.network)$name[node.index]
#   
#   # this.node.transitivity <- full.network.tr[which(names(full.network.tr)==this.node)]
#   
#   # 2. Identify all possible pairwise connections between neighbours in the containment network
#   neighbours <- neighbors(containment.network, this.node)
#   no.neighbours <- length(neighbours)
#   if(no.neighbours >= 2){
#     neighbour.pairs <- combn(neighbours$name, 2, simplify = FALSE)
#     neighbour.pairs.dt <- as.data.table(do.call(rbind, neighbour.pairs))
#     no.pairs <- nrow(neighbour.pairs.dt)
#     neighbour.pair.weigths <- ecf.network[.(neighbour.pairs.dt)]$weight
#     neighbour.pair.weigths <- fcoalesce(neighbour.pair.weigths, 0)
# 
#         # 5. Return 
#     out <- data.table(ecf = this.node, neighbour.max.cov = neighbour.pair.weigths)  
#   }
#   return(out)
# })
# 
# contain.neighbour.max.cov <- rbindlist(contain.neighbour.max.cov.list)
# 
# contain.neighbour.max.cov.mean <- contain.neighbour.max.cov[
#   , .(avg.neighbour.max.cov = mean(neighbour.max.cov, na.rm = TRUE)),
#   by = ecf
# ]
# 
# # contain.neighbour.max.cov.mean$transitivity <- sapply(1:nrow(contain.neighbour.max.cov.mean), function(k) full.network.tr[which(contain.neighbour.max.cov.mean[k]$ecf == names(full.network.tr))])
# 
# 
# 
# ecf.network[query == this.node | subject == this.node]
# 
# 
# node.index=95
# this.node <- V(containment.network)$name[node.index]
# 
# neighbourhood_nodes <- neighborhood(containment.network, order = 2, nodes = this.node)
# subgraph <- induced_subgraph(containment.network, vids = neighbourhood_nodes[[1]])
# plot(
#   subgraph,
#   vertex.color = ifelse(V(subgraph)$name == this.node, "red", "lightblue"),
#   vertex.size = 10,
#   edge.arrow.size = 0.5,
#   main = paste("2nd-Order Neighbourhood of Node:", this.node)
# )
# 
# neighbourhood_nodes_full <- neighborhood(full.network, order = 2, nodes = this.node)
# subgraph_full <- induced_subgraph(full.network, vids = neighbourhood_nodes_full[[1]])
# plot(
#   subgraph_full,
#   vertex.color = ifelse(V(subgraph_full)$name == this.node, "red", "lightblue"),
#   vertex.size = 10,
#   edge.arrow.size = 0.5,
#   main = paste("2nd-Order Neighbourhood of Node:", this.node)
# )
# 
# ecf.network[query == this.node | subject == this.node]
# 
# 
