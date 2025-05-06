### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

results.dir <- "analysis2_modularity"
make.dir(results.dir)

# Load ECF tables
ecf.table <- fread(ecf.table.output)
all.ecfs <- unique(ecf.table$ecf.name)
setkey(ecf.table, ecf.name)

# Calculate coverage
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]

# Make a plot of qcov/scov values from ECF vs ECF table
ecfs.vs.ecfs.cov.high.prob <- ecfs.vs.ecfs.cov[prob >= 99.9]
p.qscov <- ggplot(ecfs.vs.ecfs.cov.high.prob, aes(x = qcov, y = scov)) + geom_bin2d(bins = 100)
p.qscov.fn <- file.path(results.dir, "plot1_ecf-certain-qscov.png")
ggsave(p.qscov.fn, p.qscov, width = 6, height = 5)


# Create a full network of ECFs which contains information about 
# the maximum coverage per each pair of ECFs between which a connection
# was found
ecf.network <- ecfs.vs.ecfs.cov %>%
  select(query, subject, weight = cov.max)
setkey(ecf.network, query, subject)

# Create a containment network of ECFs:
# a directed graph which pairs of contained ECFs
# where 'from' ECFs are contained in 'to' ECFs 
edges.q2s <- ecfs.vs.ecfs.cov %>% select(query, subject, qcov, scov)
edges.s2q <- ecfs.vs.ecfs.cov %>% select(query = subject, subject = query, qcov = scov, scov = qcov)
all.edges <- rbind(edges.q2s, edges.s2q)
all.edges[, query.contained := ifelse(qcov >= DOMAIN.CONTAINMENT.COV.THR, TRUE, FALSE)]
containment.network <- all.edges[qcov >= DOMAIN.CONTAINMENT.COV.THR] %>% select(from = query, to = subject)
# convert into igraph object
containment.graph <- graph_from_data_frame(containment.network, directed = TRUE)
containment.graph.nodes <- V(containment.graph)$name
no.nodes.cg <- length(containment.graph.nodes)

# 1. Find all pairs of ECFs that contain an ECF contained in both of them
# 2. Calculate cov.max between all of them

contain.neighbour.cov.max.list <- lapply(1:no.nodes.cg, function(node.index){
  # cat(node.index,"\n")
  out <- NULL
  this.node <- containment.graph.nodes[node.index]
  
  # this.node.transitivity <- full.network.tr[which(names(full.network.tr)==this.node)]
  
  # 2. Identify all possible pairwise connections between neighbours in the containment network
  neighbours <- neighbors(containment.graph, this.node)
  no.neighbours <- length(neighbours)
  if(no.neighbours >= 2){
    neighbour.pairs <- combn(neighbours$name, 2, simplify = FALSE)
    neighbour.pairs.dt <- as.data.table(do.call(rbind, neighbour.pairs))
    no.pairs <- nrow(neighbour.pairs.dt)
    neighbour.pair.weigths <- ecf.network[.(neighbour.pairs.dt)]$weight
    neighbour.pair.weigths <- fcoalesce(neighbour.pair.weigths, 0)
    
    # 5. Return 
    out <- data.table(ecf = this.node, neighbour.cov.max = neighbour.pair.weigths)  
  }
  return(out)
})
contain.neighbour.cov.max <- rbindlist(contain.neighbour.cov.max.list)
contain.neighbour.cov.max.avg <- contain.neighbour.cov.max[
  , .(avg.neighbour.cov.max = mean(neighbour.cov.max, na.rm = TRUE)),
  by = ecf
]
contain.neighbour.cov.max.avg <- merge(contain.neighbour.cov.max.avg, contain.neighbour.cov.max[, .N, by = ecf], by = "ecf")

data2plot <- rbind(data.table(type = "full", value = contain.neighbour.cov.max$neighbour.cov.max),
      data.table(type = "average", value = contain.neighbour.cov.max.avg$avg.neighbour.cov.max))

# Create the histogram with facets and independent Y axes
p.qij <- ggplot(data2plot, aes(x = value)) +
  geom_histogram(bins = 30, fill = "#69b3a2", alpha = 0.7, color = "black") +  # Histogram
  facet_wrap(~ type, scales = "free_y") +  # Independent Y axis per panel
  labs(
    x = "Value",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )

p.qij.fn <- file.path(results.dir, "plot2_transitivity-contained.png")
ggsave(p.qij.fn, p.qij, width = 7, height = 3)

