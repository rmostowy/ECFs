### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Set up output directory
#################################
this.part.dir <- file.path(visual.dir, "part2_containment")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. ECF tables
### 2. ECF comparsion tables
#################################

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

# ECF comparison tables
ecfs.vs.ecfs <- fread(ecfs.vs.ecfs.output) # raw
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output) # processed

#################################
### Process loaded tables
#################################

all.ecfs <- unique(ecf.table$ecf.name)
setkey(ecf.family.table, ecf)
setkey(ecf.table, ecf.name)

ecfs.vs.ecfs[, max.cov := pmax(qcov, scov)]
setkey(ecfs.vs.ecfs, qname, sname)

ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]
eve.conf <- ecfs.vs.ecfs.cov[prob >= CONTAINMENT.PROB.THR]

# generate network table
ecf.network.ab <- eve.conf %>% select(query, subject, prob, qcov, scov)
ecf.network.ba <- eve.conf %>% select(query = subject, subject = query, prob, qcov = scov, scov = qcov)
ecf.network <- rbind(ecf.network.ab, ecf.network.ba)

#################################
### 1. Generate per-ECF 
### summary table
#################################

dt <- copy(ecf.network)
dt[, subject.ef := ecf.family.table[.(subject)]$family]
ecf.no.fams <- dt[, .(no.fam = uniqueN(subject.ef)), by = query]

ecf.connectivity.stats <- dt[, .(
  mean.qcov = mean(qcov),
  mean.scov = mean(scov),
  med.qcov = median(qcov),
  med.scov = median(scov),
  degree = .N
  ), by = query]
ecf.connectivity.stats[, family := ecf.family.table[.(ecf.connectivity.stats$query)]$family]
ecf.connectivity.stats[, length := ecf.table[.(ecf.connectivity.stats$query)]$length]
ecf.connectivity.stats[, asym_mean := mean.qcov - mean.scov]
ecf.connectivity.stats[, asym_med := med.qcov - med.scov]
ecf.connectivity.stats[, strength_mean := (mean.qcov + mean.scov)/2]
ecf.connectivity.stats[, strength_med := (med.qcov + med.scov)/2]

ecf.connectivity.stats <- merge(ecf.connectivity.stats, ecf.no.fams, by = "query")

#################################
### 2. Show qcov/scov heatmap
#################################

p.cov.heatmap <- ggplot(eve.conf, aes(x = qcov, y = scov)) +
  geom_bin2d(bins = 80) +
  theme_bw()
p.cov.heatmap.outfile <- file.path(this.part.dir, "plot1_cov-heatmap.png")
ggsave(p.cov.heatmap.outfile, p.cov.heatmap, width = 6, height = 5)

#################################
### 3. Show qcov/scov heatmap
### for per-ECF averages
#################################

p.connect.cov.heatmap <- ggplot(ecf.connectivity.stats, aes(x = mean.qcov, y = mean.scov)) +
  geom_bin2d(bins = 80) +
  theme_bw()
p.connect.cov.heatmap.outfile <- file.path(this.part.dir, "plot2_connect-cov-heatmap.png")
ggsave(p.connect.cov.heatmap.outfile, p.connect.cov.heatmap, width = 6, height = 5)

#################################
### 4. Asymmetry distribution
#################################

p.asym.dist <- ggplot(ecf.connectivity.stats, aes(x = asym_med)) +
  geom_histogram(bins = 60, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median qcov − median scov)",
    y = "Number of ECFs",
    title = "Distribution of hierarchical asymmetry"
  )

ggsave(file.path(this.part.dir, "plot3_asym-dist.png"),
       p.asym.dist, width = 6, height = 4)

#################################
### 4. Asymmetry vs. Length
#################################

p.asym.len <- ggplot(ecf.connectivity.stats,
                     aes(x = asym_med, y = length)) +
  geom_bin2d(bins = 70) +
  scale_y_continuous(trans = "log10") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median qcov − median scov)",
    y = "ECF length (aa, log10 scale)",
    title = "Asymmetry relates to fragment length"
  )

ggsave(file.path(this.part.dir, "plot4_asym-vs-length.png"),
       p.asym.len, width = 6, height = 5)

#################################
### 5. Asymmetry vs. Strength
#################################

p.asym.str <- ggplot(ecf.connectivity.stats,
                     aes(x = asym_med, y = strength_med)) +
  geom_bin2d(bins = 70) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median qcov − median scov)",
    y = "Median strength ((qcov + scov)/2)",
    title = "Hierarchical role vs match strength"
  )

ggsave(file.path(this.part.dir, "plot5_asym-vs-strength.png"),
       p.asym.str, width = 6, height = 5)

#################################
### 6. Median qcov vs. scov
### clean diagonal view
#################################

p.qs.node <- ggplot(ecf.connectivity.stats,
                    aes(x = med.qcov, y = med.scov)) +
  geom_bin2d(bins = 70) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Median qcov",
    y = "Median scov",
    title = "Per-ECF similarity geometry"
  )

ggsave(file.path(this.part.dir, "plot6_node-qcov-scov.png"),
       p.qs.node, width = 6, height = 5)


#################################
### 7. Asymmetry vs. Degree
#################################

p.asym.deg <- ggplot(ecf.connectivity.stats,
                     aes(x = asym_med, y = degree)) +
  geom_bin2d(bins = 70) +
  scale_y_continuous(trans = "log10") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median)",
    y = "Degree (log10 scale)",
    title = "Hierarchy vs connectivity"
  )

ggsave(file.path(this.part.dir, "plot7_asym-vs-degree.png"),
       p.asym.deg, width = 6, height = 5)

#################################
### 8. Asymmetry vs. Neighbour
### family count
#################################

p.asym.fam <- ggplot(ecf.connectivity.stats,
                     aes(x = asym_med, y = no.fam)) +
  geom_bin2d(bins = 60) +
  scale_y_continuous() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median)",
    y = "Number of distinct neighbor families (log10 scale)",
    title = "Hierarchy vs family diversity"
  )

ggsave(file.path(this.part.dir, "plot8_asym-vs-family-div.png"),
       p.asym.fam, width = 6, height = 5)

#################################
### 9. Strength vs. Length
#################################

p.str.len <- ggplot(ecf.connectivity.stats,
                    aes(x = strength_med, y = length)) +
  geom_bin2d(bins = 70) +
  scale_y_continuous(trans = "log10") +
  theme_bw(base_size = 14) +
  labs(
    x = "Median strength",
    y = "ECF length (log10 scale)",
    title = "Match strength vs fragment size"
  )

ggsave(file.path(this.part.dir, "plot9_strength-vs-length.png"),
       p.str.len, width = 6, height = 5)

#################################
### 10. Containment network
#################################

containment.network <- ecf.network[qcov >= DOMAIN.CONTAINMENT.COV.THR]
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
  neighbours <- neighbors(containment.graph, this.node, mode = "out")
  no.neighbours <- length(neighbours)
  if(no.neighbours >= 2){
    neighbour.pairs <- combn(neighbours$name, 2, simplify = FALSE)
    neighbour.pairs.dt <- as.data.table(do.call(rbind, neighbour.pairs))
    no.pairs <- nrow(neighbour.pairs.dt)
    neighbour.pair.weigths <- ecfs.vs.ecfs[.(neighbour.pairs.dt)]$max.cov
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

p.qij.fn <- file.path(this.part.dir, "plot10_neighbour-cov.png")
ggsave(p.qij.fn, p.qij, width = 7, height = 3)


############
# 
# 
# set.seed(1)
# all_nodes <- unique(c(ecfs.vs.ecfs.cov$query, ecfs.vs.ecfs.cov$subject))
# setkey(ecfs.vs.ecfs.cov, query, subject)
# 
# cont <- dt[qcov >= DOMAIN.CONTAINMENT.COV.THR,
#            .(ecf = query, container = subject,
#              container_fam = ecf.family.table[.(subject)]$family)]
# 
# # observed: per-ecf container lists
# cont_list <- cont[, .(containers = list(container)), by = ecf]
# cont_list[, k := lengths(containers)]
# cont_list <- cont_list[k >= 2]
# 
# get_pair_covs <- function(nodes){
#   ij <- combn(nodes, 2)
#   q <- ij[1,]; s <- ij[2,]
#   cov <- ecfs.vs.ecfs.cov[.(q, s)]$cov.max
#   # cov2 <- ecfs.vs.ecfs.cov[.(s, q)]$cov.max
#   # cov <- fifelse(is.na(cov), cov2, cov)
#   cov[is.na(cov)] <- 0
#   cov
# }
# 
# # observed distribution
# obs <- cont_list[, .(cov = get_pair_covs(containers[[1]])), by = ecf]
# obs[, type := "observed"]
# 
# # null distribution (one replicate; wrap in replicate(B, ...) if you want)
# nul <- cont_list[, {
#   nodes <- sample(all_nodes, k[1], replace = FALSE)
#   .(cov = get_pair_covs(nodes))
# }, by = ecf]
# nul[, type := "null"]
# 
# both <- rbind(obs, nul, fill=TRUE)
# 
# ggplot(both, aes(cov, fill = type)) +
#   geom_histogram(position="identity", bins=40, alpha=0.5) +
#   theme_bw()
