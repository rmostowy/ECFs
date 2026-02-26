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


results.dir <- file.path(visual.dir, "part2_modularity")
make.dir(results.dir)

# Load ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)
all.ecfs <- unique(ecf.table$ecf.name)
setkey(ecf.family.table, ecf)
setkey(ecf.table, ecf.name)

reprseq.vs.ecfs <- fread(reprseq.ecfs.output)

# Calculate coverage
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]

ecfs.vs.ecfs <- fread(ecfs.vs.ecfs.output)

# Make a plot of qcov/scov values from ECF vs ECF table
ecfs.vs.ecfs.cov.high.prob <- ecfs.vs.ecfs.cov[prob >= 99]
p.qscov <- ggplot(ecfs.vs.ecfs.cov.high.prob, aes(x = qcov, y = scov)) + geom_bin2d(bins = 100)
p.qscov.fn <- file.path(results.dir, "plot1_ecf-certain-qscov.png")
# ggsave(p.qscov.fn, p.qscov, width = 6, height = 5)


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
containment.network <- all.edges[query.contained == T] %>% select(from = query, to = subject, qcov, scov)
# convert into igraph object
containment.graph <- graph_from_data_frame(containment.network, directed = TRUE)
containment.graph.nodes <- V(containment.graph)$name
no.nodes.cg <- length(containment.graph.nodes)


###########################
## Plot one example:
## ecf_00002
###########################

this.ecf <- "ecf_00002"

# safety: ensure the node exists
stopifnot(this.ecf %in% V(containment.graph)$name)

# ego nodes = ecf plus all 1-step neighbors (in + out)
ego_nodes <- unique(c(
  this.ecf,
  neighbors(containment.graph, this.ecf, mode = "out")$name,
  neighbors(containment.graph, this.ecf, mode = "in")$name
))

subg <- induced_subgraph(containment.graph, vids = ego_nodes)

# node table
nodes_dt <- data.table(ecf = V(subg)$name)
nodes_dt[, family := ecf.family.table[.(nodes_dt$ecf), family]]
# nodes_dt[is.na(family), family := "unassigned"]

# classify nodes relative to focal ECF
out_n <- neighbors(subg, this.ecf, mode = "out")$name
in_n  <- neighbors(subg, this.ecf, mode = "in")$name
nodes_dt[, role := fifelse(ecf == this.ecf, "focal",
                           fifelse(ecf %in% out_n & ecf %in% in_n, "both",
                                   fifelse(ecf %in% out_n, "contains_focal",
                                           fifelse(ecf %in% in_n, "contained_in_focal", "other"))))]

# put node attrs into igraph
V(subg)$family <- nodes_dt$family[match(V(subg)$name, nodes_dt$ecf)]
V(subg)$role   <- nodes_dt$role[match(V(subg)$name, nodes_dt$ecf)]


p <- ggraph(subg, layout = "stress") +
  geom_edge_link(
    aes(alpha = scov),
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(3, "mm"),
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(color = family,
        shape = role,
        size  = role),
    stroke = 0.6
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3.2
  ) +
  scale_size_manual(
    values = c(
      focal = 6,
      both = 4,
      contains_focal = 3.5,
      contained_in_focal = 3.5,
      other = 3
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      focal = 21,
      both = 22,
      contains_focal = 19,
      contained_in_focal = 19,
      other = 19
    )
  ) +
  guides(
    alpha = guide_legend(title = "Containment\n(coverage)"),
    color = guide_legend(title = "EF family"),
    shape = guide_legend(title = "Node role")
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  ggtitle(sprintf("Containment neighborhood of %s", this.ecf))

p








# 1. Find all pairs of ECFs that contain an ECF contained in both of them
# 2. Calculate cov.max between all of them

contain.neighbour.cov.max.list <- lapply(1:no.nodes.cg, function(node.index){
  # cat(node.index,"\n")
  out <- NULL
  this.node <- containment.graph.nodes[node.index]
  
  # this.node.transitivity <- full.network.tr[which(names(full.network.tr)==this.node)]
  
  # 2. Identify all possible pairwise connections between neighbours in the containment network
  neighbours <- neighbors(containment.graph, this.node)
  # ecf.family.table[.(names(neighbours))]
  
  
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

# p.qij.fn <- file.path(results.dir, "plot2_transitivity-contained.png")
# ggsave(p.qij.fn, p.qij, width = 7, height = 3)


## Count the number of 

dt <- copy(all.edges)
dt[, query.contained := NULL]

ecf.connectivity.stats <- dt[, .(
  mean.qcov = mean(qcov),
  mean.scov = mean(scov),
  n = .N
  ), by = query]

ecf.connectivity.stats[mean.qcov > 0.9 & n >= 5]


# contain.neighbour.cov.max.list <- lapply(1:no.nodes.cg, function(node.index){
#   # cat(node.index,"\n")
#   out <- NULL
#   this.node <- containment.graph.nodes[node.index]
#   
#   # this.node.transitivity <- full.network.tr[which(names(full.network.tr)==this.node)]
#   
#   # 2. Identify all possible pairwise connections between neighbours in the containment network
#   neighbours <- neighbors(containment.graph, this.node)
#   ecf.family.table[.(names(neighbours))]
  
  



## plot ecf_00001 example



# ########################
# ### ARCHITECTURE BRAINSTORM
# ########################


reprseq.vs.ecfs$EF <- ecf.family.table[.(reprseq.vs.ecfs$sname)]$family
reprseqs.vs.efs <- reprseq.vs.ecfs[
  , .(EFs = list(sort(unique(EF)))), 
  by = qname
]
setkey(reprseqs.vs.efs, qname)


this.ecf <- "ecf_00001"
this.ef <- ecf.family.table[ecf == this.ecf]$family
this.ecf.vs.reprseqs.table <- reprseq.vs.ecfs[sname == this.ecf]
this.ecf.vs.reprseqs <- unique(this.ecf.vs.reprseqs.table$qname)


this.ecf.vs.reprseqs.all.efs <- reprseqs.vs.efs[.(this.ecf.vs.reprseqs)]


this.ecf.vs.reprseqs.table[qname == this.ecf.vs.reprseqs[1]]

##########

# Show EF1 is seen with many different partner EF families across proteins.

DT <- copy(reprseq.vs.ecfs)   # must already have EF assigned
EF1 <- "EF0110"

proteins <- unique(DT[EF == EF1, qname])

# per protein: how many partner families (excluding EF1)
per_protein_partners <- DT[qname %in% proteins,
                           .(partners = uniqueN(EF[EF != EF1])),
                           by = qname
]

standalone_frac <- per_protein_partners[partners == 0, .N] / nrow(per_protein_partners)

# overall partner diversity + top partners
partner_counts <- DT[qname %in% proteins & EF != EF1,
                     .(n_proteins = uniqueN(qname)),
                     by = EF
][order(-n_proteins)]

list(
  n_proteins = length(proteins),
  standalone_frac = standalone_frac,
  n_distinct_partners = nrow(partner_counts),
  top_partners = partner_counts[1:10]
)

###########

# “Is EF1 a subdomain or a whole-domain in practice?”

EF1 <- "EF0001"

cov <- DT[EF == EF1,
          .(max_qcov = max(qcov, na.rm=TRUE)),
          by = qname
]

cov[, .(
  n = .N,
  median = median(max_qcov),
  mean = mean(max_qcov),
  frac_full = mean(max_qcov > 0.9),
  frac_small = mean(max_qcov < 0.3)
)]


###########


DT <- copy(reprseq.vs.ecfs)   # must already have EF assigned
EF1 <- "EF0001"

proteins <- unique(DT[EF == EF1, qname])

# per protein: how many partner families (excluding EF1)
per_protein_partners <- DT[qname %in% proteins,
                           .(partners = uniqueN(EF[EF != EF1])),
                           by = qname
]

standalone_frac <- per_protein_partners[partners == 0, .N] / nrow(per_protein_partners)

# overall partner diversity + top partners
partner_counts <- DT[qname %in% proteins & EF != EF1,
                     .(n_proteins = uniqueN(qname)),
                     by = EF
][order(-n_proteins)]

list(
  n_proteins = length(proteins),
  standalone_frac = standalone_frac,
  n_distinct_partners = nrow(partner_counts),
  top_partners = partner_counts[1:10]
)


#########


EF1 <- "EF0001"
proteins <- unique(DT[EF == EF1, qname])

bg <- DT[, .(bg_n = uniqueN(qname)), by = EF]
obs <- DT[qname %in% proteins, .(obs_n = uniqueN(qname)), by = EF]

enr <- merge(obs, bg, by="EF", all=TRUE)
enr[is.na(obs_n), obs_n := 0L]

enr[, `:=`(
  p_obs = obs_n / length(proteins),
  p_bg  = bg_n  / uniqueN(DT$qname),
  enrich = (obs_n / length(proteins)) / (bg_n / uniqueN(DT$qname))
)]

enr[EF != EF1][order(-enrich)][1:20]




# this.ecf <- "ecf_00001"
# 
# DT <- copy(reprseq.vs.ecfs)
# # ensure EF is present
# DT[, EF := ecf.family.table[.(sname)]$family]
# 
# # 1) ecf_00001 footprint per qname
# foot <- DT[sname == this.ecf,
#            .(qname,
#              ecf_start = qstart,
#              ecf_end   = qend,
#              ecf_len   = qend - qstart + 1L)]
# setkey(foot, qname)
# 
# # 2) all other ECF hits in those proteins
# others <- DT[foot, on = "qname"][sname != this.ecf]
# 
# # 3) overlap with ecf_00001 footprint
# others[, `:=`(
#   hit_len = qend - qstart + 1L,
#   ovl = pmax(0L, pmin(qend, ecf_end) - pmax(qstart, ecf_start) + 1L)
# )]
# others[, `:=`(
#   frac_of_ecf   = ovl / ecf_len,
#   frac_of_other = ovl / hit_len,
#   jaccard       = fifelse(ovl > 0,
#                           ovl / (ecf_len + hit_len - ovl),
#                           0)
# )]
# 
# # 4) summarize per EF partner (pick thresholds you like)
# partner.summary <- others[
#   , .(
#     n_hits = .N,
#     n_proteins = uniqueN(qname),
#     n_overlap_any = uniqueN(qname[ovl > 0]),
#     n_overlap_50pc_ecf = uniqueN(qname[frac_of_ecf >= 0.5]),
#     mean_frac_ecf = mean(frac_of_ecf),
#     median_frac_ecf = median(frac_of_ecf),
#     mean_jaccard = mean(jaccard)
#   ),
#   by = EF
# ][order(-n_overlap_any, -mean_frac_ecf)]
# 
# partner.summary
# 
# per.protein <- others[
#   , .(
#     max_frac_ecf = max(frac_of_ecf),
#     n_partners = uniqueN(EF)
#   ),
#   by = qname
# ]
# 
# per.protein[, architecture :=
#               fifelse(max_frac_ecf >= 0.8, "contained",
#                       fifelse(max_frac_ecf >= 0.2, "partial",
#                               "standalone"))
# ]
# 
# focal_EF <- "EF0001"
# 
# per.protein <- others[
#   , .(
#     n_partners_all = uniqueN(EF),
#     n_partners_other = uniqueN(EF[EF != focal_EF]),
#     max_frac_ecf_other = if (.N == 0L) 0 else max(frac_of_ecf[EF != focal_EF], na.rm = TRUE)
#   ),
#   by = qname
# ]
# 
# # handle proteins where there are no "other" EFs (max() on empty -> -Inf)
# per.protein[!is.finite(max_frac_ecf_other), max_frac_ecf_other := 0]
# 
# per.protein[, architecture :=
#               fifelse(n_partners_other == 0, "standalone",
#                       fifelse(max_frac_ecf_other >= 0.8, "contained",
#                               fifelse(max_frac_ecf_other >= 0.2, "partial",
#                                       "cooccur_no_overlap")))
# ]
# 
# per.protein[, .N, by = architecture][order(-N)]
# 
# top.partner <- others[EF != focal_EF][
#   order(qname, -frac_of_ecf)
# ][
#   , .SD[1], by = qname
# ][
#   , .(qname, top_partner_EF = EF, top_partner_frac = frac_of_ecf)
# ]
# 
# per.protein <- top.partner[per.protein, on="qname"]


# ########################
# ### GPT ADD ON
# ########################
# 
# 
# ecf2fam <- ecf.family.table
# setkey(ecf2fam, ecf)
# 
# 
# # conservative coverage on each pair (mutual-ish)
# # (use your already-computed qcov/scov; if not present, compute as you do)
# ecfs.vs.ecfs.cov[, cov.min := pmin(qcov, scov)]
# ecfs.vs.ecfs.cov[, cov.max := pmax(qcov, scov)]
# 
# # fast lookup table for pair weights (undirected-ish key)
# pair_w <- ecfs.vs.ecfs.cov[, .(a = pmin(query, subject),
#                                b = pmax(query, subject),
#                                cov = cov.min,
#                                prob = prob)]
# setkey(pair_w, a, b)
# 
# get_pair_cov <- function(u, v){
#   a <- pmin(u, v); b <- pmax(u, v)
#   out <- pair_w[J(a, b), cov]
#   out[is.na(out)] <- 0
#   out
# }
# 
# ##########
# # p1
# ##########
# 
# nodes <- V(containment.graph)$name
# 
# reuse_dt <- rbindlist(lapply(nodes, function(x){
#   cont <- neighbors(containment.graph, x, mode = "out")$name  # x -> container
#   if (length(cont) < 2) return(NULL)
#   
#   fams <- ecf2fam[J(cont), family]
#   fams <- fams[!is.na(fams)]
#   if (!length(fams)) return(NULL)
#   
#   tab <- sort(table(fams), decreasing = TRUE)
#   data.table(
#     ecf = x,
#     k_containers = length(cont),
#     n_container_families = length(unique(fams)),
#     frac_top_family = as.numeric(tab[1]) / length(fams)
#   )
# }))
# 
# # Plot: how many distinct families do containers come from (conditioned on k>=2)
# p1 <- ggplot(reuse_dt, aes(x = k_containers, y = n_container_families)) +
#   geom_bin2d(bins = 80) +
#   theme_bw() +
#   labs(x = "# containers (outdegree in containment graph)",
#        y = "# distinct container families",
#        title = "Contained fragments recur across multiple ECF families")
# 
# ggsave("p1.png", p1, width = 6, height = 5)
# 
# 
# ##########
# # p1b
# ##########
# 
# reuse_dt[, fam_cat := fifelse(n_container_families == 1, "1",
#                               fifelse(n_container_families == 2, "2", "3+"))]
# 
# p1b <- ggplot(reuse_dt[, .N, by = fam_cat][, prop := N/sum(N)], aes(fam_cat, prop)) +
#   geom_col() + theme_bw() +
#   labs(x = "# container families per fragment", y = "Proportion", title = "How often a fragment is contained in multiple families")
# 
# p1b
# ggsave("p1b.png", p1b, width = 6, height = 5)
# 
# ##########
# # p2
# ##########
# 
