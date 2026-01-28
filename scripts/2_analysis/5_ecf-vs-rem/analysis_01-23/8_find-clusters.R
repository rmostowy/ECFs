### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))
# suppressWarnings(suppressPackageStartupMessages(library(gggenes)))
# suppressWarnings(suppressPackageStartupMessages(library(ggtree)))
# suppressWarnings(suppressPackageStartupMessages(library(tibble)))
# suppressWarnings(suppressPackageStartupMessages(library(ggtreeExtra)))
# suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
# suppressWarnings(suppressPackageStartupMessages(library(ape)))
# suppressWarnings(suppressPackageStartupMessages(library(phangorn)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

source("functions.R")



RE.MAX.DISTANCE=0.5
RE.MIN.CLUSTER.SIZE=5


# Define output folder
table.dir <- "tables"
output.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs")
ecf.clusters.dir <- file.path(output.dir, "ecf-clusters")
make.dir(ecf.clusters.dir)

reprseq.vs.reprseq.blast <- fread(reprseq.vs.reprseq.blast.filtered.unordered.output)
setkey(reprseq.vs.reprseq.blast, qname, sname)


all.ecfs <- list.files(ecf.subseqs.dir)
all.ecfs <- unique(substr(all.ecfs, 1, 9))


for(index in 1:length(all.ecfs)){
  ecf.name <- all.ecfs[index]
  this.ecf.dir <- file.path(ecf.clusters.dir, ecf.name)
  make.dir(this.ecf.dir)
  # if(!file.exists(ecf.tree.filename)){}
  tree_file <- file.path(ecf.subseqs.dir, sprintf("%s.aln.treefile", ecf.name))
  aln_file <- file.path(ecf.subseqs.dir, sprintf("%s.aln", ecf.name))
  aln <- readAAStringSet(aln_file)  # must be aligned
  
  # 1. Calculate distances between all RE-ECFs, considering min cov of 80%
  D <- DistanceMatrix(aln, minCoverage = 0.8, verbose = F)
  D[is.na(D)] <- 1
  hc <- hclust(as.dist(D), method = "complete")
  
  # 2. Identify clusters of sequences which have a distance less than 50%
  cl <- cutree(hc, h = RE.MAX.DISTANCE)
  cl.sizes <- table(cl)
  no.suffiently.large.clusters <- lw(cl.sizes >= RE.MIN.CLUSTER.SIZE)
  suffiently.large.clusters.exist <- no.suffiently.large.clusters>0
  if(suffiently.large.clusters.exist){
    cl.indices <- which(cl.sizes >= RE.MIN.CLUSTER.SIZE)
    
    # then, for each cluster
    ci=3
    
    # a) make an alignment
    this.cluster.index <- cl.indices[ci]
    this.cluster.name <- sprintf("clust%.*d", 2, ci)
    this.cluster.members <- names(which(cl==this.cluster.index))
    this.cluster.aln <- extract.aln.subset(aln, this.cluster.members)
    this.cluster.aln.outfile <- file.path(this.ecf.dir, sprintf("%s.aln", this.cluster.name))
    writeXStringSet(this.cluster.aln, filepath = this.cluster.aln.outfile)
    
    # b) calculate similarities between full sequences
    member.pairs <- as.data.table(t(combn(sort(this.cluster.members), 2)) )
    colnames(member.pairs) <- c("seq1", "seq2")
    member.pairs.rvr <- reprseq.vs.reprseq.blast[
      member.pairs,
      on = .(qname = seq1, sname = seq2),
      nomatch = 0L
    ]
    member.pairs.rvr[, max.cov := pmax(qcov.full, scov.full)]
    member.pairs.rvr[, sim.score := max.cov*pident/100]
    member.pairs.rvr[, .(pident, qcov.full, scov.full, max.cov, sim.score)]
    
    
    
    
    
    
  }
  
  
  largest.cluster.index <- which.max(table(cl))
  largest.cluster.size <- lw(cl == largest.cluster.index)
  largest.cl <- names(which(cl==largest.cluster.index))
  largest.cl.aln <- aln[which(names(aln) %in% largest.cl)]
  keep_cols <- colSums(as.matrix(largest.cl.aln) != "-") > 0
  largest.cl.aln.m <- as.matrix(largest.cl.aln)[,keep_cols]
  largest.cl.aln.trimmed <- AAStringSet(apply(largest.cl.aln.m, 1, paste0, collapse = ""))
  
  

  # # Compute distance = 1 - PID
  # # For big N this is heavy; see note below.
  # pid <- pid(aln, type = "PID1") / 100     # matrix in [0,1]
  # dist_mat <- as.dist(1 - pid)D
  # <- DistanceMatrix(
  # aln,
  # type = "AA",
  # correction = "none"   # raw p-distance
# )
  
  hc <- hclust(dist_mat, method = "complete")
  
  # 50% identity => distance <= 0.5
  cl <- cutree(hc, h = 0.5)
  
  
  ecf.tree.filename <- file.path(ecf.trees.dir, sprintf("%s-tree.png", ecf.name))
  if(!file.exists(ecf.tree.filename)){
    # Load phylogenetic tree
    tree_file <- file.path(ecf.subseqs.dir, sprintf("%s.aln.treefile", ecf.name))
    tree.unrooted <- read.tree(tree_file)
    tree <- midpoint(tree.unrooted)
    
    # Load ECOD domain data
    ecod_file <- file.path(ecf.subseqs.dir, sprintf("%s-ecod.txt", ecf.name))
    ecod_data <- read.csv(ecod_file, stringsAsFactors = FALSE)
    
    # Load PFAM domain data
    pfam_file <- file.path(ecf.subseqs.dir, sprintf("%s-pfam.txt", ecf.name))
    pfam_data <- read.csv(pfam_file, stringsAsFactors = FALSE)
    
    # Load conserved fragments and protein lengths
    frag_file <- file.path(ecf.subseqs.dir, sprintf("%s-frag.txt", ecf.name))
    fragment_data <- read.csv(frag_file, stringsAsFactors = FALSE)
    
    # Convert tree into a tibble for alignment
    tree_data <- tibble(label = tree$tip.label)
    
    # Ensure domain & fragment data align with tree
    pfam_data$y <- match(pfam_data$name, tree$tip.label)
    ecod_data$y <- match(ecod_data$name, tree$tip.label)
    fragment_data$y <- match(fragment_data$name, tree$tip.label)
    
    # Generate distinct colors for all unique ECOD and PFAM IDs
    set.seed(42)  # Ensures consistent color assignment
    unique_ecod_ids <- unique(ecod_data$id)
    unique_pfam_ids <- unique(pfam_data$id)
    
    # Assign colors dynamically for ECOD and PFAM **IDs**
    ecod_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique_ecod_ids)), unique_ecod_ids)
    pfam_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique_pfam_ids)), unique_pfam_ids)
    
    # Create Phylogenetic Tree Plot
    tree_plot <- ggtree(tree) +
      theme_tree2()
    
    # **PFAM domain plot with protein length as a grey line**
    pfam_plot <- ggplot(pfam_data, aes(xmin = start, xmax = end, y = name, fill = id)) +
      # Grey line for full protein length
      geom_segment(data = fragment_data, aes(x = 0, xend = length, y = name, yend = name), 
                   color = "grey40", size = 1, inherit.aes = FALSE) +
      # PFAM domains as arrows
      geom_gene_arrow() +  
      # Conserved fragment as a red segment
      geom_segment(data = fragment_data, 
                   aes(x = start, xend = end, y = name, yend = name), 
                   color = "red", size = 2, alpha = 0.8, inherit.aes = FALSE) +
      scale_fill_manual(values = pfam_colors) + 
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      labs(x = "PFAM Domains", fill = "PFAM ID")
    
    # **ECOD domain plot with protein length as a grey line**
    ecod_plot <- ggplot(ecod_data, aes(xmin = start, xmax = end, y = name, fill = id)) +
      # Grey line for full protein length
      geom_segment(data = fragment_data, aes(x = 0, xend = length, y = name, yend = name), 
                   color = "grey40", size = 1, inherit.aes = FALSE) +
      # ECOD domains as arrows
      geom_gene_arrow() +  
      # Conserved fragment as a red segment
      geom_segment(data = fragment_data, 
                   aes(x = start, xend = end, y = name, yend = name), 
                   color = "red", size = 2, alpha = 0.8, inherit.aes = FALSE) +
      scale_fill_manual(values = ecod_colors) + 
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      labs(x = "ECOD Domains", fill = "ECOD ID")
    
    # Combine tree, PFAM, and ECOD plots
    final_plot <- tree_plot + pfam_plot + ecod_plot + plot_layout(widths = c(1, 1.5, 1.5))
    ggsave(file = ecf.tree.filename, final_plot, width = 15, height = 9)  
  }
  
}
