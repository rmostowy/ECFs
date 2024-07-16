### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

######################################
##### 1. Define paths and parameters
######################################

min.qp.overlap <- 0.8 # minimum overlap between qECF and Pfam

stats.fn <- "reprseq-all-stat-mine.txt"
stats <- fread(stats.fn)

######################################
##### 2. Define functions
######################################

get.qp.overlap <- function(fragment_start, fragment_end, pfam_start, pfam_end){
  intersection_start <- max(fragment_start, pfam_start)
  intersection_end <- min(fragment_end, pfam_end)
  intersect_length <- max(0, intersection_end - intersection_start + 1)
  union_start <- min(fragment_start, pfam_start)
  union_end <- max(fragment_end, pfam_end)
  union_length <- union_end - union_start + 1
  overlap <- intersect_length/union_length
  return(overlap)
}

######################################
##### 3. Calculate overlap stats
######################################

pfam.hits <- fread(reprseq.domains.pfam.output)
pfam.hits.filtered <- pfam.hits %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)
setkey(pfam.hits.filtered, qname)

no.qecfs <- nrow(stats)

stats[, qecf.name := sprintf("%s_%s", stats$qname, stats$cluster)]

max.pfam.overlap <- sapply(1:no.qecfs, function(qecf.index){
  out <- NA
  qecf.name <- stats$qecf.name[qecf.index]
  qecf.qname <- stats$qname[qecf.index]
  qecf.start <- stats$qecf.start[qecf.index]
  qecf.end <- stats$qecf.end[qecf.index]
  query.pfams <- pfam.hits.filtered[.(qecf.qname)]
  no.pfams <- length(which(!is.na(query.pfams$qstart)))
  if(no.pfams>0){
    qecf.overlaps <- sapply(1:no.pfams, function(i) get.qp.overlap(qecf.start, qecf.end, query.pfams$qstart[i], query.pfams$qend[i]))
    max.pfam.overlap <- max(qecf.overlaps)
    out <- max.pfam.overlap
  }
  return(out)
})

stats$max.pfam.overlap <- max.pfam.overlap
stats$is.pfam.equiv <- FALSE
stats$is.pfam.equiv[which(max.pfam.overlap>=min.qp.overlap)] <- TRUE

######################################
##### 3. Visualise
######################################


stats.df <- as.data.frame(stats)
stats.df$is.pfam.equiv <- as.factor(stats.df$is.pfam.equiv)

# Assuming your data frame is called stats and has columns average.frag.overlap and is.pfam.equiv
p <- ggplot(stats.df, aes(x = is.pfam.equiv, y = average.frag.overlap, fill = is.pfam.equiv)) +
  geom_boxplot(alpha = 0.5) +  # Boxplot with semi-transparent fill for better visibility of points
  geom_jitter(width = 0.2, alpha = 0.3, color = "black", size = .5) +  # Jittered points overlay
  scale_fill_brewer(palette = "Set1") +  # Optional: Use a Color Brewer palette for fills
  labs(title = "Distribution of Average Fragment Overlap",
       x = "Is Pfam Equivalent",
       y = "Average Fragment Overlap") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide the legend if it's not needed

p.filename <- "plot1_pfam-equiv-aho.png"
ggsave(p.filename, p, width = 6, height = 5)

######################################
##### 4. Get value cut-offs
######################################

lowerbound.ac <- quantile(stats$average.coverage[stats$is.pfam.equiv], probs = c(0.05, 0.1, 0.2, 0.25))
lowerbound.aho <- quantile(stats$average.frag.overlap[stats$is.pfam.equiv], probs = c(0.05, 0.1, 0.2, 0.25))

cutoff.ac <- as.numeric(lowerbound.ac["10%"])
cutoff.aho <- as.numeric(lowerbound.aho["10%"])

# stats.qcfs <- stats[stats$]
# stats[stats$average.frag.overlap >= cutoff.aho]
stats.qcfs <- stats[stats$average.frag.overlap >= cutoff.aho & stats$average.coverage >= cutoff.ac]
nrow(stats.qcfs)



