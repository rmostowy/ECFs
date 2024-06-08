### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load ECF table
ecfs.vs.qecfs <- fread(ecf.vs.qecf.output)
ecfs.vs.qecfs$reprseq <- substr(ecfs.vs.qecfs$qecf, 1, regexpr("_", ecfs.vs.qecfs$qecf)-1)
setkey(ecfs.vs.qecfs, ecf)
all.ecfs <- unique(ecfs.vs.qecfs$ecf)

# Extract ECF lengths from ECF-ECF comparison
ecfs.vs.ecfs.all <- fread(ecfs.vs.ecfs.output)
ecf.lengths <- rbind(
  ecfs.vs.ecfs.all %>% select(name = qname, length = qlength),
  ecfs.vs.ecfs.all %>% select(name = sname, length = slength)
)
ecf.lengths <- unique(ecf.lengths)
setorder(ecf.lengths, name)
setkey(ecf.lengths, name)
outlier.ecf.lengths <- quantile(ecf.lengths$length, probs = c(0.75, 0.9, 0.95, 0.975, 0.99))
outlier99.length <- outlier.ecf.lengths[3]
outlier.ecfs <- ecf.lengths$name[ecf.lengths$length>outlier99.length]

# Load ECF classification table
ecf.classification.table <- fread(ecf.classification.table.output)
setkey(ecf.classification.table, ecf)

# Create networks
network.thresholds <- data.table(qcov = c(0.5, 0.8, .9, 0), scov = c(0.5, 0.8, 0.9, 0))
no.thresholds <- nrow(network.thresholds)

ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.lengths[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.lengths[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]
# filter outliers
ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[(!ecfs.vs.ecfs.cov$query %in% outlier.ecfs) & (!ecfs.vs.ecfs.cov$subject %in% outlier.ecfs)]

# p <- ggplot(ecfs.vs.ecfs.cov, aes(x = qlength, y = slength)) + geom_bin2d(bins = 500) +
#   geom_abline(slope = cutoff.val, intercept = 0, color = "red") +   # Line for qlength = 10 * slength
#   geom_abline(slope = 1/cutoff.val, intercept = 0, color = "red")    # Line for slength = 10 * qlength

# filter ECF pairs with an extreme length difference
cutoff.val <- 5
ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[!(ecfs.vs.ecfs.cov$qlength*cutoff.val < ecfs.vs.ecfs.cov$slength | ecfs.vs.ecfs.cov$qlength/cutoff.val > ecfs.vs.ecfs.cov$slength)]

for(index in 1:no.thresholds){
  this.qcov.thr <- network.thresholds$qcov[index]
  this.scov.thr <- network.thresholds$scov[index]
  this.ecf.network <- ecfs.vs.ecfs.cov %>% filter(qcov.min >= this.qcov.thr & scov.min >= this.scov.thr)
  this.ecf.network <- this.ecf.network %>% select(query, subject)
  singleton.network <- data.table(query = all.ecfs, subject = all.ecfs)
  this.ecf.network <- rbind(this.ecf.network, singleton.network)
  this.ecf.network.name <- sprintf("ecf-network-q%s-s%s.txt", num2s(this.qcov.thr), num2s(this.scov.thr))
  this.ecf.network.filename <- file.path(ecf.network.dir, this.ecf.network.name)
  fwrite(this.ecf.network, file = this.ecf.network.filename, sep = "\t")

}

network.thresholds.max.min <- data.table(
  max = c(.9,
          .9,
          .9,
          .95,
          .99, 
          .99), 
  min = c(0,
          .1,
          .3,
          0,
          0, 
          .1))
no.thresholds <- nrow(network.thresholds.max.min)

for(index in 1:no.thresholds){
  this.min.thr <- network.thresholds.max.min$min[index]
  this.max.thr <- network.thresholds.max.min$max[index]
  this.ecf.network <- copy(ecfs.vs.ecfs.cov)
  this.ecf.network <- this.ecf.network %>% filter(cov.min >= this.min.thr & cov.max >= this.max.thr)
  this.ecf.network <- this.ecf.network %>% select(query, subject)
  singleton.network <- data.table(query = all.ecfs, subject = all.ecfs)
  this.ecf.network <- rbind(this.ecf.network, singleton.network)
  this.ecf.network.name <- sprintf("ecf-network-min%s-max%s.txt", num2s(this.min.thr), num2s(this.max.thr))
  this.ecf.network.filename <- file.path(ecf.network.dir, this.ecf.network.name)
  fwrite(this.ecf.network, file = this.ecf.network.filename, sep = "\t")
  
}

replicative.classes <- c("DNA, RNA and nucleotide metabolism", "integration and excision", "transcription regulation")

# sanity check
t <- ecfs.vs.ecfs.cov[ecfs.vs.ecfs.cov$cov.max>=0 & ecfs.vs.ecfs.cov$cov.min>=0]
t$q.funct <- ecf.classification.table[.(t$query)]$funct 
t$q.class <- ecf.classification.table[.(t$query)]$class
t$s.funct <- ecf.classification.table[.(t$subject)]$funct 
t$s.class <- ecf.classification.table[.(t$subject)]$class

weird.edges <- (t$q.class %in% "tail" & t$s.class %in% replicative.classes) |
  (t$s.class %in% "tail" & t$q.class %in% replicative.classes)
t.weird <- t[weird.edges]
t.weird <- t.weird[t.weird$cov.max>0.9]
# 
# View(t.weird[t.weird$s.funct == "tail fiber" | t.weird$q.funct == "tail fiber"])
# 

t[(t$q.funct == 'tail assembly' & t$s.funct == "tail fiber") | 
  t$s.funct == 'tail assembly' & t$q.funct == "tail fiber"]


