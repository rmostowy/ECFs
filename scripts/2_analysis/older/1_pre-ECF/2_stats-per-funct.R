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
##### 1. Load data
######################################
reprseqs2rm <- readLines(reprseq.removed.filename)
qecf.stats <- fread(qecf.stats.filename)
qecf.stats[, qecf.name := sprintf("%s_%s", qecf.stats$qname, qecf.stats$cluster)]
qecf.stats <- qecf.stats[order(qname)]
qecf.stats[, V1 := NULL]
qecf.stats[avg_hit_coverage < 0, avg_hit_coverage := 0]
qecf.stats.filtered <- qecf.stats[qecf_span >= 20 & nb_hits >= 5 & avg_hit_coverage > 0 & avg_hits_overlap > 0]
qecf.stats.filtered <- qecf.stats.filtered[!qname %in% reprseqs2rm]

# reprseq.domains.ecod <- fread(reprseq.domains.ecod.output)
# reprseq.ecod <- reprseq.domains.ecod %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

# reprseq.domains.pfam <- fread(reprseq.domains.pfam.output)
# reprseq.pfam <- reprseq.domains.pfam %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.functions.all <- fread(reprseq.functions.all.output)
setkey(reprseq.functions.all, annotation)
top.annotations <- reprseq.functions.all[reprseq.functions.all$top.annotation] %>% select(funct, class, annotation) %>% distinct()

reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)
setorder(reprseq.ecod, qname, sname, -bitscore, eval)
reprseq.ecod <- reprseq.ecod[ , .SD[1], by = .(qname, sname)]

ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

reprseq.pfam.all <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.pfam.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)
setorder(reprseq.pfam, qname, sname, -bitscore, eval)
reprseq.pfam <- reprseq.pfam[ , .SD[1], by = .(qname, sname)]

pfam.domains <- fread(pfam.domains.output)
setkey(pfam.domains, accession)

######################################
##### 2. Calculate stats for functions
######################################

stats.per.funct.list <- lapply(1:nrow(top.annotations), function(index){
# stats.per.funct.list <- lapply(21:25, function(index){
  out <- NULL
  this.annotation <- top.annotations$annotation[index]
  this.funct <- top.annotations$funct[index]
  this.class <- top.annotations$class[index]
  this.annotation.reprseqs <- reprseq.functions.all[.(this.annotation)]$qname
  this.annotation.qecf.stats <- qecf.stats.filtered[qecf.stats.filtered$qname %in% this.annotation.reprseqs]
  no.qecfs.per.annot <- nrow(this.annotation.qecf.stats)
  if(no.qecfs.per.annot>0){
  out <- data.table(qecf.name = this.annotation.qecf.stats$qecf.name, aho = this.annotation.qecf.stats$avg_hits_overlap, 
             ac = this.annotation.qecf.stats$avg_hit_coverage, funct = this.funct, 
             class = this.class, annotation = this.annotation)
  }
  return(out)
})
stats.per.funct <- rbindlist(stats.per.funct.list)
stats.per.funct.vis <- stats.per.funct %>% select(funct, class, aho, ac)

######################################
##### 3. Plot functions
######################################

# Reordering 'funct' based on 'class'
stats.per.funct.vis[, funct := factor(funct, levels = unique(funct[order(class)]))]
# stats.per.funct.vis[, class := factor(class, levels = phrogs.class.order)]


# Plotting Aho and Ac values with color distinction for classes
p.aho <- ggplot(stats.per.funct.vis, aes(x = funct, y = aho, fill = class)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Aho Values per Funct, Ordered and Coloured by Class",
       x = "Function (sorted by Class)", y = "AHO Value") +
  scale_fill_manual(values = setNames(phrogs.col, phrogs.class))

p.ac <- ggplot(stats.per.funct.vis, aes(x = funct, y = ac, fill = class)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Aho Values per Funct, Ordered and Coloured by Class",
       x = "Function (sorted by Class)", y = "AC Value") +
  scale_fill_manual(values = setNames(phrogs.col, phrogs.class))

p.aho.filename <- "plot2a_funct-aho.pdf"
p.ac.filename <- "plot2b_funct-ac.pdf"

ggsave(file = p.aho.filename, p.aho, width = 20, height = 7)
ggsave(file = p.ac.filename, p.ac, width = 20, height = 7)

# Create a color vector for axis text
selected_functions <- c("endolysin", "minor tail", "tail fiber")
stats.per.funct.vis$a_color <- ifelse(stats.per.funct.vis$funct %in% selected_functions, "red", "black")

######################################
##### 4. Get low AHO/AC functions
######################################

mean.stats.per.funct <- stats.per.funct[, .(
  mean.aho = mean(aho, na.rm = T),
  mean.ac = mean(ac, na.rm = T)
  ), by = funct]

my.probs <- c(0.05, 0.1, 0.2, 0.25)
lowebound.aho <- quantile(stats.per.funct$aho, probs = my.probs, na.rm = T)["25%"]
lowebound.ac <- quantile(stats.per.funct$ac, probs = my.probs, na.rm = T)["25%"]

mean.stats.per.funct$low.aho <- FALSE
mean.stats.per.funct$low.ac <- FALSE
mean.stats.per.funct$low.aho[mean.stats.per.funct$mean.aho<lowebound.aho] <- TRUE
mean.stats.per.funct$low.ac[mean.stats.per.funct$mean.ac<lowebound.ac] <- TRUE

mean.stats.per.funct[mean.stats.per.funct$low.aho]
mean.stats.per.funct[mean.stats.per.funct$low.ac]


######################################
##### 5. Calculate stats for domains
######################################

# 1. Merge all fragments with ECODs by qname
merged.qecf.ecod <- qecf.stats.filtered[reprseq.ecod, on = .(qname), allow.cartesian = T]
merged.qecf.ecod <- merged.qecf.ecod[!is.na(qecf.name)]
merged.qecf.ecod[, overlap := get.fragment.domain.overlap(qecf_start, qecf_end, qstart, qend), by = .(qecf.name, sname)]
merged.qecf.ecod.filtered <- merged.qecf.ecod[order(-overlap, -bitscore), .SD[1], by = qecf.name]
merged.qecf.ecod.filtered <- merged.qecf.ecod.filtered[merged.qecf.ecod.filtered$overlap >= MIN.OVERLAP.THRESHOLD]
qecf.stats.per.ecod <- merged.qecf.ecod.filtered %>% select(qecf.name, aho = avg_hits_overlap, ac = avg_hit_coverage, sname)

qecf.stats.per.ecod$t.id <- ecod.domains[.(qecf.stats.per.ecod$sname)]$t_id
qecf.stats.per.ecod$t.name <- ecod.domains[.(qecf.stats.per.ecod$sname)]$t_name
qecf.stats.per.ecod$arch.class <- ecod.domains[.(qecf.stats.per.ecod$sname)]$arch.class

# 2. Merge all fragments with Pfam by qname
merged.qecf.pfam <- qecf.stats.filtered[reprseq.pfam, on = .(qname), allow.cartesian = T]
merged.qecf.pfam <- merged.qecf.pfam[!is.na(qecf.name)]
merged.qecf.pfam[, overlap := get.fragment.domain.overlap(qecf_start, qecf_end, qstart, qend), by = .(qecf.name, sname)]
merged.qecf.pfam.filtered <- merged.qecf.pfam[order(-overlap, -bitscore), .SD[1], by = qecf.name]
merged.qecf.pfam.filtered <- merged.qecf.pfam.filtered[merged.qecf.pfam.filtered$overlap >= MIN.OVERLAP.THRESHOLD]
qecf.stats.per.pfam <- merged.qecf.pfam.filtered %>% select(qecf.name, aho = avg_hits_overlap, ac = avg_hit_coverage, accession)

qecf.stats.per.pfam$annot <- pfam.domains[.(qecf.stats.per.pfam$accession)]$description

######################################
##### 5. Visualise domains
######################################

#### ECOD
min.number.of.qecfs.per.topology <- 100
no.qecfs.per.tid <- qecf.stats.per.ecod[, .N, by = t.id]
frequent.tids <- no.qecfs.per.tid[no.qecfs.per.tid$N >= min.number.of.qecfs.per.topology]

qecf.stats.per.ecod.freq <- qecf.stats.per.ecod[qecf.stats.per.ecod$t.id %in% frequent.tids$t.id]
stats.per.ecod.vis <- qecf.stats.per.ecod.freq %>% select(t.name, t.id, aho, ac)

# Plotting Aho and Ac values with color distinction for classes
p.de.aho <- ggplot(stats.per.ecod.vis, aes(x = t.name, y = aho)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Aho Values per ECOD topology",
       x = "T-id", y = "AHO Value")
p.de.ac <- ggplot(stats.per.ecod.vis, aes(x = t.name, y = ac)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "AC Values per ECOD topology",
       x = "T-id", y = "AC Value")

p.de.aho.filename <- "plot3a_ecod-aho.pdf"
p.de.ac.filename <- "plot3b_ecod-ac.pdf"

ggsave(file = p.de.aho.filename, p.de.aho, width = 20, height = 7)
ggsave(file = p.de.ac.filename, p.de.ac, width = 20, height = 7)

# # Plot for architecutre class
# ggplot(qecf.stats.per.ecod, aes(x = arch.class, y = aho)) + geom_boxplot()

#### PFAM
min.number.of.qecfs.per.pfam <- 100
no.qecfs.per.pfam <- qecf.stats.per.pfam[, .N, by = annot]
frequent.pfams <- no.qecfs.per.pfam[no.qecfs.per.pfam$N >= min.number.of.qecfs.per.pfam]

qecf.stats.per.pfam.freq <- qecf.stats.per.pfam[qecf.stats.per.pfam$annot %in% frequent.pfams$annot]
stats.per.pfam.vis <- qecf.stats.per.pfam.freq %>% select(annot, aho, ac)

# Plotting Aho and Ac values with color distinction for classes
p.dp.aho <- ggplot(stats.per.pfam.vis, aes(x = annot, y = aho)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "AHO Values per Pfam annotation",
       x = "Pfam annotation", y = "AHO Value")
p.dp.ac <- ggplot(stats.per.pfam.vis, aes(x = annot, y = ac)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "AC Values per Pfam annotation",
       x = "Pfam annotation", y = "AC Value")

p.dp.aho.filename <- "plot4a_pfam-aho.pdf"
p.dp.ac.filename <- "plot4b_pfam-ac.pdf"

ggsave(file = p.dp.aho.filename, p.dp.aho, width = 20, height = 7)
ggsave(file = p.dp.ac.filename, p.dp.ac, width = 20, height = 7)

