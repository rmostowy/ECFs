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
qecf.stats <- fread("reprseq-all-stat-mine.txt")
qecf.stats[, qecf.name := sprintf("%s_%s", qecf.stats$qname, qecf.stats$cluster)]

reprseq.domains.ecod <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.domains.ecod %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.domains.pfam <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.domains.pfam %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.functions.all <- fread(reprseq.functions.all.output)
# reprseq.functions.all <- reprseq.functions.all[reprseq.functions.all$qname %in% all.reprseq]
setkey(reprseq.functions.all, annotation)

top.annotations <- reprseq.functions.all[reprseq.functions.all$top.annotation] %>% select(funct, class, annotation) %>% distinct()

######################################
##### 2. Calculate overlap stats
######################################

stats.per.funct.list <- lapply(1:nrow(top.annotations), function(index){
# stats.per.funct.list <- lapply(21:25, function(index){
  out <- NULL
  this.annotation <- top.annotations$annotation[index]
  this.funct <- top.annotations$funct[index]
  this.class <- top.annotations$class[index]
  this.annotation.reprseqs <- reprseq.functions.all[.(this.annotation)]$qname
  this.annotation.qecf.stats <- qecf.stats[qecf.stats$qname %in% this.annotation.reprseqs]
  no.qecfs.per.annot <- nrow(this.annotation.qecf.stats)
  if(no.qecfs.per.annot>0){
  out <- data.table(qecf.name = this.annotation.qecf.stats$qecf.name, aho = this.annotation.qecf.stats$average.frag.overlap, 
             ac = this.annotation.qecf.stats$average.coverage, funct = this.funct, 
             class = this.class, annotation = this.annotation)
  }
  return(out)
})
stats.per.funct <- rbindlist(stats.per.funct.list)
stats.per.funct.vis <- stats.per.funct %>% select(funct, class, aho, ac)

######################################
##### 3. Plot
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

p.aho.filename <- "plot2_aho.pdf"
p.ac.filename <- "plot3_ac.pdf"

ggsave(file = p.aho.filename, p.aho, width = 20, height = 7)
ggsave(file = p.ac.filename, p.ac, width = 20, height = 7)


# Create a color vector for axis text
selected_functions <- c("endolysin", "minor tail", "tail fiber")
stats.per.funct.vis$a_color <- ifelse(stats.per.funct.vis$funct %in% selected_functions, "red", "black")

