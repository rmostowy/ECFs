suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

# Load domain and ECF dbs
reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.pfam.all <- fread(reprseq.domains.pfam.output)

reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

reprseq.ecod <- reprseq.ecod.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

reprseq.pfam <- reprseq.pfam.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

# Let us assume that two reprseqs share a domain if we find it in both of them at MIN.SCOV.FOR.DOMAIN. 

# Find all pairs of reprseq that share a Pfam domain
reprseq.pfam.minimal <- reprseq.pfam %>% select(qname, accession) %>% distinct()
pairs.that.share.pfam <- merge(reprseq.pfam.minimal, reprseq.pfam.minimal, by = "accession", allow.cartesian = T)
pairs.that.share.pfam <- pairs.that.share.pfam[qname.x != qname.y]
pairs.that.share.pfam <- pairs.that.share.pfam[qname.x < qname.y]
pairs.that.share.pfam[, pair := paste(qname.x, qname.y, sep = "-")]

# Now find all pairs of reprseq that share an ECOD domain
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)
reprseq.ecod$t_id <- ecod.domains[.(reprseq.ecod$sname)]$t_id

reprseq.ecod.minimal <- reprseq.ecod %>% select(qname, t_id) %>% distinct()
pairs.that.share.ecod <- merge(reprseq.ecod.minimal, reprseq.ecod.minimal, by = "t_id", allow.cartesian = T)
pairs.that.share.ecod <- pairs.that.share.ecod[qname.x != qname.y]
pairs.that.share.ecod <- pairs.that.share.ecod[qname.x < qname.y]
pairs.that.share.ecod[, pair := paste(qname.x, qname.y, sep = "-")]

# SO which pairs share any domain?
all.pairs <- union(unique(pairs.that.share.pfam$pair), unique(pairs.that.share.ecod$pair))
no.all.pairs <- length(all.pairs)

lw(all.pairs %in% pairs.that.share.pfam$pair)/no.all.pairs
lw(all.pairs %in% pairs.that.share.ecod$pair)/no.all.pairs


# Now with ECFs
reprseq.ecf.minimal <- reprseq.ecfs %>% select(qname, sname) %>% distinct()
pairs.that.share.ecf <- merge(reprseq.ecf.minimal, reprseq.ecf.minimal, by = "sname", allow.cartesian = T)
pairs.that.share.ecf <- pairs.that.share.ecf[qname.x != qname.y]
pairs.that.share.ecf <- pairs.that.share.ecf[qname.x < qname.y]
pairs.that.share.ecf[, pair := paste(qname.x, qname.y, sep = "-")]


lw(all.pairs %in% pairs.that.share.ecf$pair)








