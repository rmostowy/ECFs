### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# merge original PHROG tables into a single table
phrogs.functions <- fread(phrogs.annot.filename)
phrogs.table.full <- fread(phrogs.table.filename, header = F)

phrogs.custom <- as.data.table(read_excel(manual.phrog.mapping.filename, 1))
phrogs.merged <- merge.data.table(phrogs.functions, phrogs.table.full, by.x = "phrog", by.y = "V1")
phrogs.functions <- data.table(phrog = phrogs.merged$phrog,
                               funct = phrogs.merged$annot,
                               class = phrogs.merged$category,
                               annotation.orig = sprintf("%s|%s", phrogs.merged$annot, phrogs.merged$category),
                               no.seqs = phrogs.merged$V5)
phrogs.functions <- phrogs.functions[phrogs.functions$no.seqs>1]
phrogs.functions <- phrogs.functions[!is.na(phrogs.functions$funct)]
phrogs.functions <- phrogs.functions[phrogs.functions$funct != "unknown function"]

# substitute annotations into new ones
phrogs.custom$annot.orig <- sprintf("%s|%s", phrogs.custom$funct.orig, phrogs.custom$class.orig)
phrogs.custom$annot.new <- sprintf("%s|%s", phrogs.custom$funct.new, phrogs.custom$class.new)
setkey(phrogs.custom, annot.orig)

phrogs.functions$annotation.new <- phrogs.custom[.(phrogs.functions$annotation.orig)]$annot.new
phrogs.functions$funct.new <- phrogs.custom[.(phrogs.functions$annotation.orig)]$funct.new
phrogs.functions$class.new <- phrogs.custom[.(phrogs.functions$annotation.orig)]$class.new

phrogs.functions.final <- data.table(phrog = phrogs.functions$phrog,
                                     funct = phrogs.functions$funct.new,
                                     class = phrogs.functions$class.new,
                                     annotation = phrogs.functions$annotation.new,
                                     no.seqs = phrogs.functions$no.seqs,
                                     annotation.orig = phrogs.functions$annotation.orig)

# calculate function count in phrogs
funct.count.seqs <- phrogs.functions.final[, sum(no.seqs), by = .(funct, class, annotation)]
funct.count.phrogs <- phrogs.functions.final[, .N, by = .(funct, class, annotation)]
funct.count <- merge.data.table(funct.count.seqs, funct.count.phrogs, by = c("funct", "class", "annotation"))

setnames(funct.count, old = c("V1", "N"), new = c("no.seqs", "no.phrogs"))
funct.count <- funct.count[order(funct.count$class, -funct.count$no.seqs, funct.count$funct)]
funct.count <- data.table(annot.index=1:nrow(funct.count), funct.count)

# fwrite(phrogs.functions.final, file = phrogs.functions.output)
# fwrite(funct.count, file = phrogs.count.output)
