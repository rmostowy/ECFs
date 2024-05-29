root <- rprojroot::find_rstudio_root_file()

#####################
## Define wrappers
#####################

safe_source <- function(filepath) {
  if (file.exists(filepath)) {
    source(filepath)
  } else {
    warning(sprintf("The file '%s' does not exist. Please ensure the repository is set up correctly.", filepath))
  }
}

#####################
## Define main paths
#####################

main.path <- "/Users/rmostowy"
main.data.dir <- file.path(main.path, "Data", "ECFs")
main.data.tmp.dir <- file.path(main.data.dir, "tmp-ecfs")

functions.repo.path <- file.path(main.path, "Github", "proj-rmostowy", "functions")

mgg.db <- file.path(main.path, "MGG Dropbox/Databases")
phrogs.db.path <- file.path(mgg.db, "Phrogs")

mgg.projects <- file.path(main.path, "MGG Dropbox/Projects")
minimal.input.data.path <- file.path(mgg.projects, "divRBP/minimal-input-data")

ecf.explorer.path <- file.path(mgg.projects, "ECFs/ecf-explorer")
ecf.raw.data.path <- file.path(ecf.explorer.path, "ecf-explorer_2024-05-23", "output")

##############################
### Load functions
##############################
functions.basic.filename <- file.path(functions.repo.path, "basic.R")
safe_source(functions.basic.filename)


#############################
### Define PHROGs data
#############################
phrogs.annot.filename <- file.path(phrogs.db.path, "phrog_annot_v4.tsv")
phrogs.table.filename <- file.path(phrogs.db.path, "phrogs_table_almostfinal_plusGO_wNA_utf8.tsv")

#############################
### Define Smug et al. data
#############################
all.cds.filename <- file.path(minimal.input.data.path, "coding-seqs/cds-aa.fa") 
reprseq.families.raw.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "families", "dataset-full", "mcl", "repr-hits-pairwise-prob95-cov80-mcl-I20.out")
name.table.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "representative", "name-table-full.txt")
reprseq.lengths.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "representative", "repr-seqs-lengths.txt")
phrog.hhpred.hits.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "functional", "hhblits-phrogs.txt")
manual.phrog.mapping.filename <- file.path(minimal.input.data.path, "other", "v3_phrogs-table-rafal-3_12.xlsx")
antidefence.descr.filename <- file.path(minimal.input.data.path, "other", "media-4-simplified_2023_04.xlsx")
antidefence.hhpred.hits.filename <- file.path(minimal.input.data.path, "other", "table-hhr-ad-rseq.txt")

####################
### Define ECF data
####################

ecf.name.table.filename = file.path(ecf.raw.data.path,  "scan-ecf-dict.txt")
ecf.hits.filename = file.path(ecf.raw.data.path, "table-scan-reprseq_ecf_db.txt")
ecf.hits.vs.ecod.filename = file.path(ecf.raw.data.path, "table-scan-ecf_ecod.txt")
ecf.hits.vs.pfam.filename = file.path(ecf.raw.data.path, "table-scan-ecf_pfam.txt")
ecf.hits.vs.ecf.filename = file.path(ecf.raw.data.path, "table-scan-ecf_vs_ecf.txt")
reprseq.vs.reprseq.filename = file.path(ecf.raw.data.path, "reprseq-all-clust.txt")
ecod.domain.hits.filename = file.path(ecf.raw.data.path, "table-scan-reprseq_ecod.txt")
qecf.stats.filename = file.path(ecf.raw.data.path, "reprseq-all-stat.txt")


#################
## Labels
#################

phrogs.class <- c("integration and excision", "head and packaging", 
                  "transcription regulation", "connector", 
                  "tail", "lysis", 
                  "DNA, RNA and nucleotide metabolism", "other",
                  "moron, auxiliary metabolic gene and host takeover")
phrogs.col <- c("#FF7F00", "#C7E9C0", "#984EA3", "#41AB5D", "#006D2C", "#E41A1C", "#377EB8", "#999999", "#F781BF")

phrogs.class.order <- c("head and packaging", "connector", "tail", 
                        "DNA, RNA and nucleotide metabolism", "integration and excision",
                        "moron, auxiliary metabolic gene and host takeover", "transcription regulation",
                        "lysis", "other")

phrog.class.labels <- c(
  'integration and excision' = "integration\nand excision",
  'head and packaging' = "head and\npackaging",
  "connector" = "connector",
  "tail" = "tail",
  "lysis" = "lysis",
  'transcription regulation' = "transcription\nregulation",
  'DNA, RNA and nucleotide metabolism' = "DNA, RNA and\nnucleotide\nmetabolism",
  "other" = "other",
  'moron, auxiliary metabolic gene and host takeover' = "moron, auxiliary\nmetabolic gene\nand host takeover"
)

#####################
## Check file status
#####################

files.to.check <- ls(pattern = "filename$")
files.to.check.full.path <- sapply(files.to.check, function(k) get(k, envir = .GlobalEnv))
file.status.list <- sapply(files.to.check.full.path, check.file)

