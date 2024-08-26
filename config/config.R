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
mgg.db <- file.path(main.path, "MGG Dropbox/Databases")
mgg.projects <- file.path(main.path, "MGG Dropbox/Projects")

##############################
### Load functions
##############################
functions.repo.path <- file.path(main.path, "Github", "proj-rmostowy", "functions")
scripts.functions.path <- file.path(root, "scripts", "functions")

functions.basic.filename <- file.path(functions.repo.path, "basic.R")
safe_source(functions.basic.filename)

functions.custom.cov.filename <- file.path(functions.repo.path, "custom", "hhr-cov.R")
safe_source(functions.custom.cov.filename)

scripts.functions.filename <- file.path(scripts.functions.path, "helper.R")
safe_source(scripts.functions.filename)

#############################
### Define PHROGs data
#############################
phrogs.db.path <- file.path(mgg.db, "Phrogs")

phrogs.annot.filename <- file.path(phrogs.db.path, "phrog_annot_v4.tsv")
phrogs.table.filename <- file.path(phrogs.db.path, "phrogs_table_almostfinal_plusGO_wNA_utf8.tsv")
manual.phrog.mapping.filename <- file.path(phrogs.db.path, "custom", "v4_phrogs-table-2024.xlsx")

#############################
### Define ECOD data
#############################
ecod.db.path <- file.path(mgg.db, "hhpred", "ECOD", "ECOD_F70_20220613")
ecod.domains.filename <- file.path(ecod.db.path, "ecod.latest.F70.domains.txt")
ecod.domains.fasta.filename <- file.path(ecod.db.path, "ecod.latest.F70.fasta.txt")

ecod.headers <- c("uid", "ecod_domain_id","manual_rep","f_id","pdb","chain","pdb_range","seqid_range","unp_acc",
  "arch_name","x_name","h_name","t_name","f_name", "asm_status", "ligand")

#############################
### Define PFAM data
#############################
pfam.db.path <- file.path(mgg.db, "hhpred", "pfamA", "Pfam35")
pfam.domains.filename <- file.path(pfam.db.path, "pfamA.txt")
pfam.sql.filename <- file.path(pfam.db.path, "pfamA.sql")

#############################
### Define Smug et al. data
#############################
minimal.input.data.path <- file.path(mgg.projects, "divRBP/minimal-input-data")

all.cds.filename <- file.path(minimal.input.data.path, "coding-seqs/cds-aa.fa") 
reprseq.families.raw.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "families", "dataset-full", "mcl", "repr-hits-pairwise-prob95-cov80-mcl-I20.out")
name.table.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "representative", "name-table-full.txt")
reprseq.lengths.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "representative", "repr-seqs-lengths.txt")
phrog.hhpred.hits.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "functional", "hhblits-phrogs.txt")
antidefence.descr.filename <- file.path(minimal.input.data.path, "other", "media-4-simplified_2023_04.xlsx")
antidefence.hhpred.hits.filename <- file.path(minimal.input.data.path, "other", "table-hhr-ad-rseq.txt")
metadata.table.filename <- file.path(minimal.input.data.path, "metadata", "refseq_metadata_updated_v3.csv")

####################
### Define ECF data
####################
ecf.explorer.path <- file.path(mgg.projects, "ECFs/ecf-explorer")
ecf.raw.input.path <- file.path(ecf.explorer.path, "ecf-explorer_2024-05-23", "input")
ecf.raw.data.path <- file.path(ecf.explorer.path, "ecf-explorer_2024-05-23", "output")
# 
ecf.name.table.filename <- file.path(ecf.raw.data.path,  "scan-ecf-dict.txt")
ecf.hits.filename <- file.path(ecf.raw.data.path, "table-scan-reprseq_ecf_db.txt")
ecf.hits.vs.ecod.filename <- file.path(ecf.raw.data.path, "table-scan-ecf_ecod.txt")
ecf.hits.vs.pfam.filename <- file.path(ecf.raw.data.path, "table-scan-ecf_pfam.txt")
ecf.hits.vs.ecf.filename <- file.path(ecf.raw.data.path, "table-scan-ecf_vs_ecf.txt")
reprseq.vs.reprseq.filename <- file.path(ecf.raw.input.path, "table-hhr.txt")
reprseq.vs.reprseq.clust.filename <- file.path(ecf.raw.data.path, "reprseq-all-clust.txt")
ecod.domain.hits.filename <- file.path(ecf.raw.data.path, "table-scan-reprseq_ecod.txt")
pfam.domain.hits.filename <- file.path(ecf.raw.data.path, "table-scan-reprseq_pfam.txt")
qecf.stats.filename <- file.path(ecf.raw.data.path, "reprseq-all-stat-coverage.txt")
ecf.vs.qecf.filename <- file.path(ecf.raw.data.path, "ecf2qecfs.csv")

####################
### Define other data
####################

ictv.release.filename <- file.path(mgg.db, "ICTV", "VMR", "VMR_MSL38_v2.xlsx")

#######################
### Data directories
#######################
data.tmp.dir <- file.path(main.path, "Data", "tmp-ecfs")

main.data.dir <- file.path(main.path, "Dropbox", "Projects", "Data", "ECFs")
data.input.dir <- file.path(main.data.dir, "input")
data.output.dir <- file.path(main.data.dir, "output")

data.phrog.dir <- file.path(data.output.dir, "annotation", "phrogs")
data.domain.dir <- file.path(data.output.dir, "annotation", "domains")
data.ecfs.dir <- file.path(data.output.dir, "ecfs")

data.comparisons.dir <- file.path(data.output.dir, "comparisons")
ecf.network.dir <- file.path(data.comparisons.dir, "networks")

dirs2create <- c(data.tmp.dir,
                 main.data.dir, data.input.dir, data.output.dir,
                 data.phrog.dir, data.domain.dir, data.ecfs.dir,
                 data.comparisons.dir, ecf.network.dir)

make.dirs(dirs2create)

#######################
### Output files
#######################

# input data
reprseq.table.filename <- file.path(data.input.dir, "name-table.txt")
reprseq.seqs.filename <- file.path(data.input.dir, "reprseq-seq.fa")

# annotations/phrogs
phrogs.functions.output <- file.path(data.phrog.dir, "phrogs-table.txt")
phrogs.count.output <- file.path(data.phrog.dir, "phrogs-count.txt")
antidefence.table.output <- file.path(data.phrog.dir, "antidefence-table.txt")

reprseq.functions.all.output <- file.path(data.phrog.dir, "reprseq-annot-all.txt")
reprseq.functions.unique.output <- file.path(data.phrog.dir, "reprseq-annot-unique.txt")

# annotations/domains
ecod.domains.output <- file.path(data.domain.dir, "ecod-domains.txt")
pfam.domains.output <- file.path(data.domain.dir, "pfam-domains.txt")
reprseq.domains.ecod.output <- file.path(data.domain.dir, "reprseq-domains-ecod.txt") 
reprseq.domains.pfam.output <- file.path(data.domain.dir, "reprseq-domains-pfam.txt") 

# ecfs
ecf.vs.qecf.output <- file.path(data.ecfs.dir, "ecfs-vs-qecfs.txt") 

reprseq.ecfs.output <- file.path(data.ecfs.dir, "reprseq-ecfs.txt") 

ecfs.vs.ecfs.output <- file.path(data.ecfs.dir, "ecfs-vs-ecfs.txt")
ecfs.vs.ecod.output <- file.path(data.ecfs.dir, "ecfs-vs-ecod.txt")
ecfs.vs.pfam.output <- file.path(data.ecfs.dir, "ecfs-vs-pfam.txt")

ecfs.vs.ecfs.cov.output <- file.path(data.ecfs.dir, "ecfs-vs-ecfs-cov.txt")

ecf.family.table.output <- file.path(data.ecfs.dir, "ecfs-vs-EFs.txt")

# metadata
metadata.table.output <- file.path(data.input.dir, "metadata.tsv")

# comparisons
reprseq.ecf.vs.ecod.output <- file.path(data.comparisons.dir, "reprseq_ecf-vs-ecod.txt")
reprseq.ecf.vs.pfam.output <- file.path(data.comparisons.dir, "reprseq_ecf-vs-pfam.txt")
reprseq.ecf.vs.ecod.and.pfam.output <- file.path(data.comparisons.dir, "reprseq_ecf-vs-ecod-and-pfam.txt")

ecf.classification.table.output <- file.path(data.comparisons.dir, "ecf-classifications.txt")

# ecf-sharing pairs
ecf.sharing.pairs.output <- file.path(data.comparisons.dir, "ecf-sharing-pairs.txt")

#################
## Labels
#################

phrogs.class <- c("integration and excision", "head and packaging", 
                  "transcription regulation", "connector", 
                  "tail", "lysis", 
                  "DNA, RNA and nucleotide metabolism", "other",
                  "moron, auxiliary metabolic gene and host takeover",
                  "antidefence")
phrogs.col <- c("#FF7F00", "#C7E9C0", "#984EA3", "#41AB5D", "#006D2C", "#E41A1C", "#377EB8", "#999999", "#F781BF", "#111111")

phrogs.class.order <- c("head and packaging", "connector", "tail", 
                        "DNA, RNA and nucleotide metabolism", "integration and excision",
                        "moron, auxiliary metabolic gene and host takeover", "transcription regulation",
                        "lysis", "antidefence", "other")

phrog.class.labels <- c(
  'integration and excision' = "integration\nand excision",
  'head and packaging' = "head and\npackaging",
  "connector" = "connector",
  "tail" = "tail",
  "lysis" = "lysis",
  'transcription regulation' = "transcription\nregulation",
  'DNA, RNA and nucleotide metabolism' = "DNA, RNA and\nnucleotide\nmetabolism",
  "other" = "other",
  'moron, auxiliary metabolic gene and host takeover' = "moron, auxiliary\nmetabolic gene\nand host takeover",
  "antidefence" = "antidefence"
)

ecod.arch <- c("a/b barrels", "a/b three-layered sandwiches", "a+b complex topology",
               "a+b duplicates or obligate multimers", "a+b four layers", "a+b three layers",
               "a+b two layers", "alpha arrays", "alpha bundles", 
               "alpha complex topology", "alpha duplicates or obligate multimers", "alpha superhelices",
               "beta barrels", "beta complex topology", "beta duplicates or obligate multimers",
               "beta meanders", "beta sandwiches", "extended segments", 
               "few secondary structure elements", "mixed a+b and a/b")

ecod.arch.class <- c("alpha/beta", "alpha/beta", "alpha+beta", "alpha+beta",
                     "alpha+beta", "alpha+beta", "alpha+beta", "all-alpha",
                     "all-alpha", "all-alpha", "all-alpha", "all-alpha",
                     "all-beta", "all-beta", "all-beta", "all-beta",
                     "all-beta", "others", "others", "mixed alpha+beta and alpha/beta")

#####################
## Check file status
#####################

files.to.check <- ls(pattern = "filename$")
files.to.check.full.path <- sapply(files.to.check, function(k) get(k, envir = .GlobalEnv))
file.status.list <- sapply(files.to.check.full.path, check.file)

#####################
## PARAMETERS
#####################

# annotation parameters used also in Smug et al
ACCEPTED.NUM.UNKNOWN.AA <- 10
MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION <- 30
MIN.NUM.PHROG.SEQ.PER.ANNOT <- 500
MIN.NUM.PROT.PER.ANNOT <- 100 # Bogna had MIN.NUM.PROT
# MAIN.COVS.FOR.ANNOTATION = c(0.1, 0.3, 0.8)
DEFAULT.MINIMUM.COV.FOR.ANNOTATION <- 60
DEFAULT.MINIMUM.PROB.FOR.ANNOTATION <- 95
DEFAULT.MAXIMUM.EVAL.FOR.ANNOTATION <- 1e05
UNSPECIFIC.ANNOTATIONS = c("tail", "structural protein", "uknown")

# parameters related to domain detection
MIN.DOMAIN.LENGTH <- 20
MIN.QCOV.FOR.DOMAIN <- 0
MIN.SCOV.FOR.DOMAIN <- 0.7
MIN.PROB.FOR.DOMAIN <- 70
MAX.EVAL.FOR.DOMAIN <- 1e-2

# parameters ECF-domain hits
MIN.ECF.HIT.LENGTH <- 20
MIN.PROB.FOR.ECF <- 70
MAX.EVAL.FOR.ECF <- 1e-1
MIN.OVERLAP.LENGTH <- 20
MIN.PROB.FOR.ECF.QCOV.SCOV <- 95
MIN.OVERLAP.THRESHOLD <- 0.8

# parameters for comparing ECFs to ECODs:
# 1. domain equivalence
MIN.PROB.FOR.DOMAIN.EQUIVALENCE <- 95
MIN.QCOV.FOR.DOMAIN.EQUIVALENCE <- 0.9
MIN.SCOV.FOR.DOMAIN.EQUIVALENCE <- 0.9
MAX.EVAL.FOR.DOMAIN.EQUIVALENCE <- 10^(-3)
# 2. subdomain
MIN.PROB.FOR.SUBDOMAIN <- 95
MIN.QCOV.FOR.SUBDOMAIN <- 0.7
MAX.EVAL.FOR.SUBDOMAIN <- 10^(-3)
MAX.SCOV.FOR.SUBDOMAIN <- 0.6

# parameters for sECFs
ECF.LENGTH.DIFF.THR <- 5
ECF.MAX.COV.THR <- 0.9
ECF.MIN.COV.THR <- 0.1

