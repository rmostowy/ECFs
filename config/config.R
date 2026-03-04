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

#####################
## Update PATH
#####################
bio_bin_dirs <- list.dirs("/opt/bio", recursive = FALSE, full.names = TRUE)
bio_bin_dirs <- file.path(bio_bin_dirs, "bin")
bio_bin_dirs <- bio_bin_dirs[dir.exists(bio_bin_dirs)]
Sys.setenv(PATH = paste(Sys.getenv("PATH"), paste(bio_bin_dirs, collapse = ":"), sep = ":"))

##############################
### Load functions
##############################
functions.repo.path <- file.path(main.path, "Github", "proj-rmostowy", "functions")
scripts.functions.path <- file.path(root, "scripts", "functions")

functions.basic.filename <- file.path(functions.repo.path, "basic.R")
safe_source(functions.basic.filename)

functions.custom.cov.filename <- file.path(functions.repo.path, "custom", "hhr-cov.R")
safe_source(functions.custom.cov.filename)

helper.functions.filename <- file.path(scripts.functions.path, "helper.R")
safe_source(helper.functions.filename)

gpt.helper.functions.filename <- file.path(scripts.functions.path, "gpt-helper.R")
safe_source(gpt.helper.functions.filename)

gggenomes.functions.filename <- file.path(scripts.functions.path, "gggenomes.R")
safe_source(gggenomes.functions.filename)

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
ecod.db.path <- file.path(mgg.db, "hhpred", "ECOD", "ECOD_F70_20230309")
ecod.domains.filename <- file.path(ecod.db.path, "ecod.develop288.F70.domains.txt")
ecod.domains.fasta.filename <- file.path(ecod.db.path, "ecod.develop288.F70.fasta.txt")

ecod.headers <- c("uid", "ecod_domain_id","manual_rep","f_id","pdb","chain","pdb_range","seqid_range","unp_acc",
  "arch_name","x_name","h_name","t_name","f_name", "asm_status", "ligand")

#############################
### Define PFAM data
#############################
pfam.db.path <- file.path(mgg.db, "hhpred", "pfamA", "Pfam37")
pfam.domains.filename <- file.path(pfam.db.path, "pfamA.txt")
pfam.sql.filename <- file.path(pfam.db.path, "pfamA.sql")

#############################
### Define main data paths
#############################

ecf.projects.datasets <- file.path(mgg.projects, "ECFs", "datasets")
antidefence.data.path <- file.path(ecf.projects.datasets, "antidefence")
refseq.input.data.path <- file.path(ecf.projects.datasets, "refseq-input", "ecf-input-dataset")

ecf.explorer.path <- file.path(mgg.projects, "ECFs/ecf-explorer/ecf-explorer_2025-10-24")
ecf.raw.input.path <- file.path(ecf.explorer.path, "input")
ecf.raw.data.path <- file.path(ecf.explorer.path, "output")


#############################
### Define reprseq/rHMM data
#############################

all.cds.filename <- file.path(refseq.input.data.path, "refseq-aa.fa")
metadata.table.filename <- file.path(refseq.input.data.path, "refseq-metadata.tsv")
refseq.protein.stats.table.filename <- file.path(refseq.input.data.path, "refseq-protein-stats.txt")
refseq.genomes.filename <- file.path(refseq.input.data.path, "refseq-genomes.fa")
name.table.origin.filename <- file.path(ecf.raw.data.path, "prot-families", "representative", "name-table.txt")
refseq.protein.clustering.filename <- file.path(ecf.raw.data.path, "prot-families", "representative", "clustering.tsv")

reprseq.lengths.filename <- file.path(ecf.raw.data.path, "prot-families", "representative", "repr-seqs-lengths.txt")
phrog.hhpred.hits.filename <- file.path(ecf.raw.data.path, 'table-scan-rseq_phrogs.txt')
antidefence.descr.filename <- file.path(antidefence.data.path, "8_Bruria2024_seq_v2-RM.tsv")
antidefence.seqs.filename <- file.path(antidefence.data.path, "8_Bruria2024_seq.fasta")
antidefence.hhpred.hits.filename <- file.path(ecf.raw.data.path, 'table-scan-ad_rseq.txt')

# reprseq.families.raw.filename <- file.path(minimal.input.data.path, "python_pipeline", "output", "prot-families", "families", "dataset-full", "mcl", "repr-hits-pairwise-prob95-cov80-mcl-I20.out") # this needs to be recreated

####################
### Define ECF data
####################

ecf.hits.filename <- file.path(ecf.raw.data.path, "table-scan-rseq_ecf.txt")
ecf.hits.vs.ecod.filename <- file.path(ecf.raw.data.path, "table-scan-ecfs_ecod.txt")
ecf.hits.vs.pfam.filename <- file.path(ecf.raw.data.path, "table-scan-ecfs_pfam.txt")
ecf.hits.vs.ecf.filename <- file.path(ecf.raw.data.path, "table-scan-ecfs_ecf.txt")
reprseq.vs.reprseq.filename <- file.path(ecf.raw.input.path, "table-rseq-hhr.txt")
reprseq.vs.reprseq.clust.filename <- file.path(ecf.raw.data.path, "reprseq-all-clust.txt")
ecod.domain.hits.filename <- file.path(ecf.raw.data.path, "table-scan-rseq_ecod.txt")
pfam.domain.hits.filename <- file.path(ecf.raw.data.path, "table-scan-rseq_pfam.txt")
qecf.stats.filename <- file.path(ecf.raw.data.path, "reprseq-all-stat-coverage.txt")
ecf.vs.qecf.filename <- file.path(ecf.raw.data.path, "ecf2qecfs.csv")
ecfs.a3m.dir <- file.path(ecf.raw.data.path, "ecfs")

ps.vs.ecf.filename <- file.path(ecf.raw.data.path, "table-scan-ps_ecf.txt")
ps.vs.ecod.filename <- file.path(ecf.raw.data.path, "table-scan-ps_ecod.txt")
ps.vs.pfam.filename <- file.path(ecf.raw.data.path, "table-scan-ps_pfam.txt")

####################
### Define other data
####################

ictv.release.filename <- file.path(mgg.db, "ICTV", "VMR", "VMR_MSL39_v2.xlsx")
phagescope.protein.metadata.filename <- file.path(ecf.projects.datasets, "phagescope-random", "used-in-paper", "PhageScope_protein_metadata_2000.tsv")
phagescope.protein.seqs.filename <- file.path(ecf.projects.datasets, "phagescope-random", "used-in-paper", "PhageScope_protein_2000.fasta")

#######################
### Data directories
#######################
data.tmp.dir <- file.path(main.path, "Data", "tmp-ecfs")

main.data.dir <- file.path(main.path, "Dropbox", "Projects", "Data", "ECFs", "v4")
data.input.dir <- file.path(main.data.dir, "input")
data.output.dir <- file.path(main.data.dir, "output")
visual.dir <- file.path(main.data.dir, "visualisations")

data.phrog.dir <- file.path(data.output.dir, "annotation", "phrogs")
data.domain.dir <- file.path(data.output.dir, "annotation", "domains")
data.ps.dir <- file.path(data.output.dir, "annotation", "phagescope")
data.ecfs.dir <- file.path(data.output.dir, "ecfs")

data.comparisons.dir <- file.path(data.output.dir, "comparisons")
ecf.network.dir <- file.path(data.comparisons.dir, "networks")

# ECF project
ecf.proj.dir <- file.path(main.data.dir, "paper-ecf")
ecf.figures.dir <- file.path(ecf.proj.dir, "figures")
ecf.suppl.figures.dir <- file.path(ecf.figures.dir, "suppl-figs")

dirs2create <- c(data.tmp.dir,
                 main.data.dir, data.input.dir, data.output.dir, visual.dir,
                 data.phrog.dir, data.domain.dir, data.ps.dir, data.ecfs.dir,
                 data.comparisons.dir, ecf.network.dir,
                 ecf.proj.dir, ecf.figures.dir, ecf.suppl.figures.dir)

make.dirs(dirs2create)

#######################
### Output files
#######################

# input data
reprseq.table.output <- file.path(data.input.dir, "protein-name-table.txt")
reprseq.vs.reprseq.cov.output <- file.path(data.input.dir, "reprseq-vs-reprseq-cov.txt")
reprseq.seqs.filename <- file.path(data.input.dir, "reprseq-seq.fa")
# reprseq.removed.filename <- file.path(data.input.dir, "reprseq-removed.txt")

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
ecf.table.output <- file.path(data.ecfs.dir,  "ecf-table.txt")
ecf.vs.qecf.output <- file.path(data.ecfs.dir, "ecfs-vs-qecfs.txt")

reprseq.ecfs.output <- file.path(data.ecfs.dir, "reprseq-ecfs.txt")
reprseq.ecfs.all.output <- file.path(data.ecfs.dir, "reprseq-ecfs-all.txt") 

ecfs.vs.ecfs.output <- file.path(data.ecfs.dir, "ecfs-vs-ecfs.txt")
ecfs.vs.ecod.output <- file.path(data.ecfs.dir, "ecfs-vs-ecod.txt")
ecfs.vs.pfam.output <- file.path(data.ecfs.dir, "ecfs-vs-pfam.txt")
ecfs.vs.ecod.hhsenser.output <- file.path(data.ecfs.dir, "ecfs-vs-ecod-hhsenser.txt")
ecfs.vs.pfam.hhsenser.output <- file.path(data.ecfs.dir, "ecfs-vs-pfam-hhsenser.txt")

ecfs.vs.ecfs.cov.output <- file.path(data.ecfs.dir, "ecfs-vs-ecfs-cov.txt")
# ecf.lengths.output <- file.path(data.ecfs.dir, "ecf-lengths.txt")

ecf.family.table.output <- file.path(data.ecfs.dir, "ecfs-family-table.txt")
# ecf.family.network.output <- file.path(data.ecfs.dir, "ecfs-network.txt")

# metadata
metadata.table.output <- file.path(data.input.dir, "metadata.tsv")
accessions.cds.output <- file.path(data.input.dir, "ncbi-accessions-cds.tsv")

# phagescope
phagescope.sequences.output <- file.path(data.ps.dir, "phagescope-seqs.fasta")
phagescope.metadata.output <- file.path(data.ps.dir, "phagescope-metadata.txt")
phagescope.vs.ecod <- file.path(data.ps.dir, "phagescope-vs-ecod.txt")
phagescope.vs.pfam <- file.path(data.ps.dir, "phagescope-vs-pfam.txt")
phagescope.vs.ecf <- file.path(data.ps.dir, "phagescope-vs-ecf.txt")

# comparisons
reprseq.vs.reprseq.blast.output <- file.path(data.comparisons.dir, "reprseq-vs-reprseq-blast.txt")
# reprseq.vs.reprseq.blast.filtered.output <- file.path(data.comparisons.dir, "reprseq-vs-reprseq-blast-filtered.txt")
reprseq.vs.reprseq.blast.filtered.unordered.output <- file.path(data.comparisons.dir, "reprseq-vs-reprseq-blast-filt-unord.txt")
# reprseq.ecf.vs.ecod.output <- file.path(data.comparisons.dir, "reprseq_ecf-vs-ecod.txt")
# reprseq.ecf.vs.pfam.output <- file.path(data.comparisons.dir, "reprseq_ecf-vs-pfam.txt")
# reprseq.ecf.vs.ecod.and.pfam.output <- file.path(data.comparisons.dir, "reprseq_ecf-vs-ecod-and-pfam.txt")

ecf.classification.table.output <- file.path(data.comparisons.dir, "ecf-classifications.txt")

# ecf-sharing pairs
ecf.sharing.pairs.output <- file.path(data.comparisons.dir, "ecf-sharing-pairs.txt")

# REM table
rem.table.output <- file.path(data.comparisons.dir, "rem-table.txt")

# Genome all-by-all megablast table
genomes.aba.table.output <- file.path(data.comparisons.dir, "genomes-aba.txt")

#################
## ECF project dirs
#################




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
# ACCEPTED.NUM.UNKNOWN.AA <- 10
MIN.PROB.FOR.REPRSEQ.FAM <- 95
REPRSEQ.FAMILY.QCOV.THR <- 0.5
REPRSEQ.FAMILY.SCOV.THR <- 0.5
REPRSEQ.FAMILY.INFLATION.FACTOR <- 1.3

MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION <- 30
MIN.NUM.PHROG.SEQ.PER.ANNOT <- 500
MIN.NUM.PROT.PER.ANNOT <- 100 # Bogna had MIN.NUM.PROT
# MAIN.COVS.FOR.ANNOTATION = c(0.1, 0.3, 0.8)
DEFAULT.MINIMUM.COV.FOR.ANNOTATION <- 60
DEFAULT.MINIMUM.PROB.FOR.ANNOTATION <- 95
DEFAULT.MAXIMUM.EVAL.FOR.ANNOTATION <- 1e05
UNSPECIFIC.ANNOTATIONS = c("tail", "structural protein", "uknown")

# parameters related to ECF filtering
MIN.PERCENTAGE.OF.INFORMATIVE.BASES <- 0.5
MIN.HIT.LENGTH.FOR.ECF.FILTERING <- 20
MIN.SCOV.FOR.ECF.FILTERING <- 0.7
MIN.PROB.FOR.ECF.FILTERING <- 95
MAX.EVAL.FOR.ECF.FILTERING <- 1e-3
MIN.NO.REPRSEQ.PER.ECF <- 1

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
MAX.EVAL.FOR.ECF.HHSENSER <- 10
MIN.OVERLAP.LENGTH <- 20
MIN.PROB.FOR.ECF.QCOV.SCOV <- 95
MIN.OVERLAP.THRESHOLD <- 0.8

# parameters for comparing ECFs to ECOD/PFAM:
# 1. domain equivalence
MIN.PROB.FOR.DOMAIN.EQUIVALENCE <- 95
MIN.QCOV.FOR.DOMAIN.EQUIVALENCE <- 0.8
MIN.SCOV.FOR.DOMAIN.EQUIVALENCE <- 0.8
MAX.EVAL.FOR.DOMAIN.EQUIVALENCE <- 1e-3
# 2. subdomain
MIN.PROB.FOR.SUBDOMAIN <- 95
MIN.QCOV.FOR.SUBDOMAIN <- 0.7
MAX.EVAL.FOR.SUBDOMAIN <- 10^(-3)
MAX.SCOV.FOR.SUBDOMAIN <- 0.6

# parameters for sECFs
ECF.FAMILY.QCOV.THR <- 0.5
ECF.FAMILY.SCOV.THR <- 0.5
ECF.FAMILY.INFLATION.FACTOR <- 1.3

# recent recombination parameters
RECENT.REC.MIN.PIDENT <- 50       # check if still used
RECENT.REC.MAX.QCOV <- 0.7        # check if still used
RECENT.REC.MAX.SCOV <- 0.7        # check if still used
RECENT.REC.MIN.Q.HIT.LENGTH <- 20 # check if still used
RECENT.REC.MIN.S.HIT.LENGTH <- 20 # check if still used
# 
ECF.OVERLAP.SCORE.THR <- 0.9
ECF.MIN.OVERLAP.SCORE.THR <- 0.9
MIN.NO.PROTEINS.WITH.HIGH.OVERLAPPING.ECF <- 3
RE.ECF.NO.PAIRS.THR <- 20

# other parameters
GENOME.MEGABLAST.EVALUE=1e-5
DOMAIN.CONTAINMENT.COV.THR <- 0.9
CONTAINMENT.PROB.THR <- 99
RE.MAX.DISTANCE=0.3
RE.MIN.CLUSTER.SIZE=5
RE.MIN.OVERLAP=0.8
NO.GENES.TO.EXTEND <- 4L
PAD.PARM=10L
# MIN.CDS.OVERLAP.WITH.NCBI <- 0.5




