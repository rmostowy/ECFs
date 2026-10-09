### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Set up paths and directories
Sys.setenv(PATH = paste("/opt/bio/hhsuite/bin", "/opt/bio/hhsuite/scripts", Sys.getenv("PATH"), sep = ":"),
           HHLIB = "/opt/bio/hhsuite")

db.a3m <- file.path(database.dir, "a3m")
make.dir(db.a3m)

#####################################
# 1. Load and edit individual a3m files
#####################################
ecf.a3m.list <- list.files(ecfs.a3m.dir, full.names = T, recursive = F)
# ecf.full.table.list <- lapply(1:length(ecf.a3m.list), function(ecf.index){
for(ecf.index in 1:length(ecf.a3m.list)){
  this.ecf.file.path <- ecf.a3m.list[ecf.index]
  this.ecf.file.name <- basename(this.ecf.file.path)
  this.ecf.name = sub("_final\\.a3m$", "", this.ecf.file.name)
  this.ecf.a3m <- readBStringSet(this.ecf.file.path)
  names(this.ecf.a3m)[1] <- this.ecf.name
  this.ecf.a3m.outfile <- sprintf("%s.a3m", this.ecf.name)
  writeXStringSet(this.ecf.a3m, filepath = file.path(db.a3m, this.ecf.a3m.outfile), width = max(width(this.ecf.a3m)))  
}

#####################################
# 2. Create a HHsuite database
#####################################

db.name     <- "ecf_db"
db.hhsuite  <- file.path(database.dir, "hhsuitedb")   # final location (Dropbox)
make.dir(db.hhsuite)

# local scratch space for the build (outside Dropbox)
tmp.dir     <- path.expand(file.path("~/Data/tmp-ecfs", "hhsuite-build"))
tmp.a3m.dir <- file.path(tmp.dir, "a3m")
tmp.db.dir  <- file.path(tmp.dir, "db")
unlink(tmp.dir, recursive = TRUE)
dir.create(tmp.a3m.dir, recursive = TRUE)
dir.create(tmp.db.dir,  recursive = TRUE)

hhlib     <- Sys.getenv("HHLIB")
db.prefix <- file.path(tmp.db.dir, db.name)   # all build files use this prefix

# helpers
run.cmd <- function(cmd){
  message("> ", cmd)
  if(system(cmd) != 0) stop("Command failed: ", cmd, call. = FALSE)
}
db.file <- function(set, ext) shQuote(sprintf("%s_%s.%s", db.prefix, set, ext))

# sanity checks
stopifnot(all(nzchar(Sys.which(c("ffindex_build", "ffindex_apply", "ffindex_order", "hhmake", "cstranslate")))))
stopifnot(file.exists(file.path(hhlib, "data", c("cs219.lib", "context_data.lib"))))

# 2.1 copy all A3M files to the temp directory
a3m.files <- list.files(db.a3m, pattern = "\\.a3m$", full.names = TRUE)
n.a3m <- length(a3m.files)
stopifnot(n.a3m > 0)
stopifnot(all(file.copy(a3m.files, tmp.a3m.dir)))
stopifnot(length(list.files(tmp.a3m.dir)) == n.a3m)
message(n.a3m, " A3M files copied to ", tmp.a3m.dir)

# 2.2 pack the A3M files into an ffindex database
run.cmd(sprintf("ffindex_build -s %s %s %s/",
                db.file("a3m", "ffdata"), db.file("a3m", "ffindex"), shQuote(tmp.a3m.dir)))

# 2.3 build HMM profiles
run.cmd(sprintf("ffindex_apply %s %s -i %s -d %s -- hhmake -i stdin -o stdout -v 0 > /dev/null",
                db.file("a3m", "ffdata"), db.file("a3m", "ffindex"),
                db.file("hhm", "ffindex"), db.file("hhm", "ffdata")))

# 2.4 build context-state sequences for the HHblits prefilter
run.cmd(sprintf("cstranslate -f -x 0.3 -c 4 -I a3m -b -i %s -o %s -A %s -D %s > /dev/null",
                shQuote(paste0(db.prefix, "_a3m")), shQuote(paste0(db.prefix, "_cs219")),
                shQuote(file.path(hhlib, "data", "cs219.lib")),
                shQuote(file.path(hhlib, "data", "context_data.lib"))))

# 2.5 sort all three sets by length (as HHblits expects)
sorting.file <- file.path(tmp.db.dir, "sorting.dat")
run.cmd(sprintf("sort -k3 -n -r %s | cut -f1 > %s", db.file("cs219", "ffindex"), shQuote(sorting.file)))
for(set in c("a3m", "hhm", "cs219")){
  ord.data  <- sprintf("%s_%s_ord.ffdata",  db.prefix, set)
  ord.index <- sprintf("%s_%s_ord.ffindex", db.prefix, set)
  run.cmd(sprintf("ffindex_order %s %s %s %s %s",
                  shQuote(sorting.file), db.file(set, "ffdata"), db.file(set, "ffindex"),
                  shQuote(ord.data), shQuote(ord.index)))
  file.rename(ord.data,  sprintf("%s_%s.ffdata",  db.prefix, set))
  file.rename(ord.index, sprintf("%s_%s.ffindex", db.prefix, set))
}
file.remove(sorting.file)

# 2.6 check: all three sets must have as many entries as input files
n.entries <- sapply(c("a3m", "hhm", "cs219"), function(set)
  length(readLines(sprintf("%s_%s.ffindex", db.prefix, set))))
print(n.entries)
stopifnot(all(n.entries == n.a3m))

# 2.7 copy the finished database to its final location
db.files <- list.files(tmp.db.dir, pattern = paste0("^", db.name, "_.*\\.ff(data|index)$"),
                       full.names = TRUE)
stopifnot(length(db.files) == 6)
old <- list.files(db.hhsuite, pattern = paste0("^", db.name, "_"), full.names = TRUE)
if(length(old) > 0) file.remove(old)
stopifnot(all(file.copy(db.files, db.hhsuite)))

# verify the copy by comparing checksums
md5.tmp   <- tools::md5sum(db.files)
md5.final <- tools::md5sum(file.path(db.hhsuite, basename(db.files)))
stopifnot(all(unname(md5.tmp) == unname(md5.final)))
message("Database saved to ", db.hhsuite)

# 2.8 test search against the FINAL database
query <- a3m.files[1]
hhr   <- file.path(tmp.dir, "test.hhr")
run.cmd(sprintf("hhblits -i %s -d %s -o %s -n 1 > /dev/null 2>&1",
                shQuote(query), shQuote(file.path(db.hhsuite, db.name)), shQuote(hhr)))
cat(head(readLines(hhr), 15), sep = "\n")

# 2.9 clean up the scratch space (only reached if everything above succeeded)
unlink(tmp.dir, recursive = TRUE)
