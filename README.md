# ECFs: analysis code for Evolutionarily Conserved Fragments in bacteriophage proteins

This repository contains the R scripts behind the manuscript on **Evolutionarily Conserved Fragments (ECFs)** in bacteriophage proteins. The scripts take HMM-HMM comparison outputs produced by colleagues (the `ecf-explorer` run), prepare and annotate the data, build the ECF / ECF-family resource, generate the publication figures, and build the output database (an HH-suite database of ECF profiles plus an ECF metadata table).

The repository holds **code only**. All input data, intermediate tables, figures and databases live outside the repository (see [What is not in the repository](#what-is-not-in-the-repository)).

This README describes how the repository is organised. It is a description, not a tutorial. It covers the current code only.

## Citation
Bogna J. Smug\*, Krzysztof Szczepaniak\*, Jade Leconte, Stanislaw Dunin-Horkawicz & Rafal J. Mostowy. 
*Evolutionarily Conserved Fragments reveal modular architecture and expand 
functional coverage in bacteriophage proteins.* Manuscript in preparation.
(\*equal contribution)

Please also cite the archived release you used. Each release of this repository is archived on Zenodo, which assigns a DOI to that version (shown on the release page and in the Zenodo record); the "Cite this repository" button on GitHub (from `CITATION.cff`) gives the citation. See [Releases, versions and DOIs](#releases-versions-and-dois).

## Releases, versions and DOIs

Code and database are released together, in a **single version series** (`vMAJOR.MINOR.PATCH`, e.g. `v1.0.2`):

- A release is a git tag of this repository. Its automatic source archive is the **code**; from `v1.0.2` the same release also carries the **database files** as downloadable assets (`ecf-hhsuitedb.tar.gz`, `ecf-a3m.tar.gz`, `ecf-metadata.csv`, `SHA256SUMS.txt`). `v1.0.1` contained the code only.
- **DOIs come from Zenodo.** Publishing a release makes Zenodo archive the source code of the tag and issue a DOI for that version, plus one "concept" DOI that covers all versions. GitHub itself assigns no DOI, and files attached to a release are not part of the Zenodo archive.
- The database is additionally deposited at FigShare (<https://doi.org/10.6084/m9.figshare.34887162>).
- The text of each release is in `release-notes/<version>.md`. The release files are built with `scripts/4_prep-db/3_make-release.sh` and uploaded as a draft with `scripts/4_prep-db/4_push-release.sh`; the draft is then published by hand, which triggers the Zenodo archiving.
- Author and title metadata for the archive are in `CITATION.cff`.

## Terminology used in the scripts

| Term | Meaning |
|---|---|
| **reprseq** (`reprseq000123`) | Representative sequence of a cluster of RefSeq phage proteins. Its HMM is referred to as an **rHMM** in the figure scripts. Reprseqs are also grouped into protein families (`fam…`) by MCL. |
| **ECF** (`ecf_00004`) | Evolutionarily conserved fragment. Each ECF is defined by an A3M alignment (`<ecf>_final_filter.a3m`) in the ecf-explorer output; the first sequence is the ECF representative (`er.name`, `reprseqXXXXXX_suffix`). |
| **EF** (`EF0007`) | ECF family: a cluster of ECFs obtained with MCL from ECF-vs-ECF coverage. |
| **qcov / scov** | Fraction of the query / subject covered by a hit. |
| **ECOD / Pfam** | Reference domain databases (ECOD F70 release 20230309, Pfam 37). ECOD hierarchy levels are X, H, T, F. |
| **PHROGs** | Phage protein families with functional annotation; used, with an antidefence set, to annotate reprseqs by function. |
| **domain-dark** | A reprseq with ECF hits but no Pfam or ECOD hit (Fig. 5). |
| **REM table** | Table of reprseq pairs sharing partial similar segments, described in the code as putative recent recombinations (see [REM analysis](#rem-analysis-not-in-the-manuscript)). |

## Directory map

```
.
├── rafals-project.Rproj          RStudio project file; marks the project root
├── config/
│   └── config.R                  paths, output file names, parameters, labels, function loading
└── scripts/
    ├── 1_data-preparation/
    │   ├── 1_get-main/           8 scripts: raw tables -> annotated ECF / EF resource
    │   └── 2_compare/            5 scripts: classification, MMseqs2 all-vs-all, REM table,
    │                             genome megablast, protein architecture metrics
    ├── 2_analysis/
    │   ├── 1_diversity-modularity/   exploratory plots (precursors of Fig. 1-2)
    │   ├── 2_functional-cov/         exploratory plots (precursors of Fig. 3)
    │   ├── 3_novel-ecfs/             ECF network edge list, ECF-vs-domain exploration
    │   ├── 4_ecf2ecod/               ECF-vs-ECOD comparison categories
    │   └── 5_ecf-vs-rem/analysis_01-23/   REM analysis (not in the manuscript)
    ├── 3_ecf-paper/
    │   ├── make-plots-template.R     skeleton for a new figure script
    │   ├── part1_basic-prop/         Fig1.R
    │   ├── part2_ecfs-vs-domains/    Fig2.R
    │   ├── part3_funct-cov/          Fig3.R
    │   ├── part4_arch-regimes/       Fig4.R
    │   ├── part5_funct-rescue/       Fig5.R
    │   └── part6_rbp-struct/         make1-shortlist.R, make2-plot.R, make3-pymol.sh, make4-run.sh
    ├── 4_prep-db/                HH-suite database and ECF metadata table (R); release packaging and upload (bash)
    └── functions/                helper.R, gpt-helper.R, gggenomes.R
```

Every R script is standalone. It loads its libraries, locates the project root with `rprojroot::find_rstudio_root_file()`, sources `config/config.R`, reads upstream tables from disk, and writes its own outputs. The R scripts do not call each other (apart from local helper files in the REM analysis); the dependencies between them are through files on disk (see [Run order](#run-order)).

## Configuration: `config/config.R`

`config.R` is sourced by every script. It:

- defines all input locations and output file names as R variables (`*.filename` for inputs, `*.output` for tables written by the scripts), so scripts refer to files by variable rather than by path;
- defines the data directory tree and creates the directories when the file is sourced;
- appends the `bin` directory of every tool under `/opt/bio/` to `PATH`;
- sources shared functions: two files from the separate functions repository (`basic.R`, `custom/hhr-cov.R`) and `helper.R`, `gpt-helper.R`, `gggenomes.R` from `scripts/functions/`;
- holds all analysis parameters in one place (HMM hit filtering thresholds, PHROG annotation thresholds, MCL inflation factors and coverage thresholds for reprseq families and EFs, ECF-vs-domain relationship thresholds, and parameters for the REM analysis);
- holds the PHROG functional classes with their colours and display labels, and the ECOD architecture classes;
- ends by checking that every `*.filename` input exists.

### Data layout written by the scripts

All generated data go under `main.data.dir` (`<data>` below):

```
<data>/
├── input/                          protein-name-table.txt, reprseq-seq.fa, reprseq-vs-reprseq-cov.txt, metadata.tsv
├── output/
│   ├── annotation/
│   │   ├── phrogs/                 PHROG / antidefence tables, reprseq functional annotations
│   │   ├── domains/                ECOD and Pfam tables, reprseq-vs-domain hits
│   │   └── phagescope/             PhageScope sequences, metadata and hits
│   ├── ecfs/                       ECF table, hit tables, ECF family table
│   ├── comparisons/                ECF classification, all-vs-all tables, REM table, architecture metrics, networks/
│   └── database/                   a3m/, hhsuitedb/, ecf-metadata.csv
├── visualisations/                 plots from scripts/2_analysis
└── paper-ecf/
    ├── figures/                    figure1.png ... figure5.png
    │   └── suppl-figs/             suppl-<N>p<k>_*.png
    └── part*/                      intermediate plots and tables of each figure script
```

Scratch files (MCL, MMseqs2, BLAST, HH-suite builds) go to `~/Data/tmp-ecfs`.

## Data preparation (`scripts/1_data-preparation`)

In the lists below, **Reads** names inputs from outside the repository as `[raw]` and tables written by other scripts by their file name.

### `1_get-main`

**`1_get-metadata.R`** reads the RefSeq metadata table and the ICTV VMR release (`VMR_MSL39_v2.xlsx`) and writes `input/metadata.tsv`: genome-level metadata (length, GC, CDS count, host, NCBI family and genus) joined with ICTV taxonomy.

**`2_process-reprseqs.R`** reads `[raw]` the ecf-explorer name table, protein clustering table, reprseq-vs-reprseq HHR table, the RefSeq protein sequences and Prodigal statistics. It builds the reference table of all proteins and their reprseqs, extracts the reprseq sequences, computes pairwise reprseq qcov / scov, and clusters reprseqs into protein families with MCL. It writes `input/protein-name-table.txt`, `input/reprseq-seq.fa` and `input/reprseq-vs-reprseq-cov.txt`. This is the root of most downstream scripts.

**`3_process-phrogs.R`** reads `[raw]` the PHROG annotation tables, the manual PHROG mapping workbook (`v4_phrogs-table-2024.xlsx`), and the antidefence description table and FASTA. It merges the PHROG tables, removes singleton and unknown-function PHROGs, remaps functions and classes using the manual workbook, and counts sequences and PHROGs per annotation. It writes `annotation/phrogs/phrogs-table.txt`, `phrogs-count.txt` and `antidefence-table.txt`.

**`4_annotate.R`** reads `protein-name-table.txt`, the three tables from `3_process-phrogs.R`, and `[raw]` the reprseq-vs-PHROG and antidefence-vs-reprseq HH-suite hit tables. It filters hits (alignment length, probability, qcov / scov), attaches functions and classes, flags frequently used ("top") annotations, and selects the best annotation per reprseq. It writes `annotation/phrogs/reprseq-annot-all.txt` and `reprseq-annot-unique.txt`.

**`5_get-domains.R`** reads `protein-name-table.txt` and `[raw]` the ECOD F70 domain table, the Pfam 37 `pfamA.txt` / `pfamA.sql`, and the reprseq-vs-ECOD and reprseq-vs-Pfam hit tables. It adds hierarchy ids, architecture class and length to the ECOD table, simplifies the Pfam table, computes qcov / scov of the hits and filters them. It writes `annotation/domains/ecod-domains.txt`, `pfam-domains.txt`, `reprseq-domains-ecod.txt` and `reprseq-domains-pfam.txt`.

**`6_get-ecfs.R`** reads `protein-name-table.txt` and `[raw]` the ECF alignments (`ecfs/*_final_filter.a3m`) and the reprseq-vs-ECF, ECF-vs-ECF, ECF-vs-ECOD and ECF-vs-Pfam hit tables. It builds the ECF table (length, number of sequences, X-residue counts), filters the hit tables, keeps less strictly filtered (HHsenser-style) ECF-vs-ECOD / Pfam versions, and computes pairwise ECF-ECF qcov / scov. It writes to `output/ecfs/`: `ecf-table.txt`, `reprseq-ecfs.txt`, `ecfs-vs-ecfs.txt`, `ecfs-vs-ecfs-cov.txt`, `ecfs-vs-ecod.txt`, `ecfs-vs-pfam.txt`, `ecfs-vs-ecod-hhsenser.txt` and `ecfs-vs-pfam-hhsenser.txt`.

**`7_get-EFs.R`** reads `ecf-table.txt` and `ecfs-vs-ecfs-cov.txt`, clusters ECFs into families with MCL (edges weighted by qcov × scov × probability), and writes `output/ecfs/ecfs-family-table.txt` (ECF to `EF####`).

**`8_process-ps.R`** reads `[raw]` the PhageScope protein metadata and FASTA and the PhageScope-vs-ECF, -ECOD and -Pfam hit tables, plus `ecf-table.txt`. It filters the hits and writes `annotation/phagescope/phagescope-seqs.fasta`, `phagescope-metadata.txt`, `phagescope-vs-ecf.txt`, `phagescope-vs-ecod.txt` and `phagescope-vs-pfam.txt`.

### `2_compare`

**`1_classify-ecfs.R`** reads `protein-name-table.txt`, `metadata.tsv`, `ecod-domains.txt`, `reprseq-annot-unique.txt`, `ecf-table.txt`, `reprseq-ecfs.txt` and `ecfs-vs-ecod.txt`. For each ECF it assigns the PHROG function and class from the reprseqs carrying it, classifies its relationship to ECOD (`new`, `similar`, `equivalent`, `subdomain`), records the best ECOD architecture class and T-level, and marks the ECF as having function only (F), domain only (D), both (F+D) or neither. It writes `output/comparisons/ecf-classifications.txt`.

**`2_blast-reprseqs.R`** reads `reprseq-seq.fa` and `protein-name-table.txt`. It runs MMseqs2 all-vs-all on the reprseqs (scratch space under `~/Data/tmp-ecfs/recent-rec`) and recomputes query and subject coverage as the union over all alignments of a pair. It writes `comparisons/reprseq-vs-reprseq-blast.txt`.

**`3_get-rem-table.R`** reads `reprseq-vs-reprseq-blast.txt` and `protein-name-table.txt`. It removes pairs with overlapping or more than two alignments, reduces each pair to one direction, and selects pairs with a partial local similarity using the `RECENT.REC.*` parameters. It writes `comparisons/reprseq-vs-reprseq-blast-filt-unord.txt` and `comparisons/rem-table.txt`.

**`4_megablast-genomes.R`** reads `[raw]` the RefSeq genomes FASTA, builds a nucleotide BLAST database in a scratch directory, runs all-vs-all megablast, removes the scratch directory, and writes `comparisons/genomes-aba.txt`.

**`5_calc-arch-metrics.R`** reads `ecf-table.txt`, `ecfs-family-table.txt`, `reprseq-ecfs.txt` and `[raw]` the reprseq lengths (it also loads the reprseq annotation tables). Per reprseq it computes the number of ECF hits, covered fraction, number of EFs, mean effective number of EFs per position, two segregation indices (range overlap and barycentre distance) and the number of dominant-family segments. It writes `comparisons/reprseq-arch-metrics.txt`.

### Run order

The script numbering does not fully reflect dependencies, so the order is given here from which outputs each script reads. Scripts in the same step do not depend on each other; each step needs only outputs of earlier steps.

| Step | Script(s) | Needs (script outputs) |
|---|---|---|
| 1 | `1_get-main/1_get-metadata.R`, `1_get-main/2_process-reprseqs.R`, `1_get-main/3_process-phrogs.R`, `2_compare/4_megablast-genomes.R` | raw inputs only |
| 2 | `1_get-main/4_annotate.R` | `2_process-reprseqs`, `3_process-phrogs` |
| 2 | `1_get-main/5_get-domains.R`, `1_get-main/6_get-ecfs.R`, `2_compare/2_blast-reprseqs.R` | `2_process-reprseqs` |
| 3 | `1_get-main/7_get-EFs.R`, `1_get-main/8_process-ps.R` | `6_get-ecfs` |
| 3 | `2_compare/3_get-rem-table.R` | `2_blast-reprseqs` |
| 3 | `2_compare/1_classify-ecfs.R` | `1_get-metadata`, `2_process-reprseqs`, `4_annotate`, `5_get-domains`, `6_get-ecfs` (`3_process-phrogs` indirectly; not `7` or `8`) |
| 4 | `2_compare/5_calc-arch-metrics.R` | `4_annotate`, `6_get-ecfs`, `7_get-EFs` |

Downstream, the figure scripts and the metadata script use these outputs:

- **Fig. 1** reads outputs of `2_process-reprseqs` (`reprseq-seq.fa`), `5_get-domains`, `6_get-ecfs` and `7_get-EFs`.
- **Fig. 2** reads `4_annotate`, `5_get-domains`, `6_get-ecfs` and `7_get-EFs`.
- **Fig. 3** reads `1_get-metadata`, `2_process-reprseqs`, `4_annotate`, `5_get-domains`, `6_get-ecfs`, `7_get-EFs` and `8_process-ps`.
- **Fig. 4** and **part 6** read `4_annotate`, `6_get-ecfs`, `7_get-EFs` and `5_calc-arch-metrics`; part 6 also reads `2_process-reprseqs`, `5_get-domains` and the ECF-vs-ECOD hits.
- **Fig. 5** reads `2_process-reprseqs`, `4_annotate`, `5_get-domains`, `6_get-ecfs` and `7_get-EFs`.
- **`4_prep-db/2_ecf-metadata.R`** reads `4_annotate`, `5_get-domains`, `6_get-ecfs` and `7_get-EFs`.

The outputs of `2_compare/1_classify-ecfs.R`, `3_get-rem-table.R` and `4_megablast-genomes.R` are not read by any of the figure scripts.

## Exploratory analyses (`scripts/2_analysis`)

Scripts 1-4 are earlier exploratory plot scripts. They write to `<data>/visualisations/part*/`. Their content was later folded into the figure scripts, and they are not needed to produce the figures.

| Script | Content | Output folder |
|---|---|---|
| `1_diversity-modularity/1_plot-distr.R` | reprseqs per ECF, ECFs per EF, ECF / ECOD / Pfam length distributions | `part1_diversity/analysis1_distr` |
| `1_diversity-modularity/2_modularity.R` | high-confidence ECF-vs-ECF qcov / scov and containment transitivity | `part1_diversity/analysis2_modularity` |
| `2_functional-cov/1_plot-refseq.R` | ECF / ECOD / Pfam coverage of RefSeq reprseqs by PHROG function | `part2_funct-cov/analysis1_refseq` |
| `2_functional-cov/2_plot-phagescope.R` | the same comparison for the PhageScope proteins (UpSet and functional coverage) | `part2_funct-cov/analysis2-phagescope` |
| `3_novel-ecfs/1_network.R` | ECF coverage-network edge list (`networks/ecf-vis-network.txt`) | `output/comparisons/networks` |
| `3_novel-ecfs/2_new-ecfs.R` | ECF domain hits by length, diversity, function and taxonomy; writes `ecf-domain-table.txt` | `part3_novel-ecfs` |
| `4_ecf2ecod/1_ecf2ecod.R` | ECF-vs-ECOD comparison categories and ECF-linked ECOD groups | `part4_ecf2ecod` |

### REM analysis (not in the manuscript)

`5_ecf-vs-rem/analysis_01-23/` (scripts `0_` to `8_`, `functions.R`, `plot-ecf-genomes.R` and the `tables/` folder) analyses the REM table from `2_compare/3_get-rem-table.R`: relative positions of the shared segments, overlaps with ECOD / Pfam / ECF hits, ECF-based alignments and trees, and synteny of the corresponding genomic regions (using `genomes-aba.txt`). It was intended for the paper but is not part of the manuscript, and none of Figs. 1-6 depends on it. Its outputs go to `<data>/visualisations/part5_rem/`.

## Figures (`scripts/3_ecf-paper`)

Each figure script reads the tables listed under [Run order](#run-order), saves its intermediate plots to `<data>/paper-ecf/<part>/`, writes the assembled main figure to `paper-ecf/figures/figureN.png` and the supplementary figures to `paper-ecf/figures/suppl-figs/suppl-Np<k>_*.png`. The figure scripts are independent of one another, apart from the part 6 chain below. `make-plots-template.R` is a skeleton with the usual header (libraries, config, output directory, data loading) for a new figure script.

| Script | Content | Main figure | Supplementary |
|---|---|---|---|
| `part1_basic-prop/Fig1.R` | Basic properties: ECF size, reprseqs per ECF, ECFs per EF, length distributions of ECFs vs Pfam vs ECOD, relative-position maps of hits on reprseqs, EFs per reprseq, ECF coverage of domain regions | `figure1.png` (panels tagged C-E) | suppl-1p1 to 1p4 |
| `part2_ecfs-vs-domains/Fig2.R` | ECFs vs domains: nested (containment) structure of the ECF coverage network, asymmetry and a Poisson GLM; shorter "fragment" ECFs linking unrelated ECFs; ECF-to-ECOD / Pfam relationship classes by length and ECOD-Pfam agreement; ECFs bridging ECOD X-groups; worked examples | `figure2.png` (panels A-E) | suppl-2p1 to 2p7 |
| `part3_funct-cov/Fig3.R` | Functional coverage: UpSet of ECF / ECOD / Pfam coverage of reprseqs, and coverage by PHROG functional class, for RefSeq and for PhageScope | `figure3.png` (panels A-D) | suppl-3p1 |
| `part4_arch-regimes/Fig4.R` | Architecture regimes: EF diversity vs number of ECF hits, diversity vs spatial segregation of EFs, example proteins, function-labelled diversity / segregation summary; also average hit overlap (AHO) by function, correlation structure and PCA of the architecture metrics | `figure4.png` (panels A-D) | suppl-4p1 to 4p4 |
| `part5_funct-rescue/Fig5.R` | Functional rescue: functional enrichment of EFs in domain-dark proteins (scatter and dot plot) | `figure5.png` | none |

`part4_arch-regimes/Fig4.R` and the part 6 scripts also read the raw reprseq lengths file.

### Part 6: RBP structures (`part6_rbp-struct`)

A four-step chain for the receptor-binding-protein (tail fibre / tail spike) structural figure:

- **`make1-shortlist.R`** filters tail-fibre / tail-spike reprseqs by EF count, length and ECF hits, plots their ECF / EF architectures (`vis-rbp-archs/`), writes the FASTA of the set chosen for AlphaFold (`af-model-set.fasta`), and saves a combined supplementary figure (`suppl-6p1`).
- **`make2-plot.R`** draws the EF / ECOD / Pfam architecture of five focal RBPs and writes `figure6A.pdf`, `ef-coords.txt` and `pymol-domain-colors.csv` (per-residue colour table) to `paper-ecf/part6_rbp-struct/`.
- **`make3-pymol.sh`** renders AlphaFold 3 models (`fold_<reprseq>.cif`) in PyMOL, colouring ECF / ECOD / Pfam segments from `pymol-domain-colors.csv`; it supports per-segment or combined images, monomer or trimer views, and a highlight colour.
- **`make4-run.sh`** is a wrapper that calls `make3-pymol.sh` with the options used for the trimer rendering with a pink highlight; the other runs are kept as commented lines.

The AlphaFold 3 models are not in the repository.

## Database (`scripts/4_prep-db`)

**`1_process-db.R`** reads `[raw]` the ECF A3M alignments only. It renames the first sequence of each alignment to the ECF name and writes `output/database/a3m/<ecf>.a3m`. It then builds an HH-suite database named `ecf_db` (a3m, hhm and cs219 ffindex sets sorted by length) in a scratch directory under `~/Data/tmp-ecfs`, copies it to `output/database/hhsuitedb/`, verifies the copy by md5 checksum, and runs a test `hhblits` search. It uses the HH-suite installation under `/opt/bio/hhsuite`.

**`2_ecf-metadata.R`** reads `ecf-table.txt`, `ecfs-family-table.txt`, `ecod-domains.txt`, `pfam-domains.txt`, `ecfs-vs-ecod.txt`, `ecfs-vs-pfam.txt`, `reprseq-ecfs.txt` and the reprseq function annotations. It writes `output/database/ecf-metadata.csv` with one row per ECF: EF, length, number of sequences, top confident ECOD T-level hit and name, top confident Pfam hit and description, and the top PHROG function(s) of the reprseqs carrying the ECF.

The two R scripts are independent of each other.

**`3_make-release.sh`** (bash) packs the database for release: `hhsuitedb/` and `a3m/` into `ecf-hhsuitedb.tar.gz` and `ecf-a3m.tar.gz`, a copy of `ecf-metadata.csv`, and `SHA256SUMS.txt`, written to `output/release/`. It checks the inputs first and verifies the archives and checksums afterwards.

**`4_push-release.sh VERSION`** (bash) creates the GitHub release `VERSION` as a draft on the current (pushed) commit, with those four files attached and the text of `release-notes/VERSION.md`. It verifies the checksums, refuses a tag that already exists (drafts included), and refuses notes that still contain `<TODO` placeholders. It requires the GitHub CLI (`gh`).

## Shared functions (`scripts/functions`)

- `helper.R`: interval overlap helpers, mapping of ECF amino-acid coordinates onto genomes, HMM-hit table helpers, and `plot.ecf.distr.over.query()`, which draws ECF hits along a query protein.
- `gpt-helper.R`: interval range reduction, covered-length and intersection-length helpers.
- `gggenomes.R`: preparation of gene, feature and link tables and window / zoom helpers for synteny plots with gggenomes.

## External dependencies

### What is not in the repository

- The HH-suite comparison tables and ECF alignments produced by colleagues' `ecf-explorer` run (`ecf-explorer_2026-10-09`): the reprseq-vs-reprseq HHR table; the protein-family name table, clustering table and reprseq lengths; the reprseq-vs-PHROG, -antidefence, -ECF, -ECOD and -Pfam tables; the ECF-vs-ECF, -ECOD and -Pfam tables; the PhageScope-vs-ECF, -ECOD and -Pfam tables; and one A3M alignment per ECF.
- The RefSeq input dataset: protein sequences, genome sequences, genome metadata and Prodigal protein statistics.
- **PHROGs** (annotation tables v4) and the manual PHROG mapping workbook, plus the antidefence sequence set and its description table.
- **ECOD** F70 (release 20230309, `develop288`) domain table and FASTA.
- **Pfam 37** (`pfamA.txt`, `pfamA.sql`).
- **ICTV** VMR (MSL39 v2).
- **PhageScope** random protein set (2,000 proteins: metadata and FASTA).
- The AlphaFold 3 models used by part 6.
- All generated data: intermediate tables, figures, the HH-suite database and `ecf-metadata.csv`.

### Code and tools

- **Functions repository.** `config.R` sources `basic.R` and `custom/hhr-cov.R` from a separate repository, expected at `~/Github/proj-rmostowy/functions`. The scripts call functions not defined in this repository, such as `make.dir`, `make.dirs`, `delete.dir`, `check.file`, `run.mcl`, `get.all.pairs`, `get.pair.table.hits`, `lw` and `num2s`; these are expected to come from there. If these files are missing, `safe_source()` in `config.R` only issues a warning and sourcing carries on.
- **Command-line tools** under `/opt/bio/<tool>/bin` (added to `PATH` by `config.R`): **MCL** (reprseq families and EFs; `/opt/bio/mcl/bin` is also passed explicitly), **MMseqs2** (`2_blast-reprseqs.R`), **BLAST+** (`makeblastdb`, `blastn`; `4_megablast-genomes.R`), **HH-suite** (`ffindex_*`, `hhmake`, `cstranslate`, `hhblits` and the `HHLIB` data files; `1_process-db.R`), **PyMOL** (`make3-pymol.sh`), and, for the REM analysis, **Clustal Omega** and **IQ-TREE 3**.
- **R packages.** The core packages are `data.table`, `dplyr`, `ggplot2`, `Biostrings`, `stringr`, `parallel`, `readxl`, `patchwork` and `rprojroot`. The figure scripts additionally use `ggh4x`, `igraph`, `ggraph`, `broom`, `ggExtra`, `ComplexUpset`, `ggrepel`, `ggnewscale` and `scales`; the synteny helpers use `gggenomes`. The REM analysis uses `DECIPHER`, `ape`, `phangorn`, `ggtree`, `ggtreeExtra`, `gggenes`, `ggalluvial` and `ggupset`. The repository contains no lockfile (`renv` or similar) pinning versions.
