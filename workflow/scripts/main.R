# ============================================================================
# PIMBA 3.0 - r_curation v5.0
# Main
# ============================================================================

# ---------------------------
# COMMAND LINE ARGUMENTS
# ---------------------------
args <- commandArgs(TRUE)

# Parameters:
# 1  = bold_ref_file
# 2  = rdp_ref_file
# 3  = silva_ref_file
# 4  = unite_ref_file
# 5 = greengenes_ref_file
# 6 = otu_file
# 7 = fasta_file
# 8 = hits_file
# 9 = ref_tax_file
# 10 = db_format
# 11 = mode

default_ref_files <- list(
  "COI-BOLD-DB"        = args[1],
  "16S-RDP-DB"         = args[2],
  "16S-SILVA-DB"       = args[3],
  "ITS-FUNGI-UNITE-DB" = args[4],
  "16S-GREENGENES-DB"  = args[5]
)

# ---------------------------
# PACKAGE DEPENDENCIES
# ---------------------------
required_packages <- c(
  "dplyr", "tidyr", "openxlsx", "readxl", "data.table",
  "stringr", "taxizedb", "ggplot2", "rentrez"
)

sapply(required_packages, require, character.only = TRUE)

library(Biostrings)

rm(required_packages)

# ---------------------------
# CONFIGURE LOCAL TAXIZEDB CACHE
# ---------------------------
ncbi_db_dir <- "/ncbi_db"

if (dir.exists(ncbi_db_dir)) {
  taxizedb::tdb_cache$cache_path_set(full_path = ncbi_db_dir)
  message("taxizedb cache path: ", taxizedb::tdb_cache$cache_path_get())
  message("NCBI DB path: ", taxizedb::db_path("ncbi"))
} else {
  message("Warning: /ncbi_db not found. taxizedb will use its default cache path.")
}

# ---------------------------
# LOAD CUSTOM FUNCTIONS
# ---------------------------
source("/app/custom_functions.R")

otu_file     <- args[6]
fasta_file   <- args[7]
hits_file    <- args[8]
ref_tax_file <- args[9]
db_format    <- args[10]
mode         <- args[11]

# ---------------------------
# INPUT ARGUMENTS
# ---------------------------
otu_file     <- args[6]
fasta_file   <- args[7]
hits_file    <- args[8]
#----------------------
ref_tax_file <- args[9]
if (!is.null(ref_tax_file) && ref_tax_file %in% c(" ", "", "NULL", "null", "NA", "na")) {
  ref_tax_file <- NULL
}
#----------------------
db_format    <- args[10]
mode         <- args[11]

# ---------------------------
# RUN PIPELINE
# ---------------------------
tax_datasets <- combine_tax_otu(
  otu_file = otu_file,
  fasta_file = fasta_file,
  hits_file = hits_file,
  ref_tax_file = ref_tax_file,
  db_format = db_format,
  mode = mode,
  default_ref_files = default_ref_files
)

parsing_taxa(
  tax_otu_df = tax_datasets$tax_otu_table,
  raw_otus_df = tax_datasets$raw_otus,
  requested_mode = mode
)
