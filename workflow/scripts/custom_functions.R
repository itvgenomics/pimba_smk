
find_ref_tax_file <- function(
  db_format,
  ref_tax_file = NULL,
  default_ref_files
) {
  
  # -----------------------------
  # NCBI-DB does not require a reference taxonomy file
  # -----------------------------
  if (db_format == "NCBI-DB") {
    return(NULL)
  }
  
  # -----------------------------
  # CUSTOM database uses user-provided reference taxonomy
  # parsed with BOLD-like formatting
  # -----------------------------
  if (db_format == "CUSTOM") {
    
    if (is.null(ref_tax_file) ||
        ref_tax_file == "" ||
        ref_tax_file %in% c("NA", "na", "NULL", "null") ||
        !file.exists(ref_tax_file)) {
      
      stop("For db_format = 'CUSTOM', a valid ref_tax_file must be provided.")
    }
    
    message("CUSTOM database selected.")
    message("Using user-provided reference taxonomy file: ", ref_tax_file)
    message("CUSTOM reference will be parsed using COI-BOLD-DB format.")
    
    return(ref_tax_file)
  }
  
  # -----------------------------
  # Validate supported PIMBA database names
  # -----------------------------
  if (!db_format %in% names(default_ref_files)) {
    stop("Unsupported db_format: ", db_format)
  }
  
  # -----------------------------
  # Retrieve developer-defined default reference taxonomy
  # -----------------------------
  ref_tax_file <- default_ref_files[[db_format]]
  
  if (is.null(ref_tax_file) ||
      ref_tax_file == "" ||
      ref_tax_file %in% c("NA", "na", "NULL", "null") ||
      !file.exists(ref_tax_file)) {
    
    stop(
      "Default reference taxonomy file not found for db_format = '",
      db_format,
      "': ",
      ref_tax_file
    )
  }
  
  message(
    "Using default reference taxonomy file for ",
    db_format,
    ": ",
    ref_tax_file
  )
  
  return(ref_tax_file)
}


combine_tax_otu <- function(otu_file, fasta_file,
                            hits_file, ref_tax_file = NULL,
                            db_format = c(
                            "COI-BOLD-DB",
                            "16S-RDP-DB",
                            "16S-GREENGENES-DB",
                            "16S-SILVA-DB",
                            "NCBI-DB",
                            "ITS-FUNGI-UNITE-DB",
                            "CUSTOM"
                            ),
                            mode = c("multi", "single"),
                            default_ref_files,
                            chunk_size = 50000) {
  
  mode <- match.arg(mode)
  db_format <- match.arg(db_format)
  
  # =========================================================
  # INTERNAL HELPER — selective reference-tax retrieval
  # =========================================================
  format_reference_tax_subset <- function(ref_tax_file,
                                          ref_ids,
                                          db_format,
                                          chunk_size = 50000) {
    
    ref_ids <- unique(ref_ids)
    ref_ids <- ref_ids[!is.na(ref_ids) & ref_ids != ""]
    
    if (length(ref_ids) == 0) {
      stop("No valid RefID values found in BLAST hits.")
    }
    
    # -----------------------------
    # NCBI-DB: direct parsing from RefIDs
    # -----------------------------
    if (db_format == "NCBI-DB") {
      
      message(
        "NCBI-DB mode detected: formatting taxonomy directly from matched RefIDs..."
      )
      
      tmp_df <- data.frame(
        V1 = ref_ids,
        stringsAsFactors = FALSE
      )
      
      ref_tax <- parse_ncbi(tmp_df)
      
      return(
        ref_tax %>%
          dplyr::distinct(RefID, .keep_all = TRUE)
      )
    }
    
    # -----------------------------
    # Other DBs: read reference file in chunks
    # -----------------------------
    if (is.null(ref_tax_file) || !file.exists(ref_tax_file)) {
      stop("Reference taxonomy file not found.")
    }
    
    message("Reading reference taxonomy in chunks...")
    message("Target matched RefIDs: ", length(ref_ids))
    
    con <- file(ref_tax_file, open = "r")
    on.exit(close(con), add = TRUE)
    
    matched_lines <- list()
    chunk_index <- 0
    total_lines <- 0
    total_matches <- 0
    
    repeat {
      
      lines <- readLines(con, n = chunk_size, warn = FALSE)
      
      if (length(lines) == 0) break
      
      chunk_index <- chunk_index + 1
      total_lines <- total_lines + length(lines)
      
      if (db_format == "ITS-FUNGI-UNITE-DB") {
        lines <- lines[!grepl("^Feature ID\tTaxon$", lines)]
      }
      
      if (length(lines) == 0) next
      
      first_col <- vapply(
        strsplit(lines, "\t", fixed = TRUE),
        function(x) if (length(x) >= 1) x[1] else NA_character_,
        character(1)
      )
      
      keep_idx <- which(first_col %in% ref_ids)
      
      if (length(keep_idx) > 0) {
        
        matched_lines[[length(matched_lines) + 1]] <- lines[keep_idx]
        
        total_matches <- total_matches + length(keep_idx)
      }
    }
    
    message("Total reference lines scanned: ", total_lines)
    message("Total matched reference lines recovered: ", total_matches)
    
    if (length(matched_lines) == 0) {
      
      warning("No matching RefID found in reference taxonomy file.")
      
      return(
        data.frame(
          RefID = character(),
          Taxonomy = character(),
          stringsAsFactors = FALSE
        )
      )
    }
    
    matched_lines <- unlist(matched_lines, use.names = FALSE)
    
    ref_tax_raw <- data.table::fread(
      text = matched_lines,
      sep = "\t",
      header = FALSE,
      quote = "",
      fill = TRUE,
      data.table = FALSE
    )
    
    # -----------------------------
    # Apply parser according to PIMBA database name
    # -----------------------------
    if (db_format == "COI-BOLD-DB" || db_format == "CUSTOM") {
      
      ref_tax <- parse_bold(ref_tax_raw)
      
    } else if (db_format %in% c("16S-RDP-DB", "16S-GREENGENES-DB")) {
      
      ref_tax_raw <- ref_tax_raw[, 1:2, drop = FALSE]
      colnames(ref_tax_raw) <- c("RefID", "Taxonomy")
      ref_tax <- ref_tax_raw
      
    } else if (db_format == "16S-SILVA-DB") {
      
      ref_tax <- parse_silva(ref_tax_raw)
      
    } else if (db_format == "ITS-FUNGI-UNITE-DB") {
      
      ref_tax <- parse_unite(ref_tax_raw)
      
    } else {
      
      stop("Unsupported db_format: ", db_format)
    }
    
    ref_tax <- ref_tax %>%
      dplyr::distinct(RefID, .keep_all = TRUE)
    
    return(ref_tax)
  }
  
  # =========================================================
  # RESOLVE REFERENCE FILE
  # =========================================================
  ref_tax_file <- find_ref_tax_file(
    db_format = db_format,
    ref_tax_file = ref_tax_file,
    default_ref_files = default_ref_files
  )
  
  # -----------------------------
  # STEP 1 — OTU / ASV table
  # -----------------------------
  if (!file.exists(otu_file)) {
    stop("OTU/ASV table not found.")
  }
  
  OTU <- data.table::fread(
    otu_file,
    sep = "\t",
    header = TRUE
  )
  
  colnames(OTU)[1] <- "UnitID"
  
  # -----------------------------
  # STEP 2 — FASTA
  # -----------------------------
  if (!file.exists(fasta_file)) {
    stop("FASTA file not found.")
  }
  
  fasta_seq <- Biostrings::readDNAStringSet(fasta_file)
  
  raw_units <- data.frame(
    UnitID = names(fasta_seq),
    Fasta = as.character(fasta_seq),
    stringsAsFactors = FALSE
  )
  
  # -----------------------------
  # STEP 3 — Detect mode: ASV or OTU
  # -----------------------------
  prefix <- substr(OTU$UnitID[1], 1, 3)
  
  detected_mode <- ifelse(
    toupper(prefix) == "ASV",
    "ASV",
    "OTU"
  )
  
  message("Detected mode: ", detected_mode)
  
  # -----------------------------
  # STEP 4 — Validate input files
  # -----------------------------
  if (missing(hits_file) || !file.exists(hits_file)) {
    stop("Hits file not found.")
  }
  
  if (db_format != "NCBI-DB" &&
      (is.null(ref_tax_file) || !file.exists(ref_tax_file))) {
    
    stop("Reference taxonomy file not found.")
  }
  
  # -----------------------------
  # STEP 5 — Read BLAST hits
  # -----------------------------
  message("Reading BLAST hits file...")
  
  hits <- data.table::fread(
    hits_file,
    header = FALSE,
    sep = "\t"
  )
  
  # -----------------------------
  # Normalize UnitID from newer BOLD logs
  # Old format:
  #   OTU12345
  # New format:
  #   gb|OTU12345|
  # Convert both to:
  #   OTU12345
  # -----------------------------
  hits[[1]] <- sub("^.*\\|(OTU[0-9]+)\\|?$", "\\1", hits[[1]])
  
  if (ncol(hits) >= 8) {
    
    colnames(hits)[1:8] <- c(
      "UnitID",
      "X1",
      "RefID",
      "X2",
      "Description",
      "PID",
      "Coverage",
      "Evalue"
    )
    
  } else if (ncol(hits) >= 7) {
    
    colnames(hits)[1:7] <- c(
      "UnitID",
      "X1",
      "RefID",
      "X2",
      "X3",
      "PID",
      "Coverage"
    )
    
  } else {
    
    colnames(hits)[1:3] <- c(
      "UnitID",
      "X1",
      "RefID"
    )
    
    if (ncol(hits) >= 6) {
      colnames(hits)[6] <- "PID"
    }
    
    if (ncol(hits) >= 7) {
      colnames(hits)[7] <- "Coverage"
    }
  }
  
  message("Total raw hits: ", nrow(hits))
  
  # -----------------------------
  # PID
  # -----------------------------
  if ("PID" %in% colnames(hits)) {
    
    hits$PID <- round(
      as.numeric(hits$PID) / 100,
      3
    )
    
  } else {
    
    hits$PID <- NA_real_
    
    message("Warning: PID column not found")
  }
  
  # -----------------------------
  # Coverage
  # -----------------------------
  if ("Coverage" %in% colnames(hits)) {
    
    suppressWarnings({
      hits$Coverage <- as.numeric(hits$Coverage)
    })
    
    hits$Coverage <- round(
      as.numeric(hits$Coverage) / 100,
      3
    )
    
    hits$Coverage <- ifelse(
      is.na(hits$Coverage),
      "Not recovered by BLAST",
      as.character(hits$Coverage)
    )
    
  } else {
    
    hits$Coverage <- "Not recovered by BLAST"
  }
  
  # -----------------------------
  # Clean RefIDs
  # -----------------------------
  hits <- hits %>%
    dplyr::mutate(
      RefID = ifelse(
        grepl("\\|", RefID),
        stringr::str_replace(
          RefID,
          ".*\\|(.*)\\|.*",
          "\\1"
        ),
        RefID
      )
    )
  
  # -----------------------------
  # Keep best hit per UnitID and RefID
  # -----------------------------
  hits <- hits %>%
    dplyr::group_by(UnitID, RefID) %>%
    dplyr::slice_max(
      PID,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
  
  # -----------------------------
  # STEP 6 — Selective taxonomy recovery
  # -----------------------------
  message("\nRecovering only matched reference taxonomy entries...")
  
  matched_ref_ids <- unique(hits$RefID)
  
  ref_tax <- format_reference_tax_subset(
    ref_tax_file = ref_tax_file,
    ref_ids = matched_ref_ids,
    db_format = db_format,
    chunk_size = chunk_size
  )
  
  # -----------------------------
  # STEP 7 — Process taxonomy ranks
  # -----------------------------
  ref_tax <- ref_tax %>%
    tidyr::separate(
      Taxonomy,
      into = c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      sep = ";",
      fill = "right",
      extra = "drop"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        Kingdom:Species,
        ~ stringr::str_trim(
          stringr::str_replace(., "^.*__", "")
        )
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        Kingdom:Species,
        ~ ifelse(is.na(.) | . == "", "", .)
      )
    ) %>%
    dplyr::distinct(RefID, .keep_all = TRUE)
  
  # -----------------------------
  # STEP 8 — Build taxonomic table
  # -----------------------------
  message("\nJoining hits with taxonomy...")
  
  hits_ids <- unique(hits$RefID)
  ref_ids <- unique(ref_tax$RefID)
  
  matches <- sum(hits_ids %in% ref_ids)
  
  message("Unique IDs in hits: ", length(hits_ids))
  message("Unique IDs recovered in formatted taxonomy: ", length(ref_ids))
  message("Matches found: ", matches)
  
  Tax <- hits %>%
    dplyr::left_join(ref_tax, by = "RefID") %>%
    dplyr::select(
      UnitID,
      RefID,
      PID,
      Coverage,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      Species
    )
  
  # -----------------------------
  # STEP 9 — Single / multi-hit behavior
  # -----------------------------
  if (mode == "single") {
    
    message(
      "Applying SINGLE mode: keeping only best hit per UnitID"
    )
    
    Tax <- Tax %>%
      dplyr::group_by(UnitID) %>%
      dplyr::slice_max(
        PID,
        n = 1,
        with_ties = FALSE
      ) %>%
      dplyr::ungroup()
    
  } else {
    
    message(
      "Applying MULTI mode: keeping all hits per UnitID"
    )
  }
  
  # -----------------------------
  # STEP 10 — Merge taxonomy with OTU / ASV table
  # -----------------------------
  message("\nMerging with ", prefix, " table...")
  
  OTU_expanded <- OTU[
    match(Tax$UnitID, OTU$UnitID),
  ]
  
  combined_data <- cbind(
    Tax,
    OTU_expanded[, -1, drop = FALSE]
  )
  
  # -----------------------------
  # STEP 11 — Return datasets
  # -----------------------------
  message("\nProcessing completed successfully!")
  
  message(
    "Final table dimensions: ",
    nrow(combined_data),
    " rows x ",
    ncol(combined_data),
    " columns"
  )
  
  return(
    list(
      mode = detected_mode,
      tax_otu_table = combined_data,
      raw_otus = raw_units
    )
  )
}

parsing_taxa <- function(
  tax_otu_df,
  raw_otus_df,
  output_xlsx = "tax_assignments__to_validate.xlsx",
  requested_mode = NULL
) {
  
  # -----------------------------
  # Detect mode: OTU or ASV
  # -----------------------------
  first_id <- if (
    !is.null(tax_otu_df[[1]]) &&
    length(tax_otu_df[[1]]) > 0 &&
    !is.na(tax_otu_df[[1]][1])
  ) {
    as.character(tax_otu_df[[1]][1])
  } else if (
    !is.null(raw_otus_df[[1]]) &&
    length(raw_otus_df[[1]]) > 0 &&
    !is.na(raw_otus_df[[1]][1])
  ) {
    as.character(raw_otus_df[[1]][1])
  } else {
    ""
  }
  
  prefix <- if (nzchar(first_id)) substr(first_id, 1, 3) else ""
  mode <- ifelse(toupper(prefix) == "ASV", "ASV", "OTU")
  
  message("Detected mode: ", mode)
  
  id_column <- ifelse(mode == "ASV", "ASVs", "OTUs")
  
  names(tax_otu_df)[1] <- id_column
  names(raw_otus_df)[1] <- id_column
  
  # -----------------------------
  # Detect single or multi taxa
  # -----------------------------
  has_refid <- "RefID" %in% colnames(tax_otu_df)
  has_multiple <- FALSE
  
if (has_refid) {
  
  id_counts <- tax_otu_df %>%
    dplyr::group_by(!!rlang::sym(id_column)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      .groups = "drop"
    )
  
  has_multiple <- any(id_counts$n > 1)
  
  if (has_multiple) {
    
    message("Multi Taxa per ", mode)
    
  } else {
    
    message("Unique Taxa per ", mode)
    
    if (!is.null(requested_mode) &&
        requested_mode == "multi") {
      
      message(
    "Multi-hit output was requested, but only single-hit taxonomic assignments were detected. \n",
    "Generating singleSeq output instead.\n",
    "This commonly occurs when hits_per_subject was previously set to 1 by the user in earlier PIMBA workflow steps."
    )
    }
  }
  
} else {
  
  message("Unique Taxa per ", mode)
}
  
  suffix <- ifelse(has_multiple, "_multiSeq", "_singleSeq")
  
  output_xlsx <- sub(
    "\\.xlsx$",
    paste0(suffix, ".xlsx"),
    output_xlsx
  )
  
  # -----------------------------
  # Merge taxonomic table with sequences
  # -----------------------------
  raw_data_with_seqs <- merge(
    tax_otu_df,
    raw_otus_df,
    by = id_column,
    all.x = TRUE
  ) %>%
    dplyr::relocate(Fasta, .after = all_of(id_column)) %>%
    dplyr::mutate(Length = nchar(Fasta)) %>%
    dplyr::relocate(Length, .after = Fasta)
  
  # -----------------------------
  # Initial species/genus preparation
  # -----------------------------
  filtered_data_with_flags <<- raw_data_with_seqs %>%
    dplyr::mutate(
      Species = stringr::str_replace_all(Species, "_", " "),
      Genus   = stringr::str_replace_all(Genus, "_", " ")
    ) %>%
    dplyr::mutate(
      Species = mapply(
        function(sp, gn) {
          if (!is.na(sp) && !is.na(gn) && gn != "") {
            
            gn_words <- unlist(stringr::str_split(gn, "\\s+"))
            
            for (w in gn_words) {
              if (nzchar(w)) {
                sp <- stringr::str_remove_all(
                  sp,
                  stringr::regex(w, ignore_case = TRUE)
                )
              }
            }
            
            stringr::str_trim(stringr::str_squish(sp))
            
          } else {
            sp
          }
        },
        Species,
        Genus,
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
    )
  
  # -----------------------------
  # Species cleaning and validation flags
  # -----------------------------
  filtered_data_with_flags <<- filtered_data_with_flags %>%
  # Remove placeholder species "None"
  # when genus is already assigned
  dplyr::mutate(
    Species = ifelse(
      !is.na(Genus) &
        Genus != "" &
        !is.na(Species) &
        Species == "None",
      "",
      Species
    )
  ) %>%
    dplyr::mutate(
      dplyr::across(
        c(Species),
        ~ ifelse(grepl("uncultured", ., ignore.case = TRUE), "", .)
      )
    ) %>%
    dplyr::mutate(
      Species = ifelse(
        grepl("sp\\.", Species, ignore.case = TRUE) &
          (is.na(Genus) | Genus == ""),
        "",
        Species
      )
    ) %>%
    dplyr::mutate(
      Species = ifelse(
        !is.na(Species) &
          Species != "" &
          (is.na(Genus) | Genus == ""),
        "",
        Species
      )
    ) %>%
    dplyr::mutate(
      Species = ifelse(
        grepl("^sp\\.\\s*[0-9]+$", Species, ignore.case = TRUE),
        "",
        Species
      )
    ) %>%
    dplyr::mutate(
      Species = ifelse(
        stringr::str_detect(
          Species,
          "[A-Z]{2}|[A-Z][0-9]|[0-9][A-Z]|[0-9]{2}"
        ),
        stringr::str_replace_all(
          Species,
          "[^a-záéíóúâêîôûãõç. ]",
          " "
        ),
        Species
      )
    ) %>%
    dplyr::mutate(
      Species = stringr::str_replace_all(Species, "\\s+", " ")
    ) %>%
    dplyr::mutate(
      Species = ifelse(
        stringr::str_detect(Species, "^sp\\.[:alnum:.]$") |
          stringr::str_detect(Species, "^sp\\..*[^[:alnum:][:space:]]"),
        "sp.",
        Species
      )
    ) %>%
    dplyr::mutate(
      Species = stringr::str_replace_all(
        Species,
        "\\b(sp\\.?)(\\s|$)",
        ""
      )
    ) %>%
    dplyr::mutate(
      Species = stringr::str_trim(Species)
    ) %>%
    dplyr::mutate(
      Species = ifelse(nchar(Species) == 1, "", Species)
    ) %>%
    dplyr::mutate(
      Validation = dplyr::case_when(
  
  is.na(Species) | Species == "" ~ "Empty Species Field",
  
  stringr::str_count(Species, "\\S+") > 1 |
    stringr::str_detect(Species, "[^a-záéíóúâêîôûãõç]") ~ "Format Issue",
  
  TRUE ~ "Format Validated"
)
    ) %>%
    dplyr::relocate(Validation, .after = 1)
  
  # -----------------------------
  # NCBI validation section
  # -----------------------------
  validated_data <- filtered_data_with_flags %>%
    dplyr::filter(
      Validation == "Format Validated",
      !is.na(Species),
      Species != ""
    ) %>%
    dplyr::mutate(
      Full_Name = stringr::str_trim(
        paste(Genus, Species, sep = " ")
      )
    ) %>%
    dplyr::group_by(Full_Name) %>%
    dplyr::ungroup()
  
  if (nrow(validated_data) > 0) {
    
    message("Using taxizedb cache path: ", taxizedb::db_path("ncbi"))
    
    src <- taxizedb::src_ncbi()
    
    taxa_names_vec <- unique(validated_data$Full_Name)
    
taxa_names <- lapply(taxa_names_vec, function(x) {
  tryCatch(
    {
      taxizedb::classification(x, db = "ncbi")[[1]]
    },
    error = function(e) {
      message("NCBI validation failed for ", x, ": ", e$message)
      
      data.frame(
        name = NA,
        rank = NA,
        id = NA
      )
    }
  )
})
    

    names(taxa_names) <- taxa_names_vec
    
    is_dataframe <- function(obj) {
      inherits(obj, "data.frame")
    }
    
    for (i in seq_along(taxa_names)) {
      if (!is_dataframe(taxa_names[[i]])) {
        taxa_names[[i]] <- data.frame(
          name = NA,
          rank = NA,
          id = NA
        )
      }
    }
    
    taxa_names_df <- do.call(
      rbind,
      lapply(names(taxa_names), function(name) {
        df <- taxa_names[[name]]
        df$taxid <- name
        return(df)
      })
    )
    
    taxa_names_df <- taxa_names_df[
      ,
      c("taxid", names(taxa_names_df)[1:(ncol(taxa_names_df) - 1)])
    ]
    
    taxa_names_df <- taxa_names_df[, -ncol(taxa_names_df)]
    
    taxa_names_df <- taxa_names_df %>%
      dplyr::filter(
        rank %in% c(
          "kingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species",
          NA
        )
      ) %>%
      dplyr::group_by(taxid) %>%
      dplyr::distinct()
    
    identified_otus <- taxa_names_df %>%
      tidyr::pivot_wider(
        names_from = rank,
        values_from = name
      )
    
    if (!"species" %in% names(identified_otus)) {
      identified_otus$species <- NA_character_
    }
    
    names(identified_otus)[names(identified_otus) == "species"] <- "NCBI_species"
    
    # -----------------------------
    # Merge NCBI results
    # -----------------------------
    if (has_refid) {
      
      validated_data_updated <- validated_data %>%
        dplyr::left_join(
          identified_otus %>%
            dplyr::select(taxid, NCBI_species),
          by = c("Full_Name" = "taxid")
        ) %>%
        dplyr::mutate(
          Validation_new = ifelse(
            is.na(NCBI_species),
            Validation,
            "Format and DB Validated"
          )
        ) %>%
        dplyr::select(
          dplyr::all_of(id_column),
          RefID,
          Validation_new
        ) %>%
        dplyr::distinct()
      
    } else {
      
      validated_data_updated <- validated_data %>%
        dplyr::left_join(
          identified_otus %>%
            dplyr::select(taxid, NCBI_species),
          by = c("Full_Name" = "taxid")
        ) %>%
        dplyr::mutate(
          Validation_new = ifelse(
            is.na(NCBI_species),
            Validation,
            "Format and DB Validated"
          )
        ) %>%
        dplyr::select(
          dplyr::all_of(id_column),
          Validation_new
        ) %>%
        dplyr::distinct()
    }
    
    # -----------------------------
    # Update main dataset
    # -----------------------------
    if (has_refid) {
      
      filtered_data_with_flags <<- filtered_data_with_flags %>%
        dplyr::left_join(
          validated_data_updated,
          by = c(id_column, "RefID")
        ) %>%
        dplyr::mutate(
          Validation = dplyr::coalesce(
            Validation_new,
            Validation
          )
        ) %>%
        dplyr::select(-Validation_new)
      
    } else {
      
      filtered_data_with_flags <<- filtered_data_with_flags %>%
        dplyr::left_join(
          validated_data_updated,
          by = id_column
        ) %>%
        dplyr::mutate(
          Validation = dplyr::coalesce(
            Validation_new,
            Validation
          )
        ) %>%
        dplyr::select(-Validation_new)
    }
  }
  
  # -----------------------------
  # Export results
  # -----------------------------
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "raw_data_with_seqs")
  openxlsx::addWorksheet(wb, "filtered_data_with_flags")
  
  openxlsx::writeData(
    wb,
    sheet = "raw_data_with_seqs",
    raw_data_with_seqs
  )
  
  openxlsx::writeData(
    wb,
    sheet = "filtered_data_with_flags",
    filtered_data_with_flags
  )
  
  dir.create(
    dirname(output_xlsx),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  openxlsx::saveWorkbook(
    wb,
    output_xlsx,
    overwrite = TRUE
  )
  
  message("File saved as: ", output_xlsx)
}

format_reference_tax <- function(
  ref_tax_file,
  db_format = c(
    "COI-BOLD-DB",
    "16S-RDP-DB",
    "16S-GREENGENES-DB",
    "16S-SILVA-DB",
    "NCBI-DB",
    "ITS-FUNGI-UNITE-DB",
    "CUSTOM"
  )
) {
  
  db_format <- match.arg(db_format)
  
  if (!file.exists(ref_tax_file)) {
    stop("Reference taxonomy file not found.")
  }
  
  if (db_format == "CUSTOM") {
    
    message("CUSTOM database selected.")
    message("CUSTOM reference will be parsed using COI-BOLD-DB format.")
  }
  
  input_dir <- dirname(ref_tax_file)
  base_name <- tools::file_path_sans_ext(basename(ref_tax_file))
  
  output_file <- file.path(
    input_dir,
    paste0(base_name, "_for_curation.txt")
  )
  
  ref_tax <- data.table::fread(
    ref_tax_file,
    sep = "\t",
    header = FALSE
  )
  
  # ----------------------------------------------------------
  # Apply parser according to PIMBA database naming
  # ----------------------------------------------------------
  if (db_format == "COI-BOLD-DB" || db_format == "CUSTOM") {
    
    ref_tax <- parse_bold(ref_tax)
    
  } else if (db_format %in% c("16S-RDP-DB", "16S-GREENGENES-DB")) {
    
    ref_tax <- ref_tax[, 1:2]
    colnames(ref_tax) <- c("RefID", "Taxonomy")
    
  } else if (db_format == "16S-SILVA-DB") {
    
    ref_tax <- parse_silva(ref_tax)
    
  } else if (db_format == "NCBI-DB") {
    
    ref_tax <- parse_ncbi(ref_tax)
    
  } else if (db_format == "ITS-FUNGI-UNITE-DB") {
    
    ref_tax <- parse_unite(ref_tax)
  }
  
  ref_tax <- ref_tax %>%
    dplyr::distinct(RefID, .keep_all = TRUE)
  
  data.table::fwrite(
    ref_tax,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  message("Formatted reference taxonomy saved as: ", output_file)
  
  invisible(ref_tax)
}

parse_bold <- function(df) {
  
  # Assume:
  # col1 = RefID
  # col2 = taxonomy string
  
  df <- df[, 1:2]
  colnames(df) <- c("RefID", "Taxonomy")
  
  df <- df %>%
    tidyr::separate(
      Taxonomy,
      into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species_raw"),
      sep = ";",
      fill = "right"
    )
  
  # Clean
  df <- df %>%
    dplyr::mutate(
      dplyr::across(Kingdom:Genus, ~ stringr::str_trim(.)),
      Species_raw = stringr::str_trim(Species_raw)
    )
  
  # Extract species epithet
  df <- df %>%
    dplyr::mutate(
      Species = dplyr::case_when(
        is.na(Species_raw) ~ "",
        Species_raw == "" ~ "",
        TRUE ~ stringr::word(Species_raw, -1)
      )
    )
  
  # Replace NA
  df <- df %>%
    dplyr::mutate(
      dplyr::across(Kingdom:Species, ~ ifelse(is.na(.) | . == "", "", .))
    )
  
  # Build final taxonomy string
  df <- df %>%
    dplyr::mutate(
      Taxonomy = paste0(
        "k__", Kingdom, "; ",
        "p__", Phylum, "; ",
        "c__", Class, "; ",
        "o__", Order, "; ",
        "f__", Family, "; ",
        "g__", Genus, "; ",
        "s__", Species
      )
    ) %>%
    dplyr::select(RefID, Taxonomy)
  
  return(df)
}

parse_silva <- function(df) {
  
  # Assume: col1 = RefID, col2 = taxonomy string
  df <- df[, 1:2]
  colnames(df) <- c("RefID", "Taxonomy")
  
  # SEPARAR DINAMICAMENTE (máximo de colunas possível)
  df_sep <- df %>%
    tidyr::separate(
      Taxonomy,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Extra1", "Extra2"),
      sep = ";",
      fill = "right",
      extra = "merge"  # IMPORTANTE: junta o que sobra na última coluna
    )
  
  # Limpar espaços
  df_sep <- df_sep %>%
    dplyr::mutate(
      dplyr::across(Kingdom:Extra2, ~ stringr::str_trim(.))
    )
  
  # Lidar com os níveis extras do SILVA
  df_sep <- df_sep %>%
    dplyr::mutate(
      # Se tem Extra1, significa que Species está incompleto
      Species = dplyr::case_when(
        !is.na(Extra1) & Extra1 != "" ~ paste(Species, Extra1, Extra2, sep = " "),
        !is.na(Extra1) ~ Species,
        TRUE ~ Species
      ),
      # Limpar espaços múltiplos
      Species = stringr::str_squish(Species)
    ) %>%
    dplyr::select(-Extra1, -Extra2)
  
  # Substituir NA por vazio
  df_sep <- df_sep %>%
    dplyr::mutate(
      dplyr::across(Kingdom:Species, ~ ifelse(is.na(.) | . == "", "", .))
    )
  
  # Extrair epíteto específico (última palavra)
  df_sep <- df_sep %>%
    dplyr::mutate(
      Species_epithet = dplyr::case_when(
        Species == "" ~ "",
        TRUE ~ stringr::word(Species, -1)
      )
    )
  
  # Construir taxonomy no formato padrão (7 níveis)
  df_sep <- df_sep %>%
    dplyr::mutate(
      Taxonomy = paste0(
        "k__", Kingdom, "; ",
        "p__", Phylum, "; ",
        "c__", Class, "; ",
        "o__", Order, "; ",
        "f__", Family, "; ",
        "g__", Genus, "; ",
        "s__", Species_epithet
      )
    ) %>%
    dplyr::select(RefID, Taxonomy)
  
  return(df_sep)
}

parse_ncbi <- function(df) {
  
  # Assume: col1 = RefID (pode ser accession number ou string completa do BLAST)
  # Exemplo: "gi|699050543|gb|KF010490.1|" ou apenas "KF010490.1"
  
  message("Parsing NCBI taxonomy from accessions...")
  
  # Extrair apenas a primeira coluna (accessions)
  df <- df[, 1, drop = FALSE]
  colnames(df) <- c("RawID")
  
  # Extrair accession numbers usando a mesma lógica do refDB_Blast
  df <- df %>%
    dplyr::mutate(
      Accession = dplyr::case_when(
        # Formato gi|699050543|gb|KF010490.1|
        grepl("\\|", RawID) ~ stringr::str_extract(RawID, "[A-Z]{2}[0-9]+\\.?[0-9]*"),
        # Se já for accession direto
        grepl("^[A-Z]{2}[0-9]+\\.?[0-9]*$", RawID) ~ RawID,
        # Se for número puro (ex: NCBI taxid)
        grepl("^[0-9]+$", RawID) ~ RawID,
        # Caso contrário, manter original
        TRUE ~ RawID
      )
    ) %>%
    dplyr::filter(!is.na(Accession) & Accession != "")
  
  # Remover duplicatas
  df <- dplyr::distinct(df, Accession, .keep_all = TRUE)
  
  message("Fetching NCBI taxonomy for ", nrow(df), " unique accessions...")
  
  # Função para obter taxid do accession
  get_taxid <- function(accession) {
    tryCatch({
      # Se já for um taxid numérico, usar diretamente
      if (grepl("^[0-9]+$", accession)) {
        return(accession)
      }
      
      # Caso contrário, buscar taxid no nuccore
      summary <- rentrez::entrez_summary(
        db = "nuccore",
        id = accession
      )
      
      return(summary$taxid)
      
    }, error = function(e) {
      message("Error fetching taxid for ", accession, ": ", e$message)
      return(NA)
    })
  }
  
  # Fallback online: buscar classificação diretamente no NCBI Taxonomy
  get_taxonomy_from_ncbi_online <- function(taxid) {
    
    tryCatch({
      
      xml <- rentrez::entrez_fetch(
        db = "taxonomy",
        id = taxid,
        rettype = "xml",
        parsed = FALSE
      )
      
      doc <- xml2::read_xml(xml)
      
      taxon_nodes <- xml2::xml_find_all(
        doc,
        ".//LineageEx/Taxon"
      )
      
      lineage_df <- data.frame(
        name = xml2::xml_text(
          xml2::xml_find_first(taxon_nodes, "ScientificName")
        ),
        rank = xml2::xml_text(
          xml2::xml_find_first(taxon_nodes, "Rank")
        ),
        id = xml2::xml_text(
          xml2::xml_find_first(taxon_nodes, "TaxId")
        ),
        stringsAsFactors = FALSE
      )
      
      current_name <- xml2::xml_text(
        xml2::xml_find_first(doc, ".//Taxon/ScientificName")
      )
      
      current_rank <- xml2::xml_text(
        xml2::xml_find_first(doc, ".//Taxon/Rank")
      )
      
      current_id <- xml2::xml_text(
        xml2::xml_find_first(doc, ".//Taxon/TaxId")
      )
      
      tax_data <- dplyr::bind_rows(
        lineage_df,
        data.frame(
          name = current_name,
          rank = current_rank,
          id = current_id,
          stringsAsFactors = FALSE
        )
      )
      
      tax_data <- tax_data %>%
        dplyr::filter(
          !is.na(name),
          !is.na(rank),
          name != "",
          rank != ""
        )
      
      return(tax_data)
      
    }, error = function(e) {
      message("Online NCBI taxonomy fallback failed for taxid ", taxid, ": ", e$message)
      return(data.frame(name = NA, rank = NA, id = NA))
    })
  }
  
  # Obter taxids para cada accession
  taxid_list <- sapply(df$Accession, get_taxid)
  
  # Criar dataframe com Accession e taxid
  accession_taxid <- data.frame(
    Accession = df$Accession,
    taxid = as.character(taxid_list),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(taxid))
  
  if (nrow(accession_taxid) == 0) {
    stop("No valid taxids found")
  }
  
  message("Getting taxonomic classification for ", nrow(accession_taxid), " taxids...")
  
  # Primeiro tenta classificação offline no NCBI.sql via taxizedb
  tax_classification <- taxizedb::classification(
    accession_taxid$taxid,
    db = "ncbi"
  )
  
  # Processar classificação
  taxonomy_list <- list()
  
  for (i in seq_along(tax_classification)) {
    
    accession <- accession_taxid$Accession[i]
    taxid <- accession_taxid$taxid[i]
    tax_data <- tax_classification[[i]]
    
    # Se taxizedb não encontrar, usar fallback online no NCBI Taxonomy
    if (!inherits(tax_data, "data.frame")) {
      
      message(
        "TaxID ",
        taxid,
        " not found in local NCBI.sql for ",
        accession,
        ". Trying online NCBI Taxonomy..."
      )
      
      tax_data <- get_taxonomy_from_ncbi_online(taxid)
    }
    
    # Garantir que é data.frame
    if (!inherits(tax_data, "data.frame")) {
      tax_data <- data.frame(name = NA, rank = NA, id = NA)
    }
    
    # Filtrar apenas ranks padrão
    tax_data <- tax_data %>%
      dplyr::filter(
        rank %in% c(
          "kingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species"
        )
      )
    
    if (nrow(tax_data) > 0) {
      
      # Criar vetor nomeado com os ranks
      tax_vector <- setNames(tax_data$name, tax_data$rank)
      
      taxonomy_list[[accession]] <- data.frame(
        RefID = accession,
        Kingdom = ifelse("kingdom" %in% names(tax_vector), tax_vector["kingdom"], ""),
        Phylum = ifelse("phylum" %in% names(tax_vector), tax_vector["phylum"], ""),
        Class = ifelse("class" %in% names(tax_vector), tax_vector["class"], ""),
        Order = ifelse("order" %in% names(tax_vector), tax_vector["order"], ""),
        Family = ifelse("family" %in% names(tax_vector), tax_vector["family"], ""),
        Genus = ifelse("genus" %in% names(tax_vector), tax_vector["genus"], ""),
        Species_raw = ifelse("species" %in% names(tax_vector), tax_vector["species"], ""),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combinar todas as taxonomias
  if (length(taxonomy_list) == 0) {
    stop("No taxonomy data retrieved")
  }
  
  ref_tax <- dplyr::bind_rows(taxonomy_list)
  
  # Extrair epíteto específico da espécie
  ref_tax <- ref_tax %>%
    dplyr::mutate(
      Species = dplyr::case_when(
        Species_raw == "" ~ "",
        TRUE ~ stringr::word(Species_raw, -1)
      )
    ) %>%
    dplyr::select(-Species_raw)
  
  # Substituir NA por vazio
  ref_tax <- ref_tax %>%
    dplyr::mutate(
      dplyr::across(
        Kingdom:Species,
        ~ ifelse(is.na(.) | . == "", "", .)
      )
    )
  
  # Construir taxonomy string no formato padrão
  ref_tax <- ref_tax %>%
    dplyr::mutate(
      Taxonomy = paste0(
        "k__", Kingdom, "; ",
        "p__", Phylum, "; ",
        "c__", Class, "; ",
        "o__", Order, "; ",
        "f__", Family, "; ",
        "g__", Genus, "; ",
        "s__", Species
      )
    ) %>%
    dplyr::select(RefID, Taxonomy)
  
  message("Successfully parsed ", nrow(ref_tax), " NCBI taxonomy entries")
  
  return(ref_tax)
}

parse_unite <- function(df) {
  
  # Assume:
  # col1 = RefID
  # col2 = taxonomy string
  
  df <- df[, 1:2]
  colnames(df) <- c("RefID", "Taxonomy")
  
  df <- df %>%
    dplyr::filter(!(RefID == "Feature ID" & Taxonomy == "Taxon"))
  
  df <- df %>%
    tidyr::separate(
      Taxonomy,
      into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species","SH"),
      sep = ";",
      fill = "right"
    )
  
  # Clean
  df <- df %>%
    dplyr::mutate(
      dplyr::across(Kingdom:Species, ~ stringr::str_trim(.))
    )
  
  # Replace NA
  df <- df %>%
    dplyr::mutate(
      dplyr::across(Kingdom:Species, ~ ifelse(is.na(.) | . == "", "", .))
    )
  
  # Build final taxonomy string
  df <- df %>%
    dplyr::mutate(
      Taxonomy = paste0(
        "k__", stringr::str_replace(Kingdom, "^.*__", ""), "; ",
        "p__", stringr::str_replace(Phylum, "^.*__", ""), "; ",
        "c__", stringr::str_replace(Class, "^.*__", ""), "; ",
        "o__", stringr::str_replace(Order, "^.*__", ""), "; ",
        "f__", stringr::str_replace(Family, "^.*__", ""), "; ",
        "g__", stringr::str_replace(Genus, "^.*__", ""), "; ",
        "s__", stringr::str_replace(Species, "^.*__", "")
      )
    ) %>%
    dplyr::select(RefID, Taxonomy)
  
  return(df)
}