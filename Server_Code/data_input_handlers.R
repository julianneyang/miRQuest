## Lines 15- 19 of this script are taken from https://github.com/nf-core/smrnaseq/blob/master/bin/edgeR_miRBase.r (MIT license).
## If you use this software, please also cite: 10.5281/zenodo.3456879 (nf-core/smrnaseq)
# Handler for "Use Default Data" button
observeEvent(input$useDefaultData, {
  # Load miRNA data
  mirna_path <- here::here("inst", "Data", "COREAD_Downsampled_Tumor_vs_Normal_Mature_miRNA.csv")
  data <- read.csv(mirna_path, header = TRUE)
  
  df_wide <- data
  df_wide[is.na(df_wide)] <- 0
  row.names(df_wide) <- df_wide$X
  df_wide <- df_wide %>% 
    dplyr::select(-c(X))

  row_sub = apply(df_wide, 1, function(row) all(row == 0))
  df_wide <- df_wide[!row_sub, , drop = FALSE]
  
  drop_colsum_zero <- (colSums(df_wide, na.rm = TRUE) != 0)
  df_wide <- df_wide[, drop_colsum_zero]
  
  wideData(df_wide)
  
  ld <- data %>%
    pivot_longer(cols = -1, names_to = "SampleID", values_to = "Count")
  
  longData(ld)
  output$mirna_status <- renderText("Demo data loaded (COREAD)")
  
  # Load metadata
  metadata_path <- here::here("inst", "Data", "COREAD_Tumor_vs_Normal_Metadata.csv")
  metadataData <- read.csv(metadata_path, header = TRUE)
  metadata(metadataData)
  updateSelectInput(session, "designColumn", choices = names(metadataData))
  updateSelectInput(session, "groupBy", choices = names(metadataData))
  updateSelectInput(session, "groupBy_Boxplot", choices = names(metadataData))
  updateSelectInput(session, "metadataColumn", choices = names(metadataData))
  output$metadata_status <- renderText("Demo data loaded (COREAD)")
  
  # Load mRNA data
  mrna_path <- here::here("inst", "Data", "COREAD_Downsampled_Tumor_vs_Normal_mRNA.csv")
  mrna_data <- read.csv(mrna_path, header = TRUE)
  mrna_wideData(mrna_data)
  output$mrna_status <- renderText("Demo data loaded (COREAD)")
  
  # Set species to human for COREAD data
  updateCheckboxInput(session, "hsa", value = TRUE)
  updateCheckboxInput(session, "mmu", value = FALSE)
  
  # Select tissue_type in all metadata dropdowns if it exists
  if ("tissue_type" %in% names(metadataData)) {
    updateSelectInput(session, "designColumn", selected = "tissue_type")
    updateSelectInput(session, "groupBy", selected = "tissue_type")
    updateSelectInput(session, "groupBy_Boxplot", selected = "tissue_type")
    updateSelectInput(session, "metadataColumn", selected = "tissue_type")
    
    # Click the plot button to show stacked bar chart
    shinyjs::delay(500, {
      shinyjs::click("plotButton")
    })
  }
  
  showNotification("Default COREAD data loaded successfully!", type = "message", duration = 3)
  
  # Show workflow guidance
  showNotification(
    HTML("Next steps:<br>
          1. Go to 'Differential miRNA Expression Analysis' tab and click 'Run DESeq2 Analysis'<br>
          2. For correlation analysis, also run 'Differential Gene Expression Analysis'"),
    type = "message",
    duration = 10
  )
})

# # Module 0A: Read and Store MiRNA Counts Data
# observeEvent(input$countTable, {
#   req(input$countTable)
#   
#   data <- read.csv(input$countTable$datapath, header = TRUE)
#   
#   if (is.null(data)) return()
#   
#   df_wide <- data
#   df_wide[is.na(df_wide)] <- 0
#   row.names(df_wide) <- df_wide$X
#   df_wide <- df_wide %>% 
#     dplyr::select(-c(X))
#   
#   row_sub = apply(df_wide, 1, function(row) all(row == 0))
#   df_wide <- df_wide[!row_sub, , drop = FALSE]
#   
#   drop_colsum_zero <- (colSums(df_wide, na.rm = TRUE) != 0)
#   df_wide <- df_wide[, drop_colsum_zero]
#   
#   wideData(df_wide)
#   
#   ld <- data %>%
#     pivot_longer(cols = -1, names_to = "SampleID", values_to = "Count")
#   
#   longData(ld)
#   
#   # Update status to show user data is loaded
#   output$mirna_status <- renderText("User data loaded")
# })
# 
# # Module 0B: Read and Store Metadata
# observeEvent(input$metadataFile, {
#   req(input$metadataFile)
#   
#   metadataData <- read.csv(input$metadataFile$datapath, header = TRUE)
#   
#   metadata(metadataData)
#   
#   updateSelectInput(session, "designColumn", choices = names(metadataData))
#   updateSelectInput(session, "groupBy", choices = names(metadataData))
#   updateSelectInput(session, "groupBy_Boxplot", choices = names(metadataData))
#   updateSelectInput(session, "metadataColumn", choices = names(metadataData))
#   
#   # Update status to show user data is loaded
#   output$metadata_status <- renderText("User data loaded")
# })
# 
# # Module 0C: Read and Store mRNA Data
# observeEvent(input$mrna_countTable, {
#   req(input$mrna_countTable)
#   
#   data <- read.csv(input$mrna_countTable$datapath, header = TRUE)
#   
#   if (is.null(data)) return()
#   
#   mrna_wideData(data)
#   
#   # Update status to show user data is loaded
#   output$mrna_status <- renderText("User data loaded")
# })

# Module 0A: Read and Store miRNA Counts Data
observeEvent(input$countTable, {
  req(input$countTable)
  cat("DEBUG: entered miRNA upload observer\n")
  
  data <- tryCatch(
    read.csv(input$countTable$datapath,
             header = TRUE),
    error = function(e) {
      showNotification(
        paste("Error reading miRNA count file:", e$message),
        type = "error", duration = NULL
      )
      return(NULL)
    }
  )
  if (is.null(data)) {
    cat("DEBUG: miRNA data is NULL after read\n")
    return()
  }
  
  # Minimal validation – if any fail, show a notification and exit
  if (colnames(data)[1] != "X") {
    showNotification(
      "[miRNA counts] The first column header must be blank in the CSV.",
      type = "error", duration = NULL
    )
    return()
  }
  if (ncol(data) < 2) {
    showNotification(
      "[miRNA counts] File must have at least one sample column.",
      type = "error", duration = NULL
    )
    return()
  }
  
  mirna_ids <- data$X
  if (any(is.na(mirna_ids) | mirna_ids == "")) {
    showNotification(
      "[miRNA counts] The first column (miRNA IDs) contains empty or missing values.",
      type = "error", duration = NULL
    )
    return()
  }
  if (length(unique(mirna_ids)) != length(mirna_ids)) {
    showNotification(
      "[miRNA counts] The first column (miRNA IDs) contains duplicated values.",
      type = "error", duration = NULL
    )
    return()
  }
  
  sample_names <- colnames(data)[-1]
  if (any(sample_names == "" | is.na(sample_names))) {
    showNotification(
      "[miRNA counts] All sample columns must have non-empty column names (SampleIDs).",
      type = "error", duration = NULL
    )
    return()
  }
  
  # numeric coercion
  numeric_part <- suppressWarnings(as.data.frame(
    lapply(data[-1], function(x) as.numeric(as.character(x)))
  ))
  non_numeric_mask <- is.na(as.matrix(numeric_part)) & !is.na(as.matrix(data[-1]))
  if (any(non_numeric_mask)) {
    showNotification(
      "[miRNA counts] Some count values are non-numeric. Please provide raw integer counts.",
      type = "error", duration = NULL
    )
    return()
  }
  data[-1] <- numeric_part
  
  # Your original processing
  df_wide <- data
  df_wide[is.na(df_wide)] <- 0
  row.names(df_wide) <- df_wide$X
  df_wide <- df_wide %>% dplyr::select(-c(X))
  
  row_sub <- apply(df_wide, 1, function(row) all(row == 0))
  df_wide <- df_wide[!row_sub, , drop = FALSE]
  
  drop_colsum_zero <- (colSums(df_wide, na.rm = TRUE) != 0)
  df_wide <- df_wide[, drop_colsum_zero, drop = FALSE]
  
  wideData(df_wide)
  
  ld <- data %>%
    tidyr::pivot_longer(cols = -1, names_to = "SampleID", values_to = "Count")
  
  longData(ld)
  
  output$mirna_status <- renderText("miRNA count matrix successfully loaded.")
  cat("DEBUG: miRNA upload completed, status set\n")
})


# Module 0B: Read and Store Metadata
observeEvent(input$metadataFile, {
  req(input$metadataFile)
  cat("DEBUG: entered metadata upload observer\n")
  
  metadataData <- tryCatch(
    read.csv(input$metadataFile$datapath,
             header = TRUE,
             stringsAsFactors = FALSE,
             check.names = FALSE),
    error = function(e) {
      showNotification(
        paste("Error reading metadata file:", e$message),
        type = "error", duration = NULL
      )
      return(NULL)
    }
  )
  if (is.null(metadataData)) {
    cat("DEBUG: metadata is NULL after read\n")
    return()
  }
  
  # --- Minimal validation ---
  
  # 1) Required SampleID column
  if (!("SampleID" %in% names(metadataData))) {
    showNotification(
      "[Metadata] Required column 'SampleID' is missing. Please add a 'SampleID' column.",
      type = "error", duration = NULL
    )
    return()
  }
  
  sample_ids <- metadataData$SampleID
  
  # 2) SampleID values: non-empty
  if (any(is.na(sample_ids) | sample_ids == "")) {
    showNotification(
      "[Metadata] 'SampleID' column contains empty or missing values.",
      type = "error", duration = NULL
    )
    return()
  }
  
  # 3) SampleID values: no duplicates
  if (length(unique(sample_ids)) != length(sample_ids)) {
    showNotification(
      "[Metadata] 'SampleID' column contains duplicated sample IDs.",
      type = "error", duration = NULL
    )
    return()
  }
  
  # 4) SampleID character rules (same as for counts)
  invalid_samples <- sample_ids[!is_valid_identifier(sample_ids)]
  if (length(invalid_samples) > 0) {
    showNotification(
      paste0("[Metadata] Invalid 'SampleID' values: ",
             paste(unique(invalid_samples), collapse = ", "),
             ". SampleIDs may only contain letters, numbers, underscores (_) and periods (.)."),
      type = "error", duration = NULL
    )
    return()
  }
  
  # 5) If miRNA counts already loaded: all miRNA SampleIDs must be in metadata
  if (!is.null(wideData())) {
    mirna_samples   <- colnames(wideData())
    missing_in_meta <- setdiff(mirna_samples, sample_ids)
    
    if (length(missing_in_meta) > 0) {
      showNotification(
        paste0("[Metadata] The following SampleID(s) are present in the miRNA counts ",
               "but not in the metadata: ",
               paste(missing_in_meta, collapse = ", "),
               ". All SampleIDs in the miRNA counts file must be represented in the metadata."),
        type = "error", duration = NULL
      )
      return()
    }
  }
  
  # --- If we reach here, metadata is valid ---
  metadata(metadataData)
  
  updateSelectInput(session, "designColumn",      choices = names(metadataData))
  updateSelectInput(session, "groupBy",           choices = names(metadataData))
  updateSelectInput(session, "groupBy_Boxplot",   choices = names(metadataData))
  updateSelectInput(session, "metadataColumn",    choices = names(metadataData))
  
  output$metadata_status <- renderText("Metadata successfully loaded.")
  cat("DEBUG: metadata upload completed, status set\n")
})




# Module 0C: Read and Store mRNA Data
observeEvent(input$mrna_countTable, {
  req(input$mrna_countTable)
  cat("DEBUG: entered mRNA upload observer\n")
  
  data <- tryCatch(
    read.csv(input$mrna_countTable$datapath,
             header = TRUE),
    error = function(e) {
      showNotification(
        paste("Error reading mRNA count file:", e$message),
        type = "error", duration = NULL
      )
      return(NULL)
    }
  )
  if (is.null(data)) {
    cat("DEBUG: mRNA data is NULL after read\n")
    return()
  }
  
  if (colnames(data)[1] != "X") {
    showNotification(
      "[mRNA counts] The first column header must be blank in the CSV.",
      type = "error", duration = NULL
    )
    return()
  }
  if (ncol(data) < 2) {
    showNotification(
      "[mRNA counts] File must have at least one sample column.",
      type = "error", duration = NULL
    )
    return()
  }
  
  sample_names <- colnames(data)[-1]
  if (any(sample_names == "" | is.na(sample_names))) {
    showNotification(
      "[mRNA counts] All sample columns must have non-empty column names (SampleIDs).",
      type = "error", duration = NULL
    )
    return()
  }
  
  numeric_part <- suppressWarnings(as.data.frame(
    lapply(data[-1], function(x) as.numeric(as.character(x)))
  ))
  non_numeric_mask <- is.na(as.matrix(numeric_part)) & !is.na(as.matrix(data[-1]))
  if (any(non_numeric_mask)) {
    showNotification(
      "[mRNA counts] Some count values are non-numeric. Please provide raw integer counts.",
      type = "error", duration = NULL
    )
    return()
  }
  data[-1] <- numeric_part
  
  # SampleID match with miRNA (if miRNA is already loaded)
  if (!is.null(wideData())) {
    mirna_samples <- colnames(wideData())
    mrna_samples  <- sample_names
    if (!setequal(mirna_samples, mrna_samples)) {
      showNotification(
        "[mRNA counts] SampleIDs in mRNA counts must exactly match SampleIDs in miRNA counts.\n",
        type = "error", duration = NULL
      )
      return()
    }
  }
  
  mrna_wideData(data)
  output$mrna_status <- renderText("mRNA count matrix successfully loaded.")
  cat("DEBUG: mRNA upload completed, status set\n")
})


