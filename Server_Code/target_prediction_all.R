# Module 4: Predict target genes for all DE miRNA

### Helper functions ---

sanitize_entrez <- function(x) {
  x <- unique(as.character(x))
  x[!is.na(x) & nzchar(x)]
}

# Track the last valid min_src separately to avoid fighting with the slider
min_src_val <- reactiveVal(1L)

# Compute max based on current sources; at least 1 to keep slider valid
max_min_src <- reactive({
  n <- length(input$sources_agg %||% character())
  if (n < 1L) 1L else n
})

# When sources change, clamp the stored min_src to the new max
observeEvent(input$sources_agg, {
  nmax <- max_min_src()
  min_src_val(max(1L, min(min_src_val(), nmax)))
}, ignoreInit = TRUE)

# When the user moves the slider, store the new value
observeEvent(input$min_src, {
  # Only store if input$min_src exists and within current bounds
  nmax <- max_min_src()
  val <- input$min_src
  if (!is.null(val)) {
    min_src_val(max(1L, min(val, nmax)))
  }
}, ignoreInit = TRUE)

# Render the slider using the stored value and current max
output$min_src_ui <- renderUI({
  nmax <- max_min_src()
  value <- min_src_val()
  sliderInput(
    "min_src",
    "Minimum number of sources intersecting for target retrieval (2 recommended)",
    min = 1L, max = nmax, value = value, step = 1L
  )
})

########################

observeEvent(input$predict_all_genes, {
  
  #### Dev code to be commented out or deleted later ########

  # res_significant <- read.csv(here("dev_files/Significant_miRNA.csv"),row.names=1)
  # mirVector <- res_significant$miRNA
  # hsmart <- readRDS("dev_files/mart.RDS")
  # min_src <- 2
  # org_col <- "hgnc_symbol"
  # 
  # hsa_mirtarbase <- read.csv(here("data/hsa_MTI.csv"))
  # mirtarbase <- hsa_mirtarbase %>%
  #   filter(Species..Target.Gene. =="hsa") %>%
  #   dplyr::rename(entrezgene_id = Target.Gene..Entrez.ID.)
  # mirtarbase$entrezgene_id <- as.character(mirtarbase$entrezgene_id)
  # 
  # mmu_mirtarbase <- read.csv(here("data/mmu_MTI.csv")) %>%
  #   filter(Species..Target.Gene.=="mmu") %>%
  #   dplyr::rename(entrezgene_id = Target.Gene..Entrez.ID.)
  # mmu_mirtarbase$entrezgene_id <- as.character(mmu_mirtarbase$entrezgene_id)

  ###########################################################
  
  if (is.null(res_significant()) || nrow(res_significant()) == 0) {
    showNotification(
      HTML("Please run Module 3 (Differential miRNA Expression Analysis) first.<br>
            Click on the 'Differential miRNA Expression Analysis' tab and run DESeq2 analysis."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  # Require at least one source and compute bounded min_src
  sources_sel <- input$sources_agg %||% character(0)
  if (length(sources_sel) == 0) {
    showNotification("Please select at least one source.", type = "error", duration = 8)
    return()
  }
  # Use the same clamped min_src as in runAnalysis
  min_src <- min(min_src_val(), length(sources_sel))
  
  # User comfort message 
  ensembl_notification <- showNotification("Contacting Ensembl website...", type = "message", duration = NULL)
  
  # Require these values 
  req(res_significant(), mart(), org_string())
  
  # Load required values 
  mirVector <- res_significant()$miRNA
  if (length(mirVector) == 0) {
    showNotification("No DE miRNAs found.", type = "warning", duration = 8)
    return()
  }
  hsmart <- mart()
  org_col <- org_string()
  scope <- input$result_scope_all_agg
  

  if (is.null(hsmart)) {
    removeNotification(ensembl_notification)
    return()
  }
  
  removeNotification(ensembl_notification)
  
  # Map aggregated scope to base scope decisions
  # input$result_scope_all_agg in {"validated_agg","predicted_agg","both_agg","all_agg"}
  agg_scope <- input$result_scope_all_agg %||% "all_agg"
  scope <- switch(
    agg_scope,
    "validated_agg" = "validated",
    "predicted_agg" = "predicted",
    "both_agg"      = "both",
    "all_agg"       = "all",
    "all"
  )
  
  gene_notification <- showNotification("Retrieving targets for all DE miRNAs...", type = "message", duration = NULL)
  
  # Prepare collections
  empty_schema <- make_empty_schema(org_col)
  predicted_list  <- list()
  
  # For predicted warnings
  minSrcReducedPattern <- "sources which returned a target list<min_src, min_src reduced to"
  
  if (!identical(scope, "predicted")) {
    
    # Ensure org_mirtarbase() returned data
    if (is.null(org_mirtarbase())) {
      removeNotification(ensembl_notification)
      removeNotification(gene_notification)
      showNotification("No validated data, did you download and upload the miRTarBase file? Please do this or adjust the scope to 'Predicted only'.", type = "error", duration = 8)
      return()  # clean early exit
    }
    
    mirtarbase <- org_mirtarbase()
    
    # Pull all validated rows for DE miRNAs once
    vt_all <- tryCatch(
      mirtarbase %>%
        dplyr::filter(miRNA %in% mirVector) %>%
        dplyr::distinct(),
      error = function(e) empty_schema
    )
    
    if (!is_empty_df(vt_all) && nrow(vt_all) > 0) {
      # Collect unique, non-empty Entrez IDs
      vals <- unique(vt_all$entrezgene_id)
      vals <- vals[!is.na(vals) & vals != ""]
      
      # Single Ensembl mapping for all validated genes
      val_mapping <- if (length(vals) > 0) {
        getBM(
          attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_col),
          filters = 'entrezgene_id',
          values = vals,
          mart = hsmart
        )
      } else data.frame()
      if (is.null(val_mapping)) val_mapping <- data.frame()
      val_mapping$entrezgene_id <- as.factor(val_mapping$entrezgene_id)
      
      # Merge and normalize once
      vt_all <- dplyr::left_join(vt_all, val_mapping, by = "entrezgene_id") %>% dplyr::distinct()
      if (!("miRNA" %in% names(vt_all))) vt_all[["miRNA"]] <- NA_character_
      if (!(org_col %in% names(vt_all))) vt_all[[org_col]] <- NA_character_
      if (!("ensembl_gene_id" %in% names(vt_all))) vt_all[["ensembl_gene_id"]] <- NA_character_
      if (!("description" %in% names(vt_all))) vt_all[["description"]] <- NA_character_
      vt_all[["rank_final"]] <- NA_real_
      vt_all <- normalize_columns(vt_all, org_col)
      vt_all <- vt_all %>% dplyr::mutate(dplyr::across(!rank_final, as.character))
      
    } else {
      showNotification("No validated miRNA-gene interactions exist in miRTarBase for the input miRNA list.", type = "error", duration = 8)
    }
  }
  
  
  #########################################
  # PREDICTED block with single mapping   #
  #########################################
  predicted_list <- list()
  
  if (!identical(scope, "validated")) {
    withProgress(message = 'Searching databases (predicted)', value = 0, {
      n <- length(mirVector)
      for (i in seq_len(n)) {
        incProgress(1/n, detail = paste("Predicted targets for DE miRNA #", i))
        mir <- mirVector[i]
        
        # 1) Run predictions with prefix
        predictions <- NULL
        attempt1 <- captureWarnings(
          getPredictedTargets(
            mir,
            species = org_abbrev(),
            sources = sources_sel,
            method = "geom",
            min_src = min_src
          )
        )
        if (any(grepl(minSrcReducedPattern, attempt1$warnings))) {
          mir_no_prefix <- sub(paste0("^", org_abbrev(), "-"), "", mir)
          showNotification(
            paste("No predicted targets found for miRNA:", mir_no_prefix,
                  ". Suggest reducing the intersection threshold."),
            type = "error", duration = 2
          )
        } else {
          predictions <- attempt1$result
        }
        
        # 2) Retry without prefix if needed
        if (is_empty_df(predictions)) {
          mir_no_prefix <- sub(paste0("^", org_abbrev(), "-"), "", mir)
          attempt2 <- captureWarnings(
            getPredictedTargets(
              mir_no_prefix,
              species = org_abbrev(),
              sources = sources_sel,
              method = "geom",
              min_src = min_src
            )
          )
          if (any(grepl(minSrcReducedPattern, attempt2$warnings))) {
            showNotification(
              paste("No predicted targets found for miRNA:", mir_no_prefix,
                    ". Suggest reducing the intersection threshold."),
              type = "error", duration = 2
            )
          } else {
            predictions <- attempt2$result
          }
        }
        
        # 3) If we have predictions, keep only minimal columns now; mapping later
        if (!is_empty_df(predictions) && nrow(predictions) > 0) {
          df <- as.data.frame(predictions)
          df$entrezgene_id <- row.names(df)
          df$miRNA <- mir
          
          # Ensure rank_final exists; keep only fields we need now
          if (!("rank_final" %in% names(df))) df[["rank_final"]] <- NA_real_
          df <- df[, intersect(c("entrezgene_id", "miRNA", "rank_final", names(df)), names(df)), drop = FALSE]
          
          predicted_list[[mir]] <- df
        } else {
          mir_no_prefix <- sub(paste0("^", org_abbrev(), "-"), "", mir)
          showNotification(
            paste("No predicted targets found for miRNA:", mir_no_prefix,
                  ". Suggest reducing the intersection threshold."),
            type = "error", duration = 2
          )
        }
      } # end for
    }) # end withProgress
  }
  
  
  removeNotification(gene_notification)
  
  # After the loop: single mapping for all predicted Entrez IDs
  predicted_df_all <- {
    if (length(predicted_list) == 0) {
      empty_schema
    } else {
      # Row-bind all minimal predicted frames
      tmp <- tryCatch(dplyr::bind_rows(predicted_list), error = function(e) do.call(rbind, predicted_list))
      if (is.null(tmp) || nrow(tmp) == 0) {
        empty_schema
      } else {
        tmp
      }
    }
  }
  
  if (!identical(scope, "validated") && nrow(predicted_df_all) > 0) {
    # Unique IDs
    all_ids <- unique(predicted_df_all$entrezgene_id)
    all_ids <- all_ids[!is.na(all_ids) & all_ids != ""]
    
    # Single Ensembl lookup
    mapping <- if (length(all_ids) > 0) {
      getBM(
        attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_string()),
        filters = 'entrezgene_id',
        values = all_ids,
        mart = hsmart
      )
    } else data.frame()
    if (is.null(mapping)) mapping <- data.frame()
    
    # Merge mapping once
    predicted_df_all <- merge(predicted_df_all, mapping, by = "entrezgene_id", all.x = TRUE)
    
    # Ensure validated columns exist as NA and schema is consistent
    for (nm in c("miRTarBase.ID","Experiments","Support.Type","References..PMID.")) {
      if (!(nm %in% names(predicted_df_all))) predicted_df_all[[nm]] <- NA_character_
    }
    if (!("rank_final" %in% names(predicted_df_all))) predicted_df_all[["rank_final"]] <- NA_real_
    
    predicted_df_all <- normalize_columns(predicted_df_all, org_string())
  }
  
  predicted_df_all <- if (identical(scope, "validated")) empty_schema else predicted_df_all
  validated_df_all <- if (identical(scope, "predicted")) empty_schema else vt_all
  
  # Build combined based on scope logic
  key_cols <- c("entrezgene_id", "ensembl_gene_id", "miRNA", "description", org_col)
  
  validated_only <- validated_df_all %>%
    dplyr::filter(!is.na(miRTarBase.ID) & miRTarBase.ID != "") %>%
    dplyr::select(c(dplyr::all_of(key_cols), "miRTarBase.ID","Experiments","Support.Type","References..PMID."))
  
  predicted_only <- predicted_df_all %>%
    dplyr::filter(!is.na(rank_final)) %>%
    dplyr::select(c(dplyr::all_of(key_cols), "rank_final"))
  
  combined_df <- dplyr::full_join(
    predicted_only,
    validated_only,
    by = key_cols,
    relationship = "many-to-many"
  )
  
  combined_df <- normalize_columns(combined_df, org_col)
  
  # Flags
  validated_flag <- !is.na(combined_df$miRTarBase.ID) & combined_df$miRTarBase.ID != ""
  predicted_flag <- !is.na(combined_df$rank_final)
  
  # Apply aggregated scope filter (validated/predicted/both/all)
  finalgeneResults <- switch(
    scope,
    "validated" = combined_df[validated_flag, , drop = FALSE],
    "predicted" = combined_df[predicted_flag & !validated_flag, , drop = FALSE],
    "both"      = combined_df[predicted_flag & validated_flag, , drop = FALSE],
    "all"       = combined_df[predicted_flag | validated_flag, , drop = FALSE],
    combined_df[predicted_flag | validated_flag, , drop = FALSE]
  )
  finalgeneResults <- normalize_columns(finalgeneResults, org_col)
  
  # Annotate Source for transparency
  finalgeneResults$Source <- dplyr::case_when(
    (!identical(scope, "validated")) & (!is.na(finalgeneResults$rank_final)) &
      (!identical(scope, "predicted")) & (!is.na(finalgeneResults$miRTarBase.ID)) & finalgeneResults$miRTarBase.ID != "" ~ "Both",
    (!identical(scope, "validated")) & (!is.na(finalgeneResults$rank_final)) ~ "Predicted",
    (!identical(scope, "predicted")) & (!is.na(finalgeneResults$miRTarBase.ID)) & finalgeneResults$miRTarBase.ID != "" ~ "Validated",
    TRUE ~ "Unknown"
  )
  
  # if (!identical(scope, "validated_agg")) {
  #   withProgress(message = 'Searching databases', value = 0, {
  #     n <- length(mirVector)
  #     
  #     for (i in 1:n) {
  #       mir <- mirVector[i]
  #       predictions <- getPredictedTargets(mir, 
  #                                          species = org_abbrev(), 
  #                                          method = 'geom', 
  #                                          min_src = 2)
  #       
  #       if (is.null(predictions) || nrow(predictions) == 0) {
  #         mir_no_prefix <- sub(paste0("^",org_abbrev(),"-"), "", mir)
  #         message(paste("Trying again without prefix for miRNA:", mir_no_prefix))
  #         
  #         predictions <- getPredictedTargets(mir_no_prefix, 
  #                                            species = org_abbrev(), 
  #                                            method = 'geom', 
  #                                            min_src = 2)
  #         
  #         if (is.null(predictions) || nrow(predictions) == 0) {
  #           message(paste("No targets found for miRNA (after removing prefix):", mir_no_prefix))
  #           next
  #         }
  #       }
  #       
  #       df <- as.data.frame(predictions)
  #       
  #       df$entrezgene_id <- row.names(df)
  #       mapping <- getBM(
  #         attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_string()),
  #         filters = 'entrezgene_id',
  #         values = df$entrezgene_id,
  #         mart = hsmart
  #       )
  #       
  #       df <- merge(df, mapping, by="entrezgene_id")
  #       
  #       df$miRNA <- mir
  #       
  # 
  #       # Add miRTarBase 2025 results - 
  #       filtered_mirtarbase <- org_mirtarbase() %>% 
  #         filter(miRNA == mir)
  #       
  #       df <- df %>%
  #         full_join(filtered_mirtarbase, by= c("miRNA"="miRNA", "entrezgene_id" = "entrezgene_id"), relationship = "many-to-many")
  #       
  #       geneResults[[mir]] <- df[, c("entrezgene_id",
  #                                  "rank_final",
  #                                  "ensembl_gene_id",
  #                                  "description",
  #                                  org_string(),
  #                                  "miRNA",
  #                                  "miRTarBase.ID",
  #                                  "Experiments",
  #                                  "Support.Type",
  #                                  "References..PMID.")]
  #       
  #       incProgress(1/n, detail = paste("Getting predicted targets for DE miRNA #", i))
  #     }
  #   })
  # }
  # 
  # finalgeneResults <- do.call(rbind,geneResults)
  # 
  # 
  # # Flags
  # validated_flag <- !is.na(finalgeneResults$miRTarBase.ID) & finalgeneResults$miRTarBase.ID != ""
  # 
  # # rank_final is the predicted aggregate rank
  # predicted_flag <- !is.na(finalgeneResults$rank_final)
  # 
  # # annotate source for transparency
  # finalgeneResults$Source <- dplyr::case_when(
  #   predicted_flag & validated_flag ~ "Both",
  #   validated_flag ~ "Validated",
  #   predicted_flag ~ "Predicted",
  #   TRUE ~ "Unknown"
  # )
  # 
  # # Apply scope filter (ensure this matches the UI inputId)
  # scope <- input$result_scope_all_agg
  # if (identical(scope, "validated_agg")) {
  #   finalgeneResults <- finalgeneResults[validated_flag, , drop = FALSE]
  # } else if (identical(scope, "predicted_agg")) {
  #   finalgeneResults <- finalgeneResults[predicted_flag & !validated_flag, , drop = FALSE]
  # } else if (identical(scope, "both_agg")) {
  #   finalgeneResults <- finalgeneResults[predicted_flag & validated_flag, , drop = FALSE]
  # } else if (identical(scope, "all_agg")) {
  #   finalgeneResults <- finalgeneResults[predicted_flag | validated_flag, , drop = FALSE]
  # }
  
  predicted_target_reactive(finalgeneResults)
  print(head(finalgeneResults))

  output$all_genetargets <- renderDT({
    if (!is.null(finalgeneResults)) {
      dat <- predicted_target_reactive()
      datatable(dat, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))

    } else {
      DT::datatable(
        data.frame(Message = "No targets to display"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
      return()
    }
  })
  
  # Final results already stored in predicted_target_reactive
  
  ########################## DEV_ CODE ###########################
  
  # finalgeneResults <- readr::read_csv("dev_files/gene_results_2025-11-07_.csv")
  # names(finalgeneResults)
  # res_sig_df <- read.csv("dev_files/Significant_miRNA.csv", row.names=1)
  # df <- finalgeneResults %>%
  #        dplyr::select(miRNA, hgnc_symbol, entrezgene_id,
  #             description, miRTarBase.ID, Support.Type, Experiments, `References..PMID.`) %>%
  #        dplyr::filter(!is.na(miRNA), !is.na("hgnc_symbol"), !is.na(entrezgene_id)) %>%
  #        dplyr::distinct()
  # 
  # edges <- edges_df %>%
  #   dplyr::transmute(
  #     from = miRNA,
  #     to   = hgnc_symbol,
  #     title = paste0(
  #       "<b>miRNA:</b> ", miRNA,
  #       "<br><b>Gene:</b> ", hgnc_symbol,
  #       ifelse(!is.na(miRTarBase.ID), paste0("<br><b>miRTarBase:</b> ", miRTarBase.ID), ""),
  #       ifelse(!is.na(description), paste0("<br><b>description:</b> ", description), ""),
  #       ifelse(!is.na(Support.Type), paste0("<br><b>Support:</b> ", Support.Type), ""),
  #       ifelse(!is.na(Experiments), paste0("<br><b>Experiments:</b> ", Experiments), ""),
  #       ifelse(!is.na(`References..PMID.`), paste0("<br><b>PMID:</b> ", `References..PMID.`), "")
  #     )
  #   ) %>%
  #   dplyr::distinct(from, to, .keep_all = TRUE)  # ensure unique edges
  # 
  # upregulated_miRNA <- res_sig_df %>%
  #   filter(log2FoldChange > 0 ) %>%
  #   dplyr::select(miRNA, log2FoldChange)
  # upregulated_miRNA <- inner_join(upregulated_miRNA, finalgeneResults, by= "miRNA")
  # allPredictedGenes <- unique(downregulated_miRNA$entrezgene_id) %>% as.character()
  
  ######################################################3####
  
  # Determine which fold change column to use (DESeq2 vs limma)
  res_sig_df <- res_significant() %>% as.data.frame()
  fc_column <- if ("log2FoldChange" %in% colnames(res_sig_df)) "log2FoldChange" else "logFC"
  
  upregulated_miRNA <- res_sig_df %>% 
    filter(!!sym(fc_column) > 0 ) %>%
    dplyr::select(miRNA, !!fc_column)
  upregulated_miRNA <- inner_join(upregulated_miRNA, finalgeneResults, by= "miRNA")
  upregulated_miRNA_reactive(upregulated_miRNA)
  upregulated_miRNA_genes_reactive(unique(upregulated_miRNA$entrezgene_id))
  
  downregulated_miRNA <- res_sig_df %>% 
    filter(!!sym(fc_column) < 0 ) %>%
    dplyr::select(miRNA, !!fc_column)
  downregulated_miRNA <- inner_join(downregulated_miRNA, finalgeneResults, by= "miRNA")
  downregulated_miRNA_reactive(downregulated_miRNA)
  downregulated_miRNA_genes_reactive(unique(downregulated_miRNA$entrezgene_id))
})

# Helpers to get name vectors safely
up_miRNA_names <- reactive({
  out <- tryCatch(upregulated_miRNA_reactive()$miRNA, error = function(e) NULL)
  unique(out %||% character(0))
})

down_miRNA_names <- reactive({
  out <- tryCatch(downregulated_miRNA_reactive()$miRNA, error = function(e) NULL)
  unique(out %||% character(0))
})

all_miRNA_names <- reactive({
  req(predicted_target_reactive())
  sort(unique(predicted_target_reactive()$miRNA))
})

# Update the miRNA_filter choices when direction or sources change
observeEvent(
  list(input$miRNA_dir, up_miRNA_names(), down_miRNA_names(), all_miRNA_names()),
  {
    dir <- input$miRNA_dir %||% "both"
    
    choices <- switch(
      dir,
      "up"   = sort(up_miRNA_names()),
      "down" = sort(down_miRNA_names()),
      "both" = all_miRNA_names()
    )
    
    # Preserve any previously selected items that remain valid
    prev_sel <- input$miRNA_filter %||% character(0)
    new_sel  <- intersect(prev_sel, choices)
    
    updateSelectInput(
      session, "miRNA_filter",
      choices  = choices,
      selected = new_sel
    )
  },
  ignoreInit = FALSE
)

#Your edges reactive: filter by direction, then by selected miRNAs
miRNA_gene_edges_reactive <- reactive({
  req(predicted_target_reactive())
  finalgeneResults <- predicted_target_reactive()
  
  df <- finalgeneResults %>%
    dplyr::select(miRNA, !!sym(org_string()), entrezgene_id,
                  description, miRTarBase.ID, Support.Type, Experiments, `References..PMID.`) %>%
    dplyr::filter(!is.na(miRNA), !is.na(!!sym(org_string())), !is.na(entrezgene_id)) %>%
    dplyr::distinct()
  
  # Apply direction filter consistently with the dropdown
  dir <- input$miRNA_dir %||% "both"
  if (dir == "up") {
    df <- df %>% dplyr::filter(miRNA %in% up_miRNA_names())
  } else if (dir == "down") {
    df <- df %>% dplyr::filter(miRNA %in% down_miRNA_names())
  } else {
    # both: no restriction
  }
  
  # Then apply user selections within that direction
  if (!is.null(input$miRNA_filter) && length(input$miRNA_filter) > 0) {
    df <- df %>% dplyr::filter(miRNA %in% input$miRNA_filter)
  }
  
  # validate to show a clear message if empty
  validate(need(nrow(df) > 0, "No data after direction and selection filters."))
  
  df
})

observeEvent(input$showGeneNetwork, {
  req(miRNA_gene_edges_reactive(), upregulated_miRNA_reactive(), downregulated_miRNA_reactive())
  
  edges_df <- miRNA_gene_edges_reactive()
  
  # Apply direction-based gene filtering using Entrez IDs
  # These reactives should return vectors of entrezgene_id
  up_mirnas   <- tryCatch(upregulated_miRNA_reactive() %>% pull(miRNA) %>% unique(), error = function(e) NULL)
  down_mirnas <- tryCatch(downregulated_miRNA_reactive() %>% pull(miRNA) %>% unique(), error = function(e) NULL)
  
  # Apply miRNA-direction filtering
  if (!is.null(input$miRNA_dir)) {
    if (input$miRNA_dir == "up") {
      edges_df <- edges_df %>% dplyr::filter(miRNA %in% up_mirnas)
    } else if (input$miRNA_dir == "down") {
      edges_df <- edges_df %>% dplyr::filter(miRNA %in% down_mirnas)
    } else {
      # "both": no restriction by direction
      edges_df <- edges_df
    }
  }
  
  # Early exit if nothing remains
  validate(need(nrow(edges_df) > 0, "No edges after filtering by miRNA direction."))
  
  # Build edges (unique pairs are recommended to avoid inflated degree)
  edges <- edges_df %>%
    dplyr::filter(!is.na(!!sym(org_string())) & !!sym(org_string()) != "") %>%
    dplyr::transmute(
      from = miRNA,
      to   = !!sym(org_string()),
      title = paste0(
        "<b>miRNA:</b> ", miRNA,
        "<br><b>Gene:</b> ", !!sym(org_string()),
        ifelse(!is.na(miRTarBase.ID), paste0("<br><b>miRTarBase:</b> ", miRTarBase.ID), ""),
        ifelse(!is.na(description), paste0("<br><b>description:</b> ", description), ""),
        ifelse(!is.na(Support.Type), paste0("<br><b>Support:</b> ", Support.Type), ""),
        ifelse(!is.na(Experiments), paste0("<br><b>Experiments:</b> ", Experiments), ""),
        ifelse(!is.na(`References..PMID.`), paste0("<br><b>PMID:</b> ", `References..PMID.`), "")
      )
    ) %>%
    dplyr::distinct(from, to, .keep_all = TRUE)  # ensure unique edges

  
  # Nodes from current edges
  miRNA_nodes <- edges %>% dplyr::distinct(from) %>%
    dplyr::transmute(id = from, label = from, group = "miRNA", color = "#6EC5FF")
  gene_nodes <- edges %>% dplyr::distinct(to) %>%
    dplyr::transmute(id = to, label = to, group = "gene", color = "#B184F0")
  nodes <- dplyr::bind_rows(miRNA_nodes, gene_nodes)
  
  # Optional filtering by min_degree and max_nodes will change the graph; we’ll recompute degree after these filters
  # First, compute a preliminary degree to allow pruning decisions
  prelim_deg <- dplyr::bind_rows(
    edges %>% dplyr::count(from, name = "deg") %>% dplyr::rename(id = from),
    edges %>% dplyr::count(to,   name = "deg") %>% dplyr::rename(id = to)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(degree = sum(deg), .groups = "drop")
  nodes <- nodes %>% dplyr::left_join(prelim_deg, by = "id") %>%
    dplyr::mutate(degree = ifelse(is.na(degree), 0L, degree))
  
  # Filter by min_degree (uses preliminary degree just to decide pruning)
  if (!is.null(input$min_degree) && input$min_degree > 0) {
    keep_ids <- nodes %>%
      dplyr::filter(degree >= input$min_degree) %>%
      dplyr::pull(id)
    edges <- edges %>% dplyr::filter(from %in% keep_ids, to %in% keep_ids)
    nodes <- nodes %>% dplyr::filter(id %in% keep_ids)
  }
  
  # If graph still too large, keep top N nodes by preliminary degree
  if (nrow(nodes) > input$max_nodes) {
    keep <- nodes %>%
      dplyr::arrange(dplyr::desc(degree)) %>%
      dplyr::slice_head(n = input$max_nodes) %>%
      dplyr::pull(id)
    nodes <- nodes %>% dplyr::filter(id %in% keep)
    edges <- edges %>% dplyr::filter(from %in% keep, to %in% keep)
  }
  
  # IMPORTANT: recompute degree based on the final edges actually shown
  final_deg <- dplyr::bind_rows(
    edges %>% dplyr::count(from, name = "deg") %>% dplyr::rename(id = from),
    edges %>% dplyr::count(to,   name = "deg") %>% dplyr::rename(id = to)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(degree = sum(deg), .groups = "drop")
  
  nodes <- nodes %>%
    dplyr::select(-degree) %>%               # remove preliminary degree
    dplyr::left_join(final_deg, by = "id") %>%
    dplyr::mutate(degree = ifelse(is.na(degree), 0L, degree)) %>%
    dplyr::mutate(
      title = ifelse(group == "miRNA",
                     paste0("<b>miRNA:</b> ", label, "<br><b>Degree (shown after filtering):</b> ", degree),
                     paste0("<b>Gene:</b> ", label, "<br><b>Degree (shown after filtering):</b> ", degree))
    )
  
  # store exactly what will be plotted
  network_displayed(list(nodes = nodes, edges = edges))
  
  output$miRNAGeneNetwork <- visNetwork::renderVisNetwork({
    visNetwork::visNetwork(nodes, edges, height = "700px", width = "100%") %>%
      visNetwork::visNodes(shape = "dot", size = 18, font = list(size = 18)) %>%
      visNetwork::visEdges(smooth = FALSE, color = list(color = "#B3B3B3"), arrows = "none") %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                             nodesIdSelection = TRUE) %>%
      visNetwork::visLegend(addNodes = list(
        list(label = "miRNA", shape = "dot", color = "#6EC5FF"),
        list(label = "Gene",  shape = "dot", color = "#B184F0")
      ),
      useGroups = FALSE) %>%
      visNetwork::visInteraction(hover = TRUE, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
      visNetwork::visIgraphLayout(randomSeed = 123)
  })
  
  ######################## NEW 2026 ################################
  output$dl_network_nodes <- downloadHandler(
    filename = function() paste0("displayed_network_nodes_", Sys.Date(), ".csv"),
    content = function(file) {
      req(network_displayed()$nodes)
      readr::write_csv(network_displayed()$nodes, file)
    }
  )
  
  output$dl_network_edges <- downloadHandler(
    filename = function() paste0("displayed_network_edges_", Sys.Date(), ".csv"),
    content = function(file) {
      req(network_displayed()$edges)
      readr::write_csv(network_displayed()$edges, file)
    }
  )
  ###################################################################
})

observeEvent(input$runORA, {
  req((input$hsa || input$mmu), 
      upregulated_miRNA_reactive(), 
      downregulated_miRNA_reactive())
  
  pathway_notification <- showNotification("Looking for enriched pathways. This may take some time", type = "message", duration = NULL)
  
  ### TARGET SET ---
  if (input$select_up_or_down_miRNA == "Up") {
    
    #allPredictedGenes <- unique(upregulated_miRNA$entrezgene_id)
    allPredictedGenes <- upregulated_miRNA_genes_reactive()
  } else if (input$select_up_or_down_miRNA == "Down") {
    allPredictedGenes <- downregulated_miRNA_genes_reactive()
  } else {
    allPredictedGenes <- character(0)
  }
  allPredictedGenes <- sanitize_entrez(allPredictedGenes) # CHANGE: sanitize target IDs
  
  if (length(allPredictedGenes) == 0) {
    removeNotification(pathway_notification)
    showNotification("No target genes available for ORA.", type = "error", duration = 8)
    return()
  }
  
  #print(head(allPredictedGenes))
  
  
  ### BACKGROUND SELECTION ---
  scope_agg <- input$result_scope_all_agg
  
  if (identical(input$select_background, "Organism_all")) { # CHANGE: new explicit option
    background_genes <- keys(org_db(), keytype = "ENTREZID")
    
  } else if (identical(input$select_background, "Custom_list")) { # CHANGE: new custom input option
    # Expect a textInput where users paste comma-separated ENTREZ IDs, e.g., input$custom_universe_ids
    custom_str <- input$custom_universe_ids %||% ""
    # Parse comma-separated list safely
    background_genes <- if (nzchar(custom_str)) {
      trimws(unlist(strsplit(custom_str, ",")))
    } else character(0)
    
  } else if (identical(input$select_background, "Validated_universe") && identical(scope_agg, "validated_agg")) { # CHANGE: conditional validated universe
    # Build universe from all miRTarBase targets (across organism), using ENTREZ IDs
    mtb <- org_mirtarbase()
    # Ensure column name matches exactly "entrezgene_id" in your miRTarBase table
    background_genes <- mtb$entrezgene_id
  } else {
    background_genes <- keys(org_db(), keytype = "ENTREZID")
  }
  
  background_genes <- sanitize_entrez(background_genes) # CHANGE: sanitize universe
  
  # Enforce target subset of background to avoid errors
  if (!all(allPredictedGenes %in% background_genes)) {
    allPredictedGenes <- intersect(allPredictedGenes, background_genes) # CHANGE: enforce subset
  }
  
  if (length(background_genes) == 0 || length(allPredictedGenes) == 0) {
    removeNotification(pathway_notification)
    showNotification("Empty background or target set for ORA.", type = "error", duration = 8)
    return()
  }
  
  ### RUN PATHWAY ORA ---
  
  if (input$analysis_type == "GO") {
    ora_results <- enrichGO(gene = allPredictedGenes,
                            universe = background_genes,
                            OrgDb = org_db(),
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = input$ORA_all_padj_Cutoff,
                            qvalueCutoff = input$ORA_all_padj_Cutoff,
                            readable = TRUE)
  } else if (input$analysis_type == "Reactome") {
    ora_results <- enrichPathway(gene = allPredictedGenes,
                                 organism = organism_reactive(),
                                 pAdjustMethod = "BH",
                                 universe=background_genes,
                                 pvalueCutoff = input$ORA_all_padj_Cutoff,
                                 qvalueCutoff = input$ORA_all_padj_Cutoff,
                                 readable = TRUE)
  }
  
  # Store table form safely
  All_miRNA_Pathways_reactive(if (!is.null(ora_results)) as.data.frame(ora_results) else data.frame())
  
  # Handle empty results
  if (is.null(as.data.frame(ora_results)) || nrow(as.data.frame(ora_results)) == 0) {
    output$all_resultsTable <- renderDT({
      DT::datatable(
        data.frame(Message = "No enriched pathways found"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
    })
    removeNotification(pathway_notification)
    return()
  }
  
  target_genes <- ora_results@result$geneID
  pathways <- ora_results@result$ID
  pathway_description <- ora_results@result$Description
  fold_enrichment <- ora_results@result$FoldEnrichment
  q_value <- ora_results@result$qvalue
  
  target_genes_list <- strsplit(target_genes, "/")
  
  pathway_genes_df <- data.frame(
    pathway = rep(pathways, times = sapply(target_genes_list, length)),
    description = rep(pathway_description, times = sapply(target_genes_list, length)),
    enrichment = rep(fold_enrichment, times = sapply(target_genes_list, length)),
    q_value = rep(q_value, times = sapply(target_genes_list, length)),
    gene = unlist(target_genes_list)
  )
  gene_count <- pathway_genes_df %>%
    group_by(pathway) %>%
    summarise(total_genes = n_distinct(gene))
  
  if (input$select_up_or_down_miRNA == "Up") {
    mirna_gene_results <- upregulated_miRNA_reactive()
  } else if (input$select_up_or_down_miRNA == "Down") {
    mirna_gene_results <- downregulated_miRNA_reactive()
  }
  
  print(head(allPredictedGenes))
  
  miRNA_df <- as.data.frame(mirna_gene_results) %>% dplyr::select(org_string(),"miRNA")
  mapped_df <- miRNA_df %>%
    inner_join(pathway_genes_df, by = setNames("gene", org_string()), relationship = "many-to-many")
  mapped_df <- mapped_df %>%
    group_by(miRNA, pathway) %>%
    mutate(gene_count = n()) %>%
    ungroup()
  mapped_df <- mapped_df %>% 
    inner_join(gene_count, by = "pathway")
  mapped_df <- mapped_df %>%
    group_by(miRNA, pathway) %>%
    mutate(coverage = gene_count/total_genes) %>%
    ungroup()
  
  miRNA_mapped_pathways(mapped_df)
  
  output$all_resultsTable <- renderDT({
    if (!is.null(as.data.frame(ora_results))) {
      datatable(as.data.frame(ora_results), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    } else {
      DT::datatable(
        data.frame(Message = "No enriched pathways found"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
      return()
    }
  })
  
  removeNotification(pathway_notification)
  
})

# Download handlers (outside observeEvents)
output$all_downloadgeneResults <- downloadHandler(
  filename = function() {
    paste("gene_results", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(predicted_target_reactive())
    readr::write_csv(predicted_target_reactive(), file)
  }
)

output$all_downloadFinalResults <- downloadHandler(
  filename = function() {
    paste("ora_results_all_predicted_genes_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    req(All_miRNA_Pathways_reactive())
    data <- All_miRNA_Pathways_reactive()
    
    # Handle NULL or empty results
    if (is.null(data) || (is.data.frame(data) && nrow(data) == 0)) {
      # Create empty data frame with appropriate columns
      data <- data.frame(
        ID = character(),
        Description = character(),
        GeneRatio = character(),
        BgRatio = character(),
        pvalue = numeric(),
        p.adjust = numeric(),
        qvalue = numeric(),
        geneID = character(),
        Count = integer()
      )
    } else if (!is.data.frame(data)) {
      data <- as.data.frame(data)
    }
    
    write.csv(data, file, row.names = FALSE, quote = FALSE)
  }
)