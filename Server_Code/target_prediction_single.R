# Module 5: Predict target genes for a single DE miRNA
# Track the last valid min_src separately to avoid fighting with the slider
min_src_val <- reactiveVal(1L)

# Compute max based on current sources; at least 1 to keep slider valid
max_min_src <- reactive({
  n <- length(input$sources %||% character())
  if (n < 1L) 1L else n
})

# When sources change, clamp the stored min_src to the new max
observeEvent(input$sources, {
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
output$single_min_src_ui <- renderUI({
  nmax <- max_min_src()
  value <- min_src_val()
  sliderInput(
    "min_src",
    "Minimum number of sources intersecting for target retrieval (2 recommended)",
    min = 1L, max = nmax, value = value, step = 1L
  )
})

# output$debug <- renderPrint({
#   list(
#     sources = input$sources,
#     nmax = max_min_src(),
#     input_min_src = input$min_src,
#     stored_min_src = min_src_val()
#   )
# })



observeEvent(input$runAnalysis, {
  
  #### Dev code to be commented out or deleted later ########
  # 
  # res_significant <- read.csv(here("dev_files/Significant_miRNA.csv"),row.names=1)
  # mir <- res_significant$miRNA[1]
  # hsmart <- readRDS("dev_files/mart.RDS")
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
  # 
  ###########################################################
  
  if (is.null(res_significant()) || nrow(res_significant()) == 0) {
    showNotification(
      HTML("Please run Module 3 (Differential miRNA Expression Analysis) first.<br>
            Click on the 'Differential miRNA Expression Analysis' tab and run DESeq2 analysis."),
      type = "warning", duration = 10
    )
    return()
  }
  
  sources_sel <- input$sources
  if (length(sources_sel) == 0) {
    showNotification("Please select at least one source.", type = "error")
    return()
  }
  
  min_src <- min(input$min_src, length(sources_sel))
  gene_notification <- showNotification("Retrieving targets...", type = "message", duration = NULL)
  
  req(res_significant(), mart())
  
  mir <- input$DEMiRNA
  hsmart <- mart()
  org_col <- org_string()
  scope <- input$result_scope_all  # "validated", "predicted", "both", "all"
  
  # Initialize empty frames with final schema to ease binding later
  org_col <- org_string()
  empty_schema <- make_empty_schema(org_col)
  validated_df <- empty_schema
  predicted_df <- empty_schema
  
  # VALIDATED chunk: run if scope != "predicted"
  if (!identical(scope, "predicted")) {
    
    # Ensure org_mirtarbase() returned data
    if (is.null(org_mirtarbase())) {
      removeNotification(gene_notification)
      showNotification("No validated data, did you download and upload the miRTarBase file? Please do this or adjust the scope to 'Predicted only'.", type = "error", duration = 8)
      
      return(invisible(NULL))  # clean early exit
    }
    
    mirtarbase_raw <- org_mirtarbase() %>% dplyr::filter(miRNA == mir) %>% unique()
    validated_df <- mirtarbase_raw 
    
    val_mapping <- getBM(
          attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_string()),
          filters = 'entrezgene_id',
          values = validated_df$entrezgene_id,
          mart = hsmart
        ) 
    
    if (is.null(val_mapping)) val_mapping <- data.frame()
    val_mapping$entrezgene_id <- as.factor(val_mapping$entrezgene_id)
    validated_df <- left_join(validated_df, val_mapping, by = "entrezgene_id") %>% unique()
      
    # Ensure required columns and placeholders
    if (!("entrezgene_id" %in% names(validated_df))) validated_df[["entrezgene_id"]] <- NA_character_
    if (!("miRNA" %in% names(validated_df))) validated_df[["miRNA"]] <- mir
    if (!(org_string() %in% names(validated_df))) validated_df[[org_string()]] <- NA_character_
    if (!("ensembl_gene_id" %in% names(validated_df))) validated_df[["ensembl_gene_id"]] <- NA_character_
    if (!("description" %in% names(validated_df))) validated_df[["description"]] <- NA_character_
    validated_df[["rank_final"]] <- NA_real_  # no predicted rank in validated-only
    validated_df <- normalize_columns(validated_df, org_string())
    validated_df <- validated_df %>% mutate(across(!rank_final, as.character))
  }
  
  # PREDICTED chunk: run if scope != "validated"
  if (!identical(scope, "validated")) {
    attempt1 <- captureWarnings(
      getPredictedTargets(
        mir, species = org_abbrev(), sources = sources_sel,
        method = "geom", min_src = min_src
      )
    )
    minSrcReducedPattern <- "sources which returned a target list<min_src, min_src reduced to"

      if (any(grepl(minSrcReducedPattern, attempt1$warnings))) {
            # End computation and show your error notification
            mir_no_prefix <- sub(paste0("^", org_abbrev(), "-"), "", mir)
            removeNotification(gene_notification)
            showNotification(
              paste(
                "No predicted targets found for miRNA:",
                mir_no_prefix,
                ". Suggest reducing the intersection threshold."
              ),
              type = "error",
              duration = 10
            )
      } else {
            predictions <- attempt1$result
          }

    if (is_empty_df(predictions)) {
      mir_no_prefix <- sub(paste0("^", org_abbrev(), "-"), "", mir)
      attempt2 <- captureWarnings(
        getPredictedTargets(
          mir_no_prefix, species = org_abbrev(), sources = sources_sel,
          method = "geom", min_src = min_src
        )
      )
      if (any(grepl(minSrcReducedPattern, attempt2$warnings))) {
              removeNotification(gene_notification)
              showNotification(
                paste(
                  "No predicted targets found for miRNA:",
                  mir_no_prefix,
                  ". Suggest reducing the intersection threshold."
                ),
                type = "error",
                duration = 10
              )
      } else {
        predictions <- attempt2$result
      }
    }
    
    if (!is_empty_df(predictions)) {
      predicted_df <- as.data.frame(predictions)
      if (nrow(predicted_df) > 0) {
        predicted_df$entrezgene_id <- row.names(predicted_df)
        mapping <- getBM(
          attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_string()),
          filters = 'entrezgene_id',
          values = predicted_df$entrezgene_id,
          mart = hsmart
        )
        if (is.null(mapping)) mapping <- data.frame()
        predicted_df <- merge(predicted_df, mapping, by = "entrezgene_id", all.x = TRUE)
        predicted_df$miRNA <- mir
        # Ensure validated-specific columns exist as NA
        for (nm in c("miRTarBase.ID","Experiments","Support.Type","References..PMID.")) {
          if (!(nm %in% names(predicted_df))) predicted_df[[nm]] <- NA_character_
        }
        if (!("rank_final" %in% colnames(predicted_df))) predicted_df[["rank_final"]] <- NA_real_
        predicted_df <- normalize_columns(predicted_df, org_string())
      } else {
        showNotification(
          paste(
            "No predicted targets found for miRNA:",
            mir_no_prefix,
            ". Suggest reducing the intersection threshold."
          ),
          type = "error",
          duration = 10
        )
        predicted_df <- empty_schema
      }
    } else {
      showNotification(
        paste(
          "No predicted targets found for miRNA:",
          mir_no_prefix,
          ". Suggest reducing the intersection threshold."
        ),
        type = "error",
        duration = 10
      )
      predicted_df <- empty_schema
    }
  }
  
  # Combine what was computed
  key_cols <- c("entrezgene_id", "ensembl_gene_id", "miRNA","description",org_string())
  
  validated_only <- validated_df %>%
    filter(!is.na(miRTarBase.ID) & miRTarBase.ID != "") %>%
    dplyr::select(c(all_of(key_cols),  "miRTarBase.ID","Experiments","Support.Type","References..PMID."))
  
  predicted_only <- predicted_df %>%
    filter(!is.na(rank_final)) %>%
    dplyr::select(c(all_of(key_cols), "rank_final"))

  
  combined_df <- full_join(
      predicted_only,
      validated_only,
      by = key_cols,
      relationship = "many-to-many"
    )
  
  combined_df <- normalize_columns(combined_df, org_string())

  # Flags for filtering
  validated_flag <- !is.na(combined_df$miRTarBase.ID) & combined_df$miRTarBase.ID != ""
  predicted_flag <- !is.na(combined_df$rank_final)
  
  # Apply scope
  finalgeneResults <- switch(
    scope,
    "validated" = combined_df[validated_flag, , drop = FALSE],
    "predicted" = combined_df[predicted_flag, , drop = FALSE],
    "both"      = combined_df[predicted_flag & validated_flag, , drop = FALSE],
    "all"       = combined_df[predicted_flag | validated_flag, , drop = FALSE],
    combined_df[predicted_flag | validated_flag, , drop = FALSE]
  )
  
  # Ensure exact columns and order
  finalgeneResults <- normalize_columns(finalgeneResults, org_string())
  
  displayed_gene_results(finalgeneResults)
  
  output$genetargets <- renderDT({
    dat <- displayed_gene_results()
    if (is.null(dat) || nrow(dat) == 0) {
      DT::datatable(
        data.frame(Message = "No targets to display"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
    } else {
      DT::datatable(dat, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    }
  })
  
  removeNotification(gene_notification)
})


observeEvent(input$runORA_single_miRNA, {
  req(displayed_gene_results())
  
  pathway_notification <- showNotification("Looking for enriched pathways...", type = "message", duration = NULL)
  
  Single_miRNA_Predicted_Genes <- unique(displayed_gene_results()$entrezgene_id)
  
  ### BACKGROUND SELECTION ---
  scope <- input$result_scope_all
  
  if (identical(input$select_background_single, "Organism_all")) { # CHANGE: new explicit option
    background_genes <- keys(org_db(), keytype = "ENTREZID")
    
  } else if (identical(input$select_background_single, "Custom_list")) { # CHANGE: new custom input option
    # Expect a textInput where users paste comma-separated ENTREZ IDs, e.g., input$custom_universe_ids
    custom_str <- input$single_custom_universe_ids %||% ""
    # Parse comma-separated list safely
    background_genes <- if (nzchar(custom_str)) {
      trimws(unlist(strsplit(custom_str, ",")))
    } else character(0)
    
  } else if (identical(input$select_background_single, "Validated_universe") && identical(scope, "validated")) { # CHANGE: conditional validated universe
    # Build universe from all miRTarBase targets (across organism), using ENTREZ IDs
    mtb <- org_mirtarbase()
    # Ensure column name matches exactly "entrezgene_id" in your miRTarBase table
    background_genes <- mtb$entrezgene_id
  } else {
    background_genes <- keys(org_db(), keytype = "ENTREZID")
  }
  
  background_genes <- sanitize_entrez(background_genes) 
  
  # Enforce target subset of background to avoid errors
  if (!all(Single_miRNA_Predicted_Genes %in% background_genes)) {
    Single_miRNA_Predicted_Genes <- intersect(Single_miRNA_Predicted_Genes, background_genes) # CHANGE: enforce subset
  }
  
  if (length(background_genes) == 0 || length(Single_miRNA_Predicted_Genes) == 0) {
    removeNotification(pathway_notification)
    showNotification("Empty background or target set for ORA.", type = "error", duration = 8)
    return()
  }
  
  if (input$single_miRNA_analysis_type == "GO") {
    ora_results <- enrichGO(gene = Single_miRNA_Predicted_Genes,
                            universe = background_genes,
                            OrgDb = org_db(),
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = input$ORA_single_padj_Cutoff,
                            qvalueCutoff = input$ORA_single_padj_Cutoff,
                            readable = TRUE)
  } else if (input$single_miRNA_analysis_type == "Reactome") {
    ora_results <- enrichPathway(gene = Single_miRNA_Predicted_Genes,
                                 organism = organism_reactive(),
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = input$ORA_single_padj_Cutoff,
                                 qvalueCutoff = input$ORA_single_padj_Cutoff,
                                 readable = TRUE)
  }
  
  Single_miRNA_Pathways_reactive(ora_results)
  
  output$resultsTable <- renderDT({
    if (is.null(ora_results) || nrow(ora_results) == 0) {
      DT::datatable(
        data.frame(Message = "No targets to display"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
    } else {
      ora_results <-as.data.frame(ora_results) 
      DT::datatable(ora_results, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    }
  })
  
  removeNotification(pathway_notification)
})

observeEvent(input$plot_pathways_single_miRNA, {
  
  ora_results <- Single_miRNA_Pathways_reactive()
  
  if (is.null(ora_results) || nrow(ora_results) == 0) {
    output$ora_dotplot <- renderPlot({
      plot.new(); text(0.5, 0.5, "no significant features to plot")
    })
    output$ora_barplot <- renderPlot({
      plot.new(); text(0.5, 0.5, "no significant features to plot")
    })
    showNotification("No pathways to plot.", type = "error", duration = 8)
    dotplot_reactive(NULL)
    barplot_reactive(NULL)
    return()
  }
  
  output$ora_dotplot <- renderPlot({
    if (!is.null(ora_results)) {
      plot <- dotplot(ora_results, showCategory = 10) +
        ggtitle("Dotplot of Overrepresented Terms")
      dotplot_reactive(plot)
      plot
    } else 
      showNotification("No pathways to plot.", type = "error", duration = 8)
  })
  
  output$ora_barplot <- renderPlot({
    if (!is.null(ora_results)) {
      plot <- barplot(ora_results, showCategory = 10) +
        ggtitle("Barplot of Overrepresented Terms")
      barplot_reactive(plot)
      plot
    } else 
      showNotification("No pathways to plot.", type = "error", duration = 8)
    
  })
})

# # Download handlers (outside observeEvents)  
# output$downloadgeneResults <- downloadHandler(
#   filename = function() {
#     req(input$DEMiRNA)
#     paste("miRNA", input$DEMiRNA, "predicted_gene_results", Sys.Date(), ".csv", sep = "_")
#   },
#   content = function(file) {
#     req(input$DEMiRNA)
#     mir <- input$DEMiRNA
#     
#     hsmart <- mart()
#     predictions <- getPredictedTargets(mir, species = org_abbrev(), method = 'geom', min_src = 2)
#     
#     if (is.null(predictions) || nrow(predictions) == 0) {
#       mir_no_prefix <- sub(paste0("^",org_abbrev(),"-"), "", mir)
#       predictions <- getPredictedTargets(mir_no_prefix, species = org_abbrev(), method = 'geom', min_src = 2)
#     }
#     
#     if (!is.null(predictions) && nrow(predictions) > 0) {
#       df <- as.data.frame(predictions)
#       df$entrezgene_id <- row.names(df)
#       mapping <- getBM(
#         attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_string()), 
#         filters = 'entrezgene_id',
#         values = df$entrezgene_id,
#         mart = hsmart
#       )
#       df <- merge(df, mapping, by = "entrezgene_id")
#       df$miRNA <- mir
#       finalgeneResults <- df[, c("entrezgene_id", "rank_final", "ensembl_gene_id", "description", org_string(), "miRNA")]
#       write.csv(finalgeneResults, file, row.names = FALSE, quote = FALSE)
#     }
#   }
# )


# Download handlers (outside observeEvents)
output$downloadgeneResults <- downloadHandler(
  filename = function() {
    req(input$DEMiRNA, displayed_gene_results())
    paste("miRNA", input$DEMiRNA, "predicted_gene_results", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    # Get whatever is currently displayed (filtered or not)
    dat <- displayed_gene_results()
    req(dat)
    
    # If you need to ensure column order, specify it here.
    # Otherwise, write as-is.
    # Example: keep the same order as your finalgeneResults
    cols <- c(
      "entrezgene_id",
      "rank_final",
      "ensembl_gene_id",
      "description",
      org_string(),
      "miRNA",
      "miRTarBase.ID",
      "Experiments",
      "Support.Type",
      "References..PMID."
    )
    common_cols <- intersect(cols, colnames(dat))
    dat_out <- dat[, common_cols, drop = FALSE]
    
    write.csv(dat_out, file, row.names = FALSE, quote = FALSE)
  }
)

output$downloadFinalResults <- downloadHandler(
  filename = function() {
    req(input$DEMiRNA)
    paste("miRNA", input$DEMiRNA, "pathway_results", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(Single_miRNA_Pathways_reactive())
    ora_results <- Single_miRNA_Pathways_reactive()
    
    # Handle different types of results
    if (is.null(ora_results)) {
      # Create empty data frame with appropriate columns
      ora_df <- data.frame(
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
    } else if (inherits(ora_results, "enrichResult")) {
      # Extract result slot from S4 object
      ora_df <- ora_results@result
    } else {
      # Already a data frame
      ora_df <- as.data.frame(ora_results)
    }
    
    write.csv(ora_df, file, row.names = FALSE, quote = FALSE)
  }
)