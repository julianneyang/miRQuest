## Lines 250-286 have been adapted from https://rdrr.io/bioc/anamiR/src/R/differExp_discrete.R (AnamiR, GPL-2 license)
## Therefore, if you use this software, please also cite: Wang TT, Lee CY, Lai LC, Tsai MH, Lu TP, Chuang EY. anamiR: integrated analysis of MicroRNA and gene expression profiling. BMC Bioinformatics. 2019 May 14;20(1):239. doi: 10.1186/s12859-019-2870-x. PMID: 31088348; PMCID: PMC6518761.
# Module 3: Do Differential miRNA Expression Analysis with DESEQ2


# --- Populate UI choices based on metadata ---

# ADDED: update design covariate choices when metadata is ready
observeEvent(metadata(), {
  req(metadata())
  cols <- setdiff(names(metadata()), "SampleID")
  updateSelectizeInput(session, "designColumns", choices = cols, server = TRUE)
})

# ADDED: keep VOI in sync with selected design covariates
observeEvent(input$designColumns, {
  updateSelectInput(session, "voi", choices = input$designColumns)
})

# ADDED: show reference level selector if VOI is a factor
output$voiRefLevelUI <- renderUI({
  req(metadata(), input$voi)
  md <- metadata()
  if (input$voi %in% names(md)) {
    vals <- md[[input$voi]]
    lv <- levels(factor(vals))
    if (length(lv) > 0) {
      selectInput("voiRef", "Reference level for VOI:", choices = lv, selected = lv[1])
    } else {
      # Non-factor (numeric/character without clear levels): no reference level UI
      NULL
    }
  } else {
    NULL
  }
})


# Helper: check if the design matrix is full rank
is_full_rank <- function(colData_df, design_formula) {
  mm <- tryCatch(
    stats::model.matrix(design_formula, data = as.data.frame(colData_df)),
    error = function(e) NULL
  )
  if (is.null(mm)) return(FALSE)
  rk <- tryCatch(qr(mm)$rank, error = function(e) NA)
  if (is.na(rk)) return(FALSE)
  rk == ncol(mm)
}

# Reactive values for storing heatmap and volcano plot functions
significant_heatmap_reactive <- reactiveVal()
volcano_plot_reactive <- reactiveVal()

observeEvent(input$runDESeq, {
  req(wideData(), metadata(),input$designColumns, input$voi)
  
  notification_id <- showNotification("Computing, please wait...", type = "message", duration = NULL)
  
  countData <- wideData()
  metadataData <- metadata()
  
  metadataData <- metadataData %>% filter(SampleID %in% colnames(countData))
  
  countData <- countData[,metadataData$SampleID]
  metadataData$SampleID <- factor(metadataData$SampleID, levels =colnames(countData))
  metadataData <- metadataData[order(metadataData$SampleID), ]
  all(metadataData$SampleID==names(countData))
  
  if (!all(metadataData$SampleID == names(countData))) {
    showNotification("Sample IDs do not match between metadata and count data", type = "error")
    removeNotification(notification_id)  # ADDED
    return()
  }
  
  # ADDED: coerce VOI to factor and set reference if applicable
  if (input$voi %in% names(metadataData)) {
    metadataData[[input$voi]] <- factor(metadataData[[input$voi]])
    if (!is.null(input$voiRef) && input$voiRef %in% levels(metadataData[[input$voi]])) {
      metadataData[[input$voi]] <- stats::relevel(metadataData[[input$voi]], ref = input$voiRef)
    }
  }
  
  # ADDED: coerce other selected design covariates (categoricals as factors, numeric left as-is)
  other_covs <- setdiff(input$designColumns, input$voi)
  for (nm in other_covs) {
    if (is.character(metadataData[[nm]])) {
      metadataData[[nm]] <- factor(metadataData[[nm]])
    }
  }
  
  
  # CHANGED: compute and display sample sizes per covariate before modeling
  # 1) Variable of interest (VOI) group counts
  voi_counts <- table(metadataData[[input$voi]])  # CHANGED
  message("Sample sizes for variable of interest (", input$voi, "):")  # CHANGED
  print(voi_counts)  # CHANGED
  
  # Optional: show in UI
  output$voi_sample_sizes <- renderText({  # CHANGED: add a textOutput("voi_sample_sizes") in UI
    paste0("VOI sample sizes (", input$voi, "): ",
           paste(names(voi_counts), voi_counts, sep = "=", collapse = "; "))
  })
  
  # 2) Other covariates (each factor’s level counts; numeric shown as N and summary)
  cov_size_summaries <- list()  # CHANGED
  for (nm in other_covs) {  # CHANGED
    x <- metadataData[[nm]]
    if (is.factor(x)) {
      cov_size_summaries[[nm]] <- paste(names(table(x)), as.integer(table(x)), sep = "=", collapse = "; ")
      message("Sample sizes for covariate ", nm, ": ", cov_size_summaries[[nm]])  # CHANGED
    } else {
      cov_size_summaries[[nm]] <- paste0("N=", length(x),
                                         ", mean=", round(mean(as.numeric(x), na.rm = TRUE), 3),
                                         ", sd=", round(sd(as.numeric(x), na.rm = TRUE), 3))
      message("Summary for numeric covariate ", nm, ": ", cov_size_summaries[[nm]])  # CHANGED
    }
  }
  
  # Optional: show covariate summaries in UI
  output$covariate_sample_sizes <- renderUI({  # CHANGED: add a uiOutput("covariate_sample_sizes") in UI
    if (length(cov_size_summaries) == 0) return(NULL)
    tags$div(
      tags$strong("Covariate sample sizes:"),
      tags$ul(
        lapply(names(cov_size_summaries), function(nm) {
          tags$li(paste(nm, "-", cov_size_summaries[[nm]]))
        })
      )
    )
  })
  
  # CHANGED: build simple multi-covariate design formula (main effects only)
  # designFormula <- as.formula(paste("~", input$designColumn))
  designFormula <- as.formula(paste("~", paste(input$designColumns, collapse = " + ")))  # CHANGED
  
  dds <- tryCatch(
    DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                   colData = metadataData,
                                   design = designFormula),
    error = function(e) {
      removeNotification(notification_id)
      showNotification(
        paste0("Could not create DESeqDataSet: ", e$message,
               " Tip: remove redundant covariates or levels."),
        type = "error", duration = NULL
      )
      return(NULL)
    }
  )
  req(dds)  # abort if NULL
  
  dds <- tryCatch(
    DESeq2::DESeq(dds),
    error = function(e) {
      # Detect the full rank error specifically
      msg <- e$message
      if (grepl("not full rank|full rank|linear combinations", msg, ignore.case = TRUE)) {
        showNotification(
          "Model is not full rank. Remove redundant covariates or collapse levels (e.g., empty or single-sample levels).",
          type = "error", duration = NULL
        )
      } else {
        removeNotification(notification_id)
        showNotification(paste("DESeq error:", msg), type = "error")
        
      }
      return(NULL)
    }
  )
  req(dds)
  
  DESEQ_obj(dds)
  
  vst_cts <- vst(dds, blind = TRUE, nsub = sum( rowMeans( counts(dds, normalized=TRUE)) > 5 )) %>% assay()
  vst_cts_reactive(vst_cts)
  
  heatmap_annotation <- as.data.frame(colData(dds)) %>%
    dplyr::select(c(input$voi))
  heatmap_annotation_reactive(heatmap_annotation)
  #print(rownames(vst_cts))
  
  res <- DESeq2::results(dds, contrast = c(input$voi, setdiff(levels(metadataData[[input$voi]]), input$voiRef)[1], input$voiRef))  # CHANGED
  
  res <- as.data.frame(res)
  DE_miRNA_results_reactive(res)
  
  res_significant_data <- res %>% 
    filter(padj< input$DEM_padj_filter) %>%
    rownames_to_column("miRNA") %>% 
    arrange(padj)
  
  res_significant(res_significant_data)
  
  output$resSignificantTable <- DT::renderDataTable({
    req(res_significant())
    
    filtered_data <- res_significant() %>%
      dplyr::filter(padj < input$DEM_padj_filter)
    
    # Update count text
    output$dem_count <- renderText({
      paste(
        "Number of DE MiRNAs meeting the threshold criteria:",
        nrow(filtered_data)
      )
    })
    
    # Handle empty result
    if (nrow(filtered_data) == 0) {
      output$NoDEGmessage <- renderText("No DE MiRNAs meet the filtering criteria.")
      DT::datatable(
        data.frame(Message = "No significant miRNA to display"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
      return(data.frame())
    } else {
      output$NoDEGmessage <- renderText("")  # clear any prior message
      DT::datatable(
        filtered_data,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          lengthMenu = c(10, 25, 50, 100),
          autoWidth = TRUE
        )
      )
    }
  })
  
  removeNotification(notification_id)
  
  # Download handler for differential miRNA
  output$download_significant_miRNA <- downloadHandler(
    filename = function() {
      paste("Differentially_Expressed_miRNA", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(as.data.frame(res_significant()), file, row.names = FALSE, quote = FALSE)
    }
  )
  
  ### Prepare the counts for negative correlation module ---
  
  unique_groups <- unique(metadataData[[input$voi]])
  
  if (length(unique_groups) < 2) {
    showNotification("Not enough unique groups in the selected design column to compare.", type = "error")
    return()
  }
  
  gp1 <- which(metadataData[[input$voi]] == unique_groups[1])
  gp2 <- which(metadataData[[input$voi]] == unique_groups[2])
  
  p_value <- res[["pvalue"]]
  p_adjust <- res[["padj"]]
  FC <- res[["log2FoldChange"]]
  
  idx <- which(p_adjust < input$DEM_padj_filter)
  
  norm_cts <- counts(dds, normalized = TRUE) %>% as.data.frame
 
  #DE_data <- countData[idx, ]
  DE_data <- norm_cts[idx, ]
  
  mean_gp1 <- if (length(gp1) == 1) {
    norm_cts[, gp1]
  } else {
    apply(norm_cts[, gp1], 1, mean)
  }
  
  mean_gp2 <- if (length(gp2) == 1) {
    norm_cts[, gp2]
  } else {
    apply(norm_cts[, gp2], 1, mean)
  }
  
  DE_data <- cbind(DE_data, FC[idx], p_value[idx], p_adjust[idx], mean_gp1[idx], mean_gp2[idx])
  
  len_col <- ncol(DE_data)
  colnames(DE_data)[(len_col - 4):len_col] <- c("log_ratio", "P-Value", "P-adjust", "mean_case", "mean_control")
  
  FC_rows <- abs(DE_data[, len_col - 4])
  DE_data <- DE_data[FC_rows >= 0.5, ]
  
  gene_names <- row.names(DE_data)
  mirna_data <- as.data.frame(DE_data)
  row.names(mirna_data) <- gene_names
  
  mirna_data_reactive(mirna_data)
  
  # future_expansion_code: DE miRNA data table for debugging/correlation analysis
  # output$mirna_DE_table <- renderTable({
  #   req(mirna_data_reactive())
  #   as.data.frame(mirna_data_reactive())
  # })
  # 
  # # Download handler for miRNA DE table (for correlation analysis)
  # output$download_mirna_DE_table <- downloadHandler(
  #   filename = function() {
  #     paste("miRNA_DE_data_for_correlation_", Sys.Date(), ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     write.csv(as.data.frame(mirna_data_reactive()), file, row.names = TRUE)
  #   }
  # )
  
  ### Show a Select DE miRNA as a Boxplot ---
  
  updateSelectInput(session, "selectedMiRNA", choices = res_significant()$miRNA)
  updateSelectInput(session, "DEMiRNA", choices = res_significant()$miRNA)
  
  grouping_variable <- input$voi
  levels_group <- levels(as.factor(metadataData[[grouping_variable]]))
  
  if (length(levels_group) < 2) {
    showNotification("Not enough levels in the selected grouping variable for comparison.", type = "error")
    return()
  }
  
  output$boxplot <- renderPlot({
    req(input$groupBy_Boxplot, 
        input$selectedMiRNA, 
        DESEQ_obj())
    grouping_variable <- input$groupBy_Boxplot
    
    dds <- DESEQ_obj()
    selected_miRNA <- input$selectedMiRNA
    
    create_boxplot(dds, selected_miRNA, grouping_variable)
  })
  
  output$downloadBoxplot <- downloadHandler(
    filename = function() {
      fmt <- input$plot_format
      paste0("boxplot_", input$selectedMiRNA, ".", fmt)
    },
    content = function(file) {
      open_graphics_device(
        file,
        format    = input$plot_format,
        width_px  = input$plot_width,
        height_px = input$plot_height,
        res       = 96
      )
      dds <- DESEQ_obj()
      selected_miRNA <- input$selectedMiRNA
      grouping_variable <- input$groupBy_Boxplot
      
      print(create_boxplot(dds, selected_miRNA, grouping_variable))
      dev.off()
    }
  )
  
  output$download_plotCounts_selected <- downloadHandler(
    filename = function() paste0("plotCounts_", input$selectedMiRNA, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(DESEQ_obj(), input$selectedMiRNA, input$groupBy_Boxplot)
      df <- DESeq2::plotCounts(
        dds = DESEQ_obj(),
        gene = input$selectedMiRNA,
        intgroup = input$groupBy_Boxplot,
        returnData = TRUE
      )
      utils::write.csv(df, file, row.names = FALSE, quote = FALSE)
    }
  )
  
  
  
  # Data processing msg:
  output$boxplotMessage <- renderUI({
    req(DESEQ_obj(), metadata(), input$groupBy_Boxplot, input$selectedMiRNA)
    
    dds <- DESEQ_obj()
    metadataData <- metadata()
    grouping_variable <- input$groupBy_Boxplot
    selected_miRNA <- input$selectedMiRNA
    
    # Get per-sample values shown in the boxplot (same as plotCounts with returnData=TRUE)
    df <- DESeq2::plotCounts(
      dds = dds,
      gene = selected_miRNA,
      intgroup = grouping_variable,
      returnData = TRUE
    )
    
    # Group sample sizes for the selected grouping variable
    group_counts <- df %>%
      dplyr::count(!!rlang::sym(grouping_variable), name = "n")
    
    group_text <- paste(
      paste0(group_counts[[grouping_variable]], " (n = ", group_counts$n, ")"),
      collapse = " | "
    )
    
    # Compose an easy-to-understand description
    transform_text <- HTML(paste0(
      "<b> Data processing:</b><br/>",
      "1) Counts are first adjusted by DESEQ2 size factors to correct for sequencing depth and sample composition.<br/>",
      "2) The distribution of normalized counts for <code>", selected_miRNA, "</code> across samples are shown.<br/>",
      "- For statistical significance and effect size, please use the DESeq2 differential expression results.<br/>"
      
    ))
    
    HTML(paste0(
      "<div>",
      "<b>Selected miRNA:</b> ", selected_miRNA, "<br/>",
      "<b>Group sample sizes:</b> ", group_text, "<br/><br/>",
      transform_text, "<br/>",
      "</div>"
    ))
  })
  
})



# Populate the selectable miRNA list from normalized counts (safer for heatmap rows)
observeEvent(vst_cts_reactive(), {
  req(vst_cts_reactive())
  vst_cts <- vst_cts_reactive()
  updateSelectizeInput(
    session,
    inputId = "heatmap_select_rows",
    choices = rownames(vst_cts),
    server = TRUE
  )
})

observeEvent(input$significant_heatmap_button, {
  req(vst_cts_reactive(), res_significant(),heatmap_annotation_reactive())
  
  res_significant_data <- res_significant()
  vst_cts <- vst_cts_reactive()
  heatmap_annotation <- heatmap_annotation_reactive()
  
  # If user specified rows, use them; otherwise fall back to "significant" logic
  user_selected <- input$heatmap_select_rows
  valid_user_selected <- intersect(user_selected %||% character(0), rownames(vst_cts))
  
  # Warn about any invalid entries
  if (!is.null(user_selected) && length(setdiff(user_selected, valid_user_selected)) > 0) {
    output$heatmap_row_warning <- renderText(
      paste("Some miRNAs not found in normalized counts:",
            paste(setdiff(user_selected, valid_user_selected), collapse = ", "))
    )
  } else {
    output$heatmap_row_warning <- renderText("")
  }
  
  # Determine which rows to plot
  if (length(valid_user_selected) > 0) {
    selected_present <- valid_user_selected
  } else {
    # Fall back to significant-based selection when no user override provided
    if (is.null(res_significant_data) || nrow(res_significant_data) == 0) {
      output$significant_miRNA_heatmap <- renderPlot({
        plot.new(); text(0.5, 0.5, "no significant features to plot")
      })
      significant_heatmap_reactive(NULL)
      return()
    }
    
    if ("log2FoldChange" %in% colnames(res_significant_data)) {
      select_significant <- res_significant_data %>%
        dplyr::mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
        dplyr::filter(abs_log2FoldChange >= input$DEM_log2FoldChange_filter) %>%
        dplyr::pull(miRNA)
    } else if ("logFC" %in% colnames(res_significant_data)) {
      select_significant <- res_significant_data %>%
        dplyr::mutate(abs_log2FoldChange = abs(logFC)) %>%
        dplyr::filter(abs_log2FoldChange >= input$DEM_log2FoldChange_filter) %>%
        dplyr::pull(miRNA)
    } else {
      select_significant <- res_significant_data$miRNA
    }
    
    if (length(select_significant) == 0) {
      output$significant_miRNA_heatmap <- renderPlot({
        plot.new(); text(0.5, 0.5, "no significant features to plot")
      })
      significant_heatmap_reactive(NULL)
      return()
    }
    
    selected_present <- intersect(rownames(vst_cts), select_significant)
  }
  
  # Guard for empty selection
  if (length(selected_present) == 0) {
    output$significant_miRNA_heatmap <- renderPlot({
      plot.new(); text(0.5, 0.5, "no matching miRNAs to plot")
    })
    significant_heatmap_reactive(NULL)
    return()
  }
  
  subset_vst_cts <- vst_cts[selected_present, , drop = FALSE]
  
  do_cluster_rows <- nrow(subset_vst_cts) >= 2
  
  heatmap_significant <- pheatmap::pheatmap(
    subset_vst_cts,
    annotation_col = heatmap_annotation,
    scale = "row",
    cluster_rows = do_cluster_rows,
    show_rownames = input$significant_heatmap_showRowLabels,
    show_colnames = input$significant_heatmap_showColLabels,
    cluster_cols = TRUE
  )
  
  output$significant_miRNA_heatmap <- renderPlot({ heatmap_significant })
  
  significant_heatmap_reactive(heatmap_significant)
  
  # NEW 2026 enable data download 
  subset_vst_cts_scaled <- t(scale(t(subset_vst_cts)))
  heatmap_mat_plotted_reactive(as.data.frame(subset_vst_cts_scaled))
  
  output$download_heatmap_matrix <- downloadHandler(
    filename = function() paste0("significant_miRNA_heatmap_matrix_as_plotted_", Sys.Date(), ".csv"),
    content = function(file) {
      req(heatmap_mat_plotted_reactive())
      out <- heatmap_mat_plotted_reactive() |> tibble::rownames_to_column("miRNA")
      write.csv(out, file, row.names = FALSE, quote = FALSE)
    }
  )

})

  


# After DE_miRNA_results_reactive() becomes available, populate label choices
observeEvent(DE_miRNA_results_reactive(), {
  req(DE_miRNA_results_reactive())
  res <- DE_miRNA_results_reactive()
  updateSelectizeInput(
    session,
    inputId = "volcano_select_labels",
    choices = rownames(res),
    server = TRUE
  )
})

observeEvent(input$Volcano_Plot_Button, {  
  req(DE_miRNA_results_reactive())

  res <- DE_miRNA_results_reactive()
  
  # Validate user selections against available rownames
  selected_labels <- input$volcano_select_labels
  valid_labels <- intersect(selected_labels %||% character(0), rownames(res))
  
  # Warn if any entered labels aren't found
  if (!is.null(selected_labels) && length(setdiff(selected_labels, valid_labels)) > 0) {
    output$volcano_label_warning <- renderText(
      paste("Some labels were not found in results:",
            paste(setdiff(selected_labels, valid_labels), collapse = ", "))
    )
  } else {
    output$volcano_label_warning <- renderText("")
  }
  
  # Toggle: NULL for all labels, vector for only selected labels
  select_lab_arg <- if (length(valid_labels) == 0) NULL else valid_labels
  
  output$volcanoPlot <- renderPlot({
    EnhancedVolcano(
      res,
      lab = rownames(res),
      x = "log2FoldChange",
      y = "padj",
      ylab = bquote(~-Log[10]~ "(p-adjusted)"),
      pCutoff = input$DEM_padj_filter,
      FCcutoff = input$DEM_log2FoldChange_filter,
      labSize = input$volcano_plot_label_size,
      selectLab = select_lab_arg,   # core logic
      max.overlaps = Inf,
    ) +
      theme_cowplot(12) +
      theme(
        legend.position = "top",
        legend.justification = "center",
        plot.subtitle = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5)
      )
  })
  
  # Keep the stored plot in sync (for export)
  volcano_plot_reactive(function() {
    EnhancedVolcano(
      res,
      lab = rownames(res),
      x = "log2FoldChange",
      y = "padj",
      ylab = bquote(~-Log[10]~ "(p-adjusted)"),
      pCutoff = input$DEM_padj_filter,
      FCcutoff = input$DEM_log2FoldChange_filter,
      labSize = input$volcano_plot_label_size,
      selectLab = select_lab_arg,
      max.overlaps = Inf,
      drawConnectors = TRUE
    ) +
      theme_cowplot(12) +
      theme(
        legend.position = "top",
        legend.justification = "center",
        plot.subtitle = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5)
      )
  })
})

### Do Differential miRNA Expression Analysis with limma ---

# future_expansion_code: limma::voom analysis
# observeEvent(input$runVoom, {
#   req(wideData(), 
#       metadata(), 
#       input$designColumn)
#   
#   countData <- wideData()
#   metadataData <- metadata()
#   
#   metadataData <- metadataData %>% filter(SampleID %in% colnames(countData))
#   countData <- countData[, metadataData$SampleID]
#   metadataData$SampleID <- factor(metadataData$SampleID, levels = colnames(countData))
#   metadataData <- metadataData[order(metadataData$SampleID), ]
#   
#   countData <- DGEList(counts=countData,genes=rownames(countData))
#   
#   design <- model.matrix(~ metadataData[[input$designColumn]])
#   
#   v <- voom(countData, design, plot = TRUE)
#   
#   fit <- lmFit(v, design)
#   
#   fit <- eBayes(fit)
#   
#   res <- topTable(fit, adjust = "BH", sort.by = "P", number = Inf)
#   
#   res_significant_data <- res %>% 
#     filter(adj.P.Val < 0.05) %>%
#     rownames_to_column("miRNA")
#   
#   res_significant(res_significant_data)
#   
#   output$resSignificantTable <- renderTable({
#     req(res_significant())
#     as.data.frame(res_significant())
#   })
#   
#   updateSelectInput(session, "selectedMiRNA", choices = res_significant()$miRNA)
#   
#   output$volcanoPlot <- renderPlot({
#     EnhancedVolcano(res,
#                     lab = rownames(res),
#                     labSize = 6,
#                     x = 'logFC',
#                     y = 'adj.P.Val',
#                     ylab = bquote(~-Log[10]~ '(adjusted p-value)'),
#                     pCutoff = 0.05,
#                     title = "Differential Expression Analysis",
#                     subtitle = "limma::voom",
#                     FCcutoff = 1) +
#       theme_cowplot(12) +
#       theme(legend.position = "top",
#             plot.subtitle = element_text(hjust = 0.5),
#             legend.justification = "center") +
#       theme(plot.title = element_text(hjust = 0.5))
#   })
#   
#   # Download handler for limma volcano plot
#   output$downloadVolcanoLimma <- downloadHandler(
#     filename = function() {
#       paste("Volcano_Plot_limma_voom_", Sys.Date(), ".png", sep = "")
#     },
#     content = function(file) {
#       png(file, width = 800, height = 600)
#       print(
#         EnhancedVolcano(res,
#                         lab = rownames(res),
#                         labSize = 6,
#                         x = 'logFC',
#                         y = 'adj.P.Val',
#                         ylab = bquote(~-Log[10]~ '(adjusted p-value)'),
#                         pCutoff = 0.05,
#                         title = "Differential Expression Analysis",
#                         subtitle = "limma::voom",
#                         FCcutoff = 1) +
#           theme_cowplot(12) +
#           theme(legend.position = "top",
#                 plot.subtitle = element_text(hjust = 0.5),
#                 legend.justification = "center") +
#           theme(plot.title = element_text(hjust = 0.5))
#       )
#       dev.off()
#     }
#   )
# })
