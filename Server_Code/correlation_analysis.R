# Module 7: miRNA-mRNA Correlation Analysis

mrna_no_zero_reactive <- reactive({
  req(mrna_data_reactive())
  mrna_data <- as.data.frame(mrna_data_reactive())
  last_five_columns <- (ncol(mrna_data) - 4):ncol(mrna_data)
  mrna_data[rowSums(mrna_data[, -last_five_columns, drop = FALSE]) != 0, ] %>% 
    filter(abs(log_ratio) >= input$filter_DEG_correlation)
})

mirna_no_zero_reactive <- reactive({
  req(mirna_data_reactive())
  mirna_data <- as.data.frame(mirna_data_reactive())
  last_five_columns <- (ncol(mirna_data) - 4):ncol(mirna_data)
  mirna_data[rowSums(mirna_data[, -last_five_columns, drop = FALSE]) != 0, ] %>% 
    filter(abs(log_ratio) >= input$filter_DEM_correlation)
})

output$DEG_DEM_filter_message <- renderText({
  mrna_no_zero <- mrna_no_zero_reactive()
  mirna_no_zero <- mirna_no_zero_reactive()
  
  num_mrna_remaining <- nrow(mrna_no_zero)
  num_mirna_remaining <- nrow(mirna_no_zero)
  
  paste("This will correlate", num_mrna_remaining, "DEGs with",
        num_mirna_remaining, "DE miRNAs")
})

observeEvent(input$runCorrelationAnalysis, {
  if (is.null(mirna_data_reactive()) || is.null(mrna_data_reactive())) {
    showNotification(
      HTML("Please complete the following steps first:<br>
            1. Run Module 3 (Differential miRNA Expression Analysis)<br>
            2. Run Module 6 (Differential Gene Expression Analysis)<br>
            Both analyses must be completed before running correlation analysis."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  req(mirna_no_zero_reactive(), mrna_no_zero_reactive())
  
  mrna_no_zero <- mrna_no_zero_reactive()
  mirna_no_zero <- mirna_no_zero_reactive()
  
  cor <- negative_cor(mrna_data = mrna_no_zero, 
                      mirna_data = mirna_no_zero, 
                      cut.off = 0, 
                      method = "spearman")
  cor <- as.data.frame(cor)
  cor$Correlation <- as.numeric(cor$Correlation)
  

  cor$P_adjust_corr <- stats::p.adjust(cor$`P-value(correlation)`)
 
  # Store full correlation results for later use
  cor_results_reactive(cor)
  
  neg_cor <- cor
  neg_cor <- dplyr::rename(neg_cor,ensembl_gene_id = Gene)
  neg_cor$ensembl_gene_id <- gsub("\\..*","",neg_cor$ensembl_gene_id)
  
  hsmart <- mart()
  
  mapping <- getBM(
    attributes = c('ensembl_gene_id',org_string(),'description', 'entrezgene_id'), 
    filters = 'ensembl_gene_id',
    values = neg_cor$ensembl_gene_id,
    mart = hsmart
  )
  
  neg_cor <- merge(neg_cor, mapping, by="ensembl_gene_id")
  neg_cor_reactive(neg_cor)
  
  output$num_neg_correlations <- renderText({
    if (!is.null(neg_cor)) {
      paste("The number of significant miRNA-gene correlations is",nrow(neg_cor %>% filter(P_adjust_corr < input$corr_padj_cutoff)))
    } else {
      paste("No miRNA-gene correlations to show!")
      return()
    }
  })
  
  output$NegativeCorrelationTable <- renderDT({
    if (!is.null(neg_cor)) {
      datatable(neg_cor %>%
                  filter(P_adjust_corr < input$corr_padj_cutoff), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    } else {
      return()
    }
  })
  
  # Static Table Version for Paper
  # output$NegativeCorrelationTable <- renderTable({
  #   if (!is.null(neg_cor)) {
  #     head(neg_cor %>% filter(P_adjust_corr < input$corr_padj_cutoff))
  #   } else {
  #     return()
  #   }
  # })
  
  # Clear any previous permutation result 
  global_cor_perm_result(NULL)
  
})

############################# 2026 NEWLY ADDED ###########################################

observeEvent(input$runGlobalPermutation, {
  if (is.null(cor_results_reactive())) {
    showNotification("Please run the correlation analysis first.", type = "warning")
    return()
  }
  
  req(cor_results_reactive(), mrna_no_zero_reactive(), mirna_no_zero_reactive())
  
  cor          <- neg_cor_reactive()
  mrna_no_zero <- mrna_no_zero_reactive()
  mirna_no_zero <- mirna_no_zero_reactive()
  
  # Observed statistic
  sig_neg <- cor %>%
    dplyr::filter(Correlation < 0, P_adjust_corr < input$corr_padj_cutoff)
  observed_stat <- nrow(sig_neg)
  
  n_perm <- input$n_perm_correlation
  
  # Expression matrices 
  ncol_mrna <- ncol(mrna_no_zero)
  ncol_mir  <- ncol(mirna_no_zero)
  
  mrna_expr  <- as.matrix(mrna_no_zero[, 1:(ncol_mrna - 5), drop = FALSE])
  mirna_expr <- as.matrix(mirna_no_zero[, 1:(ncol_mir - 5), drop = FALSE])
  
  set.seed(123)
  perm_stats <- numeric(n_perm)
  
  withProgress(message = "Running global permutation test...", value = 0, {
    for (b in seq_len(n_perm)) {
      incProgress(1 / n_perm, detail = paste("Permutation", b, "of", n_perm))
      
      perm_idx <- sample(ncol(mirna_expr))
      mirna_perm <- mirna_expr[, perm_idx, drop = FALSE]
      
      mrna_perm_df <- mrna_no_zero
      mirna_perm_df <- mirna_no_zero
      mrna_perm_df[, 1:(ncol_mrna - 5)]  <- mrna_expr
      mirna_perm_df[, 1:(ncol_mir - 5)]  <- mirna_perm
      
      cor_perm <- negative_cor(
        mrna_data = mrna_perm_df,
        mirna_data = mirna_perm_df,
        cut.off    = 0,
        method     = "spearman"
      )
      cor_perm <- as.data.frame(cor_perm)
      cor_perm$Correlation   <- as.numeric(cor_perm$Correlation)
      cor_perm$P_adjust_corr <- stats::p.adjust(cor_perm$`P-value(correlation)`, method = "BH")
      
      sig_neg_perm <- cor_perm %>%
        dplyr::filter(Correlation < 0, P_adjust_corr < input$corr_padj_cutoff)
      
      perm_stats[b] <- nrow(sig_neg_perm)
    }
  })
  
  global_p_perm <- (sum(perm_stats >= observed_stat) + 1) / (n_perm + 1)
  
  global_cor_perm_result(list(
    observed_stat = observed_stat,
    perm_stats    = perm_stats,
    n_perm        = n_perm,
    p_value       = global_p_perm
  ))
  
  output$GlobalPermutationSummary <- renderText({
    res <- global_cor_perm_result()
    if (is.null(res)) return("Run the global permutation test to see results.")
    
    interp <- if (res$p_value < 0.05) {
      "There is a statistically significant enrichment of negative miRNA–mRNA correlations beyond random expectation."
    } else {
      "The overall number of negative miRNA–mRNA correlations could be explained by random association; individual pairs may still be of interest."
    }
    
    null_vals <- res$perm_stats
    null_min  <- min(null_vals)
    null_max  <- max(null_vals)
    null_mean <- mean(null_vals)
    null_med  <- stats::median(null_vals)
    
    paste0(
      "Global permutation test (", res$n_perm, " permutations)\n",
      "Observed number of significantly negative miRNA–gene correlations ",
      "(Correlation < 0 & BH-adjusted p < ", input$corr_padj_cutoff, "): ",
      res$observed_stat, "\n",
      "Null distribution (permuted data): min = ", null_min,
      ", max = ", null_max,
      ", mean = ", signif(null_mean, 3),
      ", median = ", null_med, "\n",
      "Global permutation p-value: ", signif(res$p_value, 3), "\n",
      "Interpretation: ", interp
    )
  })
  
  
  
  
  
})

#############################################################################

observeEvent(input$show_database_support, {
  
  
  if (is.null(neg_cor_reactive()) || is.null(predicted_target_reactive())) {
    showNotification(
      HTML("Please complete the following steps first:<br>
            1. Run 'Retrieve Targeted Genes for all DE miRNA' in the 'All DE miRNA Target Retrieval and Pathway Enrichment' tab. <br>
            2. Run 'Correlation Analysis' in the current tab. <br>
            These analyses must be completed before showing which of the miRNA-gene pairs are supported by existing databases."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  req(neg_cor_reactive(),predicted_target_reactive())
  
  predicted <- as.data.frame(predicted_target_reactive()) %>%
    as_tibble() %>%
    dplyr::select(ensembl_gene_id,!!sym(org_string()),miRNA, entrezgene_id) %>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id))
  neg_cor_pairs <- as.data.frame(neg_cor_reactive()) %>% 
    as_tibble() %>%
    dplyr::select(ensembl_gene_id,!!sym(org_string()), miRNA, entrezgene_id)%>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id))
  
  identical_rows <- inner_join(predicted,neg_cor_pairs, relationship = "many-to-many") %>% unique()
  supported_neg_cor_reactive(identical_rows)
  
  if ( nrow(identical_rows) == 0) {
    output$degMessage <- renderText("None of the miRNA-gene pairs are supported by database predictions")
  } else {
    output$degMessage <- renderText(paste("The number of supported miRNA-gene pairs is", nrow(identical_rows)))
  }
  
  output$SupportedNegativeCorrelationTable <- renderDT({
    if (!is.null(identical_rows)) {
      datatable(identical_rows, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    } else {
      return()
    }
  })
  
  # Static table version for paper 
  # output$SupportedNegativeCorrelationTable <- renderTable({
  #   if (!is.null(identical_rows)) {
  #    head(identical_rows)
  #   } else {
  #     return()
  #   }
  # })
  
  
  
})

observeEvent(input$CorrelationPlotButton, {
  
  plot_notif <- showNotification(HTML("Plotting, please wait..."),type = "warning")
  
  req(neg_cor_reactive(),org_string())
  
  if (is.null(neg_cor_reactive()) || is.null(org_string())) {
    showNotification(
      HTML("Please complete the following steps first:<br>
            1. Run 'Retrieve Targeted Genes for all DE miRNA' in the 'All DE miRNA Target Retrieval and Pathway Enrichment tab. <br>
            2. Run 'Correlation Analysis' in the current tab. <br>
            These analyses must be completed first."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  neg_cor <- neg_cor_reactive()

  
  ######### DEV CODE ############
  #neg_cor <- read.csv(here("dev_files/neg_cor_2025-11-14.csv"))
  # supported_neg_cor <- read.csv(here("dev_files/Supported_Negative_Correlations_2025-11-14_.csv"))
  ###############################

  
  # Optionally restrict to supported interactions only
  if (isTRUE(input$select_supported_only)) {
    supported_neg_cor <- supported_neg_cor_reactive()  # get current value from reactiveVal
    
    if (is.null(supported_neg_cor) || nrow(supported_neg_cor) == 0) {
      showNotification(
        "Supported interactions data is empty or unavailable; skipping supported-only filter.",
        type = "error",
        duration = 5
      )
      neg_cor <- neg_cor
    } else {
      neg_cor <- neg_cor %>%
        dplyr::semi_join(supported_neg_cor, by = c("ensembl_gene_id", org_string(), "miRNA"))
    }
  } else {
    neg_cor <- neg_cor
  }
  
  print("successfully joined plot")

  
  miRNAs_with_high_correlation <- neg_cor %>%
    mutate(abs_Correlation = abs(Correlation)) %>%
    filter(abs_Correlation >= abs(input$correlationCutoff)) %>%
    filter(P_adjust_corr < input$correlation_padj_Cutoff) %>%
    mutate(abs_logratio_miRNA = abs(as.numeric(logratio_miRNA))) %>%
    filter(abs_logratio_miRNA >= input$logratio_miRNA_Cutoff) %>%
    pull(miRNA) %>%
    unique()
  
  genes_filtered <- neg_cor %>%
    mutate(abs_logratio_gene = abs(as.numeric(logratio_gene))) %>%
    filter(abs_logratio_gene >= input$logratio_gene_Cutoff) %>%
    pull(!!sym(org_string()))
  
  filtered_neg_cor <- neg_cor %>%
    filter(miRNA %in% miRNAs_with_high_correlation) %>% 
    filter(!!sym(org_string()) %in% genes_filtered)
  
  mean_correlations <- filtered_neg_cor %>%
    group_by(miRNA) %>%
    summarise(Mean_Correlation = mean(Correlation)) %>%
    arrange(Mean_Correlation)
  
  filtered_neg_cor$miRNA <- factor(filtered_neg_cor$miRNA, levels = mean_correlations$miRNA)
  
  print(names(filtered_neg_cor))
  
  plot <- ggplot(filtered_neg_cor,
                 aes(
                   x = !!sym(org_string()),
                   y = miRNA,
                   size = abs(Correlation),
                   color = Correlation,
                   text = paste0(
                     "miRNA: ", miRNA,
                     "<br>Gene: ", !!sym(org_string()),
                     "<br>Correlation: ", round(Correlation, 3),
                     "<br>description: ", description
                   )
                 )) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "Negative Correlation between miRNA and Gene",
         x = "Gene",
         y = "miRNA",
         size = "Correlation (absolute value)",
         color = "Correlation Value") +
    theme_cowplot(12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  removeNotification(plot_notif)
  
  output$correlationPlot <- plotly::renderPlotly({       
    plotly::ggplotly(plot, tooltip = c("x", "y", "color", "size"))  
  })
  
  # correlation_plot_reactive(plot)                     
  correlation_plot_reactive(plot)      

})


# Data processing details 
output$correlationMessage <- renderUI({
  req(mirna_data_reactive(), mrna_data_reactive(),
      input$filter_DEM_correlation, input$filter_DEG_correlation,
      input$correlationCutoff, input$correlation_padj_Cutoff,
      input$logratio_miRNA_Cutoff, input$logratio_gene_Cutoff)
  
  # Base inputs
  mirna_raw <- as.data.frame(mirna_data_reactive())
  mrna_raw  <- as.data.frame(mrna_data_reactive())
  
  # Count how many features enter from the significant sets
  n_mirna_sig <- nrow(mirna_raw)
  n_mrna_sig  <- nrow(mrna_raw)
  
  # Apply your zero-row filtering and log-ratio prefilters (mirror reactives)
  last5_mi <- (ncol(mirna_raw) - 4):ncol(mirna_raw)
  last5_mr <- (ncol(mrna_raw)  - 4):ncol(mrna_raw)
  
  mirna_kept <- mirna_raw[rowSums(mirna_raw[, -last5_mi, drop = FALSE]) != 0, ]
  mirna_kept <- dplyr::filter(mirna_kept, abs(log_ratio) >= input$filter_DEM_correlation)
  
  mrna_kept <- mrna_raw[rowSums(mrna_raw[, -last5_mr, drop = FALSE]) != 0, ]
  mrna_kept <- dplyr::filter(mrna_kept, abs(log_ratio) >= input$filter_DEG_correlation)
  
  n_mirna_pre  <- nrow(mirna_kept)
  n_mrna_pre   <- nrow(mrna_kept)
  
  # Build a concise description of downstream filtering used in the plot step
  downstream_filters <- paste0(
    "- Correlation filter: |Spearman r| ≥ ", signif(abs(input$correlationCutoff), 3), ".<br/>",
    "- Multiple-testing: adjusted p-value (FDR) ≤ ", signif(input$correlation_padj_Cutoff, 3), ".<br/>",
    "- miRNA effect size: |log ratio| ≥ ", signif(input$logratio_miRNA_Cutoff, 3), " (for plotting).<br/>",
    "- Gene effect size: |log ratio| ≥ ", signif(input$logratio_gene_Cutoff, 3), " (for plotting).<br/>",
    if (isTRUE(input$select_supported_only)) {
      "- Supported only: keeping pairs present in external evidence (annotation list).<br/>"
    } else { "" }
  )
  
  HTML(paste0(
    "<div>",
    "<b>Data processing:</b><br/>",
    "1) The correlation analysis starts from DESeq2 normalized counts for features marked as significant ",
    "and which also have |log2 fold change| ≥ 0.5.<br/>",
    "2) Apply prefilter on effect size: keep features with |log ratio| ≥ the user-selected thresholds ",
    "(miRNA: ", signif(input$filter_DEM_correlation, 3), "; gene: ", signif(input$filter_DEG_correlation, 3), ").<br/>",
    "3) Compute Spearman correlations between each remaining miRNA and mRNA across matched samples.<br/>",
    "4) Adjust correlation p-values for multiple testing (FDR).<br/><br/>",
    
    "<b>Current dataset sizes:</b><br/>",
    "- Significant miRNAs input: ", n_mirna_sig, "; after prefilter: ", n_mirna_pre, ".<br/>",
    "- Significant genes input: ", n_mrna_sig,  "; after  prefilter: ", n_mrna_pre,  ".<br/><br/>",
    
    "<b>Plot-level filters applied:</b><br/>",
    downstream_filters,
    "<br/>",
    
    "<b>Interpretation notes:</b><br/>",
    "- DESeq2 size-factor normalization reduces depth/composition effects before ranking (Spearman) which helps address technical artifacts.<br/>",
    "- Prefiltering small effects stabilizes ranks and reduces spurious correlations.<br/>",
    "</div>"
  ))
})



observeEvent(input$showGeneNetwork_Experimental, {
  req(neg_cor_reactive())
  
  if (is.null(neg_cor_reactive()) || is.null(org_string())) {
    showNotification(
      HTML("Please complete the following steps first:<br>
            1. Run 'Retrieve Targeted Genes for all DE miRNA' in the 'All DE miRNA Target Retrieval and Pathway Enrichment' tab. <br>
            2. Run 'Correlation Analysis' in the current tab. <br>
            These analyses must be completed first."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  neg_cor <- neg_cor_reactive()
  
  
  ################ DEV CODE #####################
  
  # neg_cor <- readr::read_csv(here("dev_files/neg_cor_2025-11-14.csv"))
  # supported_neg_cor <- readr::read_csv(here("dev_files/Supported_Negative_Correlations_2025-11-14_.csv"))
  # 
  # neg_cor <- neg_cor %>%
  #   semi_join(supported_neg_cor, by = c("ensembl_gene_id", "hgnc_symbol", "miRNA"))
  # filtered_neg_cor <- neg_cor
  ###############################################
  
  # Optionally restrict to supported interactions only
  if (isTRUE(input$select_supported_only)) {
    supported_neg_cor <- supported_neg_cor_reactive()
    
    if (is.null(supported_neg_cor) || nrow(supported_neg_cor) == 0) {
      showNotification(
        "Supported interactions data is empty or unavailable; skipping supported-only filter.",
        type = "error",
        duration = 5
      )
      neg_cor <- neg_cor
    } else {
      neg_cor <- neg_cor %>%
        dplyr::semi_join(supported_neg_cor, by = c("ensembl_gene_id", org_string(), "miRNA"))
    }
  } else {
    neg_cor <- neg_cor
  }
  
  # Filter by thresholds similar to your plot
  miRNAs_with_high_correlation <- neg_cor %>%
    mutate(abs_Correlation = abs(Correlation)) %>%
    filter(abs_Correlation >= abs(input$correlationCutoff)) %>%
    filter(P_adjust_corr < input$correlation_padj_Cutoff) %>%
    mutate(abs_logratio_miRNA = abs(as.numeric(logratio_miRNA))) %>%
    filter(abs_logratio_miRNA >= input$logratio_miRNA_Cutoff) %>%
    pull(miRNA) %>%
    unique()
  
  genes_filtered <- neg_cor %>%
    mutate(abs_logratio_gene = abs(as.numeric(logratio_gene))) %>%
    filter(abs_logratio_gene >= input$logratio_gene_Cutoff) %>%
    pull(!!sym(org_string()))
  
  filtered_neg_cor <- neg_cor %>%
    filter(miRNA %in% miRNAs_with_high_correlation) %>%
    filter(!!sym(org_string()) %in% genes_filtered)
  
  # If no edges, show notification and clear output
  if (nrow(filtered_neg_cor) == 0) {
    showNotification("No interactions match the selected thresholds.", type = "warning", duration = 6)
    output$miRNAGeneNetwork <- visNetwork::renderVisNetwork({
      visNetwork::visNetwork(data.frame(), data.frame())
    })
    # Store empty for downloads
    vis_nodes_reactive(NULL)
    vis_edges_reactive(NULL)
    return()
  }
  
  
  # Build edges (unique pairs are recommended to avoid inflated degree)
  edges <- filtered_neg_cor %>%
    dplyr::filter(!is.na(!!sym(org_string())) & !!sym(org_string()) != "") %>%
    dplyr::transmute(
      from = miRNA,
      to   = !!sym(org_string()),
      title = paste0(
        "<b>miRNA:</b> ", miRNA,
        "<br><b>Gene:</b> ", !!sym(org_string()),
        ifelse(!is.na(Correlation), paste0("<br><b>Correlation:</b> ", Correlation), ""),
        ifelse(!is.na(P_adjust_corr), paste0("<br><b>Adjusted p-value:</b> ", P_adjust_corr), ""),
        ifelse(!is.na(logratio_miRNA), paste0("<br><b>Log-Ratio of miRNA:</b> ", logratio_miRNA), ""),
        ifelse(!is.na(logratio_gene), paste0("<br><b>Log-Ratio of Gene:</b> ", logratio_gene), "")
      )
    ) %>%
    dplyr::distinct(from, to, .keep_all = TRUE)  # ensure unique edges
  

  up_mirnas   <- tryCatch(upregulated_miRNA_reactive() %>% pull(miRNA) %>% unique(), error = function(e) NULL)
  down_mirnas <- tryCatch(downregulated_miRNA_reactive() %>% pull(miRNA) %>% unique(), error = function(e) NULL)
  
  # Apply miRNA-direction filtering
  if (!is.null(input$corr_miRNA_dir)) {
    if (input$corr_miRNA_dir == "up") {
      edges <- edges %>% dplyr::filter(from %in% up_mirnas)
    } else if (input$corr_miRNA_dir == "down") {
      edges <- edges %>% dplyr::filter(from %in% down_mirnas)
    } else {
      # "both": no restriction by direction
      edges <- edges
    }
  }
  
  # Early exit if nothing remains
  validate(need(nrow(edges) > 0, "No interactions after filtering by miRNA direction."))
  
  
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
  if (!is.null(input$min_degree_experimental) && input$min_degree_experimental > 0) {
    keep_ids <- nodes %>%
      dplyr::filter(degree >= input$min_degree_experimental) %>%
      dplyr::pull(id)
    edges <- edges %>% dplyr::filter(from %in% keep_ids, to %in% keep_ids)
    nodes <- nodes %>% dplyr::filter(id %in% keep_ids)
  }
  
  # If graph still too large, keep top N nodes by preliminary degree
  if (nrow(nodes) > input$max_nodes_experimental) {
    keep <- nodes %>%
      dplyr::arrange(dplyr::desc(degree)) %>%
      dplyr::slice_head(n = input$max_nodes_experimental) %>%
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
  corr_network_displayed(list(nodes = nodes, edges = edges))
  
  output$miRNAGeneNetwork_Experimental <- visNetwork::renderVisNetwork({
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
  
  ################################# NEW 2026 ######################################
  output$corr_dl_network_nodes <- downloadHandler(
    filename = function() paste0("corr_displayed_network_nodes_", Sys.Date(), ".csv"),
    content = function(file) {
      req(corr_network_displayed()$nodes)
      readr::write_csv(corr_network_displayed()$nodes, file)
    }
  )
  
  output$corr_dl_network_edges <- downloadHandler(
    filename = function() paste0("corr_displayed_network_edges_", Sys.Date(), ".csv"),
    content = function(file) {
      req(corr_network_displayed()$edges)
      readr::write_csv(corr_network_displayed()$edges, file)
    }
  )
  ##################################################################################
})


observeEvent(input$Pathway_ORA_Experimental, {
  
  if (is.null(neg_cor_reactive()) || is.null(supported_neg_cor_reactive())) {
    showNotification(
      HTML("Please complete the following steps first:<br>
            1. Run 'Retrieve Targeted Genes for all DE miRNA' in the 'All DE miRNA Target Retrieval and Pathway Enrichment' tab. <br>
            2. Run 'Correlation Analysis' in the current tab. <br>
            3. Run 'See Correlations with Database Support' in the current tab. <br>
            These analyses must be completed first."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  req(mrna_data_reactive(),upregulated_miRNA_reactive(),downregulated_miRNA_reactive(), supported_neg_cor_reactive(), neg_cor_reactive())
  
  
  pathway_notification <- showNotification("Looking for enriched pathways. This may take some time", type = "message", duration = NULL)

  # Load required items
  neg_cor <- neg_cor_reactive() %>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id))
  background_genes <- row.names(mrna_data_reactive()) %>% unique()
  
  print(paste("The length of background genes is ", length(background_genes)))
  print(head(background_genes))
  
  up_mirnas   <- tryCatch(upregulated_miRNA_reactive() %>% pull(miRNA) %>% unique(), error = function(e) NULL)
  down_mirnas <- tryCatch(downregulated_miRNA_reactive() %>% pull(miRNA) %>% unique(), error = function(e) NULL)


  ### TARGET SET ---
  
  # Optional: restrict to supported interactions only
  if (!is.null(input$corr_select_supported)) {
    if (input$corr_select_supported == "Restrict") {
      supported_neg_cor <- supported_neg_cor_reactive() 
      neg_cor <- neg_cor %>%
        semi_join(supported_neg_cor, by = c("ensembl_gene_id", org_string(), "miRNA", "entrezgene_id"))
    } else if (input$corr_select_supported == "Relax"){
      neg_cor <- neg_cor
    } else if (input$corr_select_supported == "RestrictByPadj"){
      neg_cor <- neg_cor %>% 
        filter(P_adjust_corr < input$corr_padj_threshold)
    }
  }

  # Apply miRNA-direction filtering
  if (!is.null(input$corr_select_up_or_down_miRNA)) {
    if (input$corr_select_up_or_down_miRNA == "Up") {
      neg_cor <- neg_cor %>% 
        filter(miRNA %in% up_mirnas)
      target_genes <- neg_cor %>% 
        pull(entrezgene_id) %>% 
        unique()
    } else if (input$corr_select_up_or_down_miRNA == "Down") {
      neg_cor <- neg_cor %>% dplyr::filter(miRNA %in% down_mirnas)
      target_genes <- neg_cor %>% 
        pull(entrezgene_id) %>% 
        unique()
    } 
  }
  
  # Catch empty target gene list 
  if (length(target_genes) == 0) {
    removeNotification(pathway_notification)
    showNotification("No target genes available for pathway ORA.", type = "error", duration = 8)
    return()
  }
  
  print(paste("The number of target genes is",length(target_genes)))
  
  
  ### BACKGROUND SELECTION ---
  
  if (identical(input$corr_select_background, "Organism_all")) { # CHANGE: new explicit option
    background_genes <- keys(org_db(), keytype = "ENTREZID")
    
  } else if (identical(input$corr_select_background, "Custom_list")) { # CHANGE: new custom input option
    # Expect a textInput where users paste comma-separated ENTREZ IDs, e.g., input$custom_universe_ids
    custom_str <- input$custom_universe_ids %||% ""
    # Parse comma-separated list safely
    background_genes <- if (nzchar(custom_str)) {
      trimws(unlist(strsplit(custom_str, ",")))
    } else character(0)
  } else {
    
    hsmart <- mart()
    
    mapping <- getBM(
      attributes = c('ensembl_gene_id','entrezgene_id'), 
      filters = 'ensembl_gene_id',
      values = gsub("\\..*","",background_genes),
      mart = hsmart
    )
    
    
    write.csv(as.data.frame(mapping),here("dev_files/corr_mapping.csv"))
    background_genes <- mapping %>% 
      mutate(entrezgene_id = as.character(entrezgene_id)) %>%
      pull(entrezgene_id ) %>%
      unique()
  }

  ################ DEV CODE #####################

   #neg_cor <- readr::read_csv(here("dev_files/neg_cor_2025-12-01.csv"))
  #hsmart <- readRDS(here("dev_files/mart.RDS"))
  # supported_neg_cor <- readr::read_csv(here("dev_files/Supported_Negative_Correlations_2025-12-01_.csv"))
  # #
  # neg_cor <- neg_cor %>%
  #    semi_join(supported_neg_cor, by = c("ensembl_gene_id", "hgnc_symbol", "miRNA", "entrezgene_id"))
  # filtered_neg_cor <- neg_cor
  ###############################################

  ### RUN PATHWAY ORA ---

  if (input$corr_analysis_type == "GO") {
    ora_results <- enrichGO(gene = target_genes,
                            universe = background_genes,
                            OrgDb = org_db(),
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = input$corr_pathway_padj_cutoff,
                            qvalueCutoff = input$corr_pathway_padj_cutoff,
                            readable = TRUE)
  } else if (input$corr_analysis_type == "Reactome") {
    ora_results <- enrichPathway(gene = target_genes,
                                 organism = organism_reactive(),
                                 pAdjustMethod = "BH",
                                 universe=background_genes,
                                 pvalueCutoff = input$corr_pathway_padj_cutoff,
                                 qvalueCutoff = input$corr_pathway_padj_cutoff,
                                 readable = TRUE)
  }

  # Handle empty results
  if (is.null(as.data.frame(ora_results)) || nrow(as.data.frame(ora_results)) == 0) {
    output$corr_pathway_ORA_table <- renderDT({
      DT::datatable(
        data.frame(Message = "No enriched pathways found"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
    })
    removeNotification(pathway_notification)
    return()
  }
  
  # Store table form safely
  corr_miRNA_Pathways_reactive(if (!is.null(ora_results)) as.data.frame(ora_results) else data.frame())

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

  if (input$corr_select_up_or_down_miRNA == "Up") {
    mirna_gene_results <- upregulated_miRNA_reactive()
  } else if (input$corr_select_up_or_down_miRNA == "Down") {
    mirna_gene_results <- downregulated_miRNA_reactive()
  }

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
  
  corr_miRNA_mapped_pathways_reactive(mapped_df)
  
  
    
  output$corr_pathway_ORA_table <- renderDT({
    if (!is.null(as.data.frame(ora_results))) {
      removeNotification(pathway_notification)
      datatable(as.data.frame(ora_results), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    } else {
      removeNotification(pathway_notification)
      DT::datatable(
        data.frame(Message = "No enriched pathways found"),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 5)
      )
      return()
    }
  })
  
  # Static table version for paper
  # output$corr_pathway_ORA_table <- renderTable({
  #   if (!is.null(as.data.frame(ora_results))) {
  #     removeNotification(pathway_notification)
  #     head(as.data.frame(ora_results))
  #   } else {
  #     removeNotification(pathway_notification)
  #     DT::datatable(
  #       data.frame(Message = "No enriched pathways found"),
  #       rownames = FALSE,
  #       options = list(scrollX = TRUE, pageLength = 5)
  #     )
  #     return()
  #   }
  # })

})

# Download handlers (outside observeEvents)
output$downloadNegCor <- downloadHandler(
  filename = function() {
    paste("neg_cor_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    req(neg_cor_reactive())
    write.csv(neg_cor_reactive(), file, row.names = FALSE)
  }
)

output$downloadCorrelationPlot <- downloadHandler(
  filename = function() 
    paste0("correlation_plot_", Sys.Date(), ".png"),
  content = function(file) {
    plt_gg <- correlation_plot_reactive()               
    png(file, width = 800, height = 600)
    print(heatmap)
  }
)

output$downloadsupportednegcor <- downloadHandler(
  filename = function() {
    paste("Supported_Negative_Correlations", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(supported_neg_cor_reactive())
    write.csv(supported_neg_cor_reactive(), file, row.names = FALSE, quote = FALSE)
  }
)

# Pathway Visualization - Chord Plot and Network Plot

observeEvent(input$corr_showChordPlot,{
  req(corr_miRNA_mapped_pathways_reactive(),corr_miRNA_Pathways_reactive())
  
  chord_data <- corr_miRNA_mapped_pathways_reactive() %>% 
    dplyr::select(c("miRNA", "description", "enrichment", "coverage")) %>%
    unique()
  
  top_pathways <-  corr_miRNA_Pathways_reactive() %>% 
    as.data.frame() %>% 
    arrange(qvalue) %>% 
    slice_head(n = input$corr_num_pathways) %>%
    pull(Description)
  
  filtered_data <- chord_data %>%
    filter(description %in% top_pathways) %>% 
    filter(coverage >= input$corr_min_coverage)
  
  filtered_data$description <- stringr::str_wrap(filtered_data$description, width=30)
  filtered_data$enrichment <- as.numeric(filtered_data$enrichment)
  
  mat <- as.matrix(xtabs(enrichment ~ miRNA + description, data = filtered_data))
  
  chord_plot <- function() {
    circos.clear()
    chordDiagram(mat, 
                 transparency = 0.5, 
                 grid.col = NULL,
                 annotationTrack = "grid")
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, 
                  CELL_META$ylim[1], 
                  CELL_META$sector.index,
                  facing = "clockwise", 
                  niceFacing = TRUE, 
                  adj = c(0.8),
                  cex = 1)
    }, bg.border = NA)
  }
  
  output$corr_chordPlot <- renderPlot({
    chord_plot()
  })
  
  corr_chord_plot_reactive(chord_plot)
  
  #################### NEW 2026 #####################
  corr_chord_data_reactive(filtered_data)
  
  output$corr_downloadChordData <- downloadHandler(
    filename = function() {
      paste0("corr_chordplot_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- corr_chord_data_reactive()
      req(df)
      readr::write_csv(df, file)
    }
  )
  
  ###################################################
  
})


observeEvent(input$corr_showNetworkPlot, {
  req(corr_miRNA_mapped_pathways_reactive(),corr_miRNA_Pathways_reactive())
  
  # Prepare data
  chord_data <- corr_miRNA_mapped_pathways_reactive() %>%
    dplyr::select(miRNA, description, enrichment, coverage) %>%
    dplyr::distinct()
  
  top_pathways <- corr_miRNA_Pathways_reactive() %>%
    as.data.frame() %>%
    dplyr::arrange(qvalue) %>%
    dplyr::slice_head(n = input$corr_num_pathways) %>%
    dplyr::pull(Description)
  
  filtered_data <- chord_data %>%
    dplyr::filter(description %in% top_pathways) %>%
    dplyr::filter(coverage >= input$corr_min_coverage)
  
  # Build unique edges
  edges <- filtered_data %>%

    dplyr::transmute(
      from = miRNA,
      to   = description,
      title = paste0(
        "<b>miRNA:</b> ", miRNA,
        "<br><b>Pathway:</b> ", description,
        ifelse(!is.na(enrichment), paste0("<br><b>Enrichment:</b> ", enrichment), ""),
        ifelse(!is.na(coverage),   paste0("<br><b>Coverage:</b> ", coverage), "")
      )
    ) %>%
    dplyr::distinct(from, to, .keep_all = TRUE)
  
  # Nodes: miRNA and pathways
  miRNA_nodes <- edges %>%
    dplyr::distinct(from) %>%
    dplyr::transmute(
      id = from, label = from,
      group = "miRNA", color = "#6EC5FF"
    )
  
  pathway_nodes <- edges %>%
    dplyr::distinct(to) %>%
    dplyr::transmute(
      id = to, label = to,
      group = "pathway", color = "#B184F0"
    )
  
  nodes <- dplyr::bind_rows(miRNA_nodes, pathway_nodes)
  
  # Compute degree for tooltips and optional pruning
  deg <- dplyr::bind_rows(
    edges %>% dplyr::count(from, name = "deg") %>% dplyr::rename(id = from),
    edges %>% dplyr::count(to,   name = "deg") %>% dplyr::rename(id = to)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(degree = sum(deg), .groups = "drop")
  
  nodes <- nodes %>%
    dplyr::left_join(deg, by = "id") %>%
    dplyr::mutate(degree = ifelse(is.na(degree), 0L, degree)) %>%
    dplyr::mutate(
      title = ifelse(group == "miRNA",
                     paste0("<b>miRNA:</b> ", label, "<br><b>Degree:</b> ", degree),
                     paste0("<b>Pathway:</b> ", label, "<br><b>Degree:</b> ", degree))
    )
  
  # Optional: prune small nodes to keep it readable
  if (!is.null(input$corr_min_degree) && input$corr_min_degree > 0) {
    keep_ids <- nodes %>% dplyr::filter(degree >= input$corr_min_degree) %>% dplyr::pull(id)
    edges <- edges %>% dplyr::filter(from %in% keep_ids, to %in% keep_ids)
    nodes <- nodes %>% dplyr::filter(id %in% keep_ids)
  }
  
  # Recompute degree after pruning (for accurate tooltips)
  final_deg <- dplyr::bind_rows(
    edges %>% dplyr::count(from, name = "deg") %>% dplyr::rename(id = from),
    edges %>% dplyr::count(to,   name = "deg") %>% dplyr::rename(id = to)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(degree = sum(deg), .groups = "drop")
  
  nodes <- nodes %>%
    dplyr::select(-degree) %>%
    dplyr::left_join(final_deg, by = "id") %>%
    dplyr::mutate(degree = ifelse(is.na(degree), 0L, degree)) %>%
    dplyr::mutate(
      title = ifelse(group == "miRNA",
                     paste0("<b>miRNA:</b> ", label, "<br><b>Degree (shown):</b> ", degree),
                     paste0("<b>Pathway:</b> ", label, "<br><b>Degree (shown):</b> ", degree))
    )
  
  # store exactly what will be plotted
  path_corr_network_displayed(list(nodes = nodes, edges = edges))
  
  # Render interactive network
  output$corr_networkPlot <- visNetwork::renderVisNetwork({
    visNetwork::visNetwork(nodes, edges, height = "700px", width = "100%") %>%
      visNetwork::visNodes(shape = "dot", size = 18, font = list(size = 16)) %>%
      visNetwork::visEdges(smooth = FALSE, color = list(color = "#B3B3B3"), arrows = "none") %>%
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = TRUE
      ) %>%
      visNetwork::visLegend(
        addNodes = list(
          list(label = "miRNA",   shape = "dot", color = "#6EC5FF"),
          list(label = "Pathway", shape = "dot", color = "#B184F0")
        ),
        useGroups = FALSE
      ) %>%
      visNetwork::visInteraction(hover = TRUE, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
      visNetwork::visIgraphLayout(randomSeed = 123)  # FR-like layout
  })
  
  ################################# NEW 2026 ######################################
  output$path_corr_dl_network_nodes <- downloadHandler(
    filename = function() paste0("path_corr_displayed_network_nodes_", Sys.Date(), ".csv"),
    content = function(file) {
      req(path_corr_network_displayed()$nodes)
      readr::write_csv(path_corr_network_displayed()$nodes, file)
    }
  )
  
  output$path_corr_dl_network_edges <- downloadHandler(
    filename = function() paste0("path_corr_displayed_network_edges_", Sys.Date(), ".csv"),
    content = function(file) {
      req(path_corr_network_displayed()$edges)
      readr::write_csv(path_corr_network_displayed()$edges, file)
    }
  )
  ##################################################################################
})
