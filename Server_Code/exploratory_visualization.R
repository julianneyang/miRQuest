# Module 1: Exploratory Visualization

observeEvent(input$plotButton, {
  req(longData(),metadata())
  
  metadataData <- metadata()
  ld <- longData()
  print(names(ld))
  
  ld <- ld %>% 
    group_by(SampleID) %>%
    mutate(Percentage = Count / sum(Count) * 100)
  
  ld <- merge(metadataData, ld, by= "SampleID")
  print(dim(ld))
  
  ld <- ld %>% 
    group_by(SampleID) %>%
    filter(Percentage >= input$percentageFilter) %>%
    mutate(Percentage = Count / sum(Count) * 100)
  print(dim(ld))
  
  cols_large <- unlist(lapply(1:nrow(palettes_d_names), function(i) {
    package_name <- palettes_d_names$package[i]
    palette_name <- palettes_d_names$palette[i]
    num_colors <- palettes_d_names$length[i]
    
    palette_colors <- paletteer_d(paste0(package_name, "::", palette_name), num_colors)
    
    return(palette_colors)
  }))
  
  scrambled_cols <- unique(sample(cols_large))
  
  if (input$meanToggle) {
    req(input$groupBy)
    print(names(ld))
    
    ld_mean <- ld %>%
      group_by(!!sym(input$groupBy), X) %>%
      summarise(MeanPercentage = mean(Percentage)) %>% 
      ungroup() %>% 
      group_by(!!sym(input$groupBy)) %>% 
      mutate(Percentage=MeanPercentage/sum(MeanPercentage)*100)
    print(dim(ld_mean))
    
    output$miRNAPlot <- renderPlot({
      ggplot(ld_mean, aes(x = get(input$groupBy), y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1))
    })
    
    # Store plot for download
    stacked_plot_reactive(function() {
      ggplot(ld_mean, aes(x = get(input$groupBy), y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1))
    })
  } else {
    req(input$groupBy)
    
    output$miRNAPlot <- renderPlot({
      ggplot(ld, aes(x = SampleID, y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
    
    # Store plot for download
    stacked_plot_reactive(function() {
      ggplot(ld, aes(x = SampleID, y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
  }
})

output$facetMessage <- renderUI({
  req(longData(), metadata(), input$groupBy)
  
  # Build the same data as used for plotting (through your existing steps)
  metadataData <- metadata()
  ld <- longData()
  
  # Initial percentage within SampleID
  ld <- ld %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
    dplyr::ungroup()
  
  # Merge metadata
  ld <- merge(metadataData, ld, by = "SampleID")
  
  # Filter by threshold and recompute percentages
  ld <- ld %>%
    dplyr::group_by(SampleID) %>%
    dplyr::filter(Percentage >= input$percentageFilter) %>%
    dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
    dplyr::ungroup()
  
  # Counts per facet group (distinct samples per group in the filtered data)
  facet_counts <- ld %>%
    dplyr::distinct(SampleID, !!rlang::sym(input$groupBy)) %>%
    dplyr::count(!!rlang::sym(input$groupBy), name = "n")
  
  counts_text <- paste(
    paste0(facet_counts[[input$groupBy]], " (n = ", facet_counts$n, ")"),
    collapse = " | "
  )
  
  # Transformation description depends on meanToggle
  if (isTRUE(input$meanToggle)) {
    # Build ld_mean like your plot branch
    ld_mean <- ld %>%
      dplyr::group_by(!!rlang::sym(input$groupBy), X) %>%
      dplyr::summarise(MeanPercentage = mean(Percentage), .groups = "drop") %>%
      dplyr::group_by(!!rlang::sym(input$groupBy)) %>%
      dplyr::mutate(Percentage = MeanPercentage / sum(MeanPercentage) * 100) %>%
      dplyr::ungroup()
    
    transform_text <- HTML(paste0(
      "<b>Data processing (Mean mode):</b><br/>",
      "1) For each SampleID, Percentage = Count / sum(Count) * 100 is computed. <br/>",
      "2) MiRNA meeting the threshold Percentage ≥ ", input$percentageFilter, " are shown. Percentage is then recomputed within each SampleID to scale to 100.<br/>",
      "3) For each level of ", "<code>", input$groupBy, "</code>",
      " and miRNA (X), compute MeanPercentage = mean(Percentage) across samples in that group.<br/>",
      "4) Within each group, Percentage is then recomputed to scale to 100. <br/>"
    ))
  } else {
    transform_text <- HTML(paste0(
      "<b>Data processing (Per-sample mode):</b><br/>",
      "1) For each SampleID, Percentage = Count / sum(Count) * 100 is computed. <br/>",
      "2) MiRNA meeting the threshold Percentage ≥ ", input$percentageFilter, " are shown. Percentage is then recomputed within each SampleID to scale to 100.<br/>",
      "3) Plot per-sample stacked bars (no averaging).<br/>"
    ))
  }
  
  HTML(paste0(
    "<div>",
    "<b>Facet group sample sizes:</b> ", counts_text, "<br/><br/>",
    transform_text,
    "</div>"
  ))
})

observeEvent(input$MDSButton, {
  
  df_wide <- wideData() 
  metadataData <- metadata()
  
  dataDGE<-DGEList(counts=df_wide,genes=rownames(df_wide))
  o <- order(rowSums(dataDGE$counts), decreasing=TRUE)
  dataDGE <- dataDGE[o,]
  
  dataNorm <- calcNormFactors(dataDGE)
  MDSdata <- plotMDS(dataNorm)
  
  MDSxy = data.frame(x=MDSdata$x, y=MDSdata$y)
  colnames(MDSxy) = c(paste(MDSdata$axislabel, '1'), paste(MDSdata$axislabel, '2'))
  mds_xy <- MDSxy
  mds_xy$Dim1 <- mds_xy$`Leading logFC dim 1`
  mds_xy$Dim2 <- mds_xy$`Leading logFC dim 2`
  mds_xy$SampleID <- gsub("X","",colnames(df_wide))
  mds_xy <- merge(mds_xy, metadataData,by="SampleID")
  
  output$MDSPlot <- renderPlot({
    req(input$groupBy)
    grouping_variable <- input$groupBy
    
    plot <- ggplot2::ggplot(mds_xy, aes(x=Dim1, y=Dim2, colour= get(grouping_variable))) + 
      geom_point(aes(fill= get(grouping_variable)), colour="black", pch=21, size=3) +
      scale_fill_viridis_d(name="") +
      ggtitle("Comparison by MDS") +
      xlab("Dim1") +
      ylab("Dim2") +
      cowplot::theme_cowplot(12) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position="top", legend.justification = "center") 
    
    if (input$showSampleID) {
      plot <- plot + geom_text(aes(label=SampleID), vjust=-1, size=3) +
        scale_fill_viridis_d(name="")
    }
    
    print(plot)
    
    # Store in reactive for download handler
    mds_xy_reactive(mds_xy)
    
    # Store plot for download
    mds_plot_reactive(function() {
      ggplot2::ggplot(mds_xy, aes(x=Dim1, y=Dim2, colour= get(input$groupBy))) + 
        geom_point(aes(fill= get(input$groupBy)), colour="black", pch=21, size=3) +
        scale_fill_viridis_d(name="") +
        ggtitle("Comparison by MDS") +
        xlab("Dim1") +
        ylab("Dim2") +
        cowplot::theme_cowplot(16) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="top", legend.justification = "center") +
        (if (input$showSampleID) geom_text(aes(label=SampleID), vjust=-1, size=3) else NULL)
    })
  })
})

# Place this in server, outside observeEvent:
output$MDSMessage <- renderUI({
  req(wideData(), metadata(), input$groupBy)
  
  df_wide <- wideData()
  metadataData <- metadata()
  
  # Build the same MDS data you use for plotting
  dataDGE <- edgeR::DGEList(counts = df_wide, genes = rownames(df_wide))
  o <- order(rowSums(dataDGE$counts), decreasing = TRUE)
  dataDGE <- dataDGE[o, ]
  
  dataNorm <- edgeR::calcNormFactors(dataDGE)
  MDSdata <- limma::plotMDS(dataNorm)
  
  MDSxy <- data.frame(x = MDSdata$x, y = MDSdata$y)
  colnames(MDSxy) <- c(paste(MDSdata$axislabel, "1"), paste(MDSdata$axislabel, "2"))
  mds_xy <- MDSxy
  mds_xy$Dim1 <- mds_xy$`Leading logFC dim 1`
  mds_xy$Dim2 <- mds_xy$`Leading logFC dim 2`
  
  # Derive SampleID from df_wide column names (match plotting code)
  mds_xy$SampleID <- gsub("X", "", colnames(df_wide))
  
  # Merge metadata
  mds_xy <- merge(mds_xy, metadataData, by = "SampleID")
  
  # Counts per group (distinct samples per level of input$groupBy in the MDS data)
  facet_counts <- mds_xy %>%
    dplyr::distinct(SampleID, !!rlang::sym(input$groupBy)) %>%
    dplyr::count(!!rlang::sym(input$groupBy), name = "n")
  
  counts_text <- paste(
    paste0(facet_counts[[input$groupBy]], " (n = ", facet_counts$n, ")"),
    collapse = " | "
  )
  
  # Describe the steps used to generate the MDS
  transform_text <- HTML(paste0(
    "<b>Data processing:</b><br/>",
    "1) Counts are first normalized by Trimmed Mean of M-values (TMM) implemented in edge R to correct for sequencing depth and composition <br/>",
    "2) Samples are then visualized using multidimensional scaling implemented in limma, based on the top 500 genes that most distinguish samples <br/>"
  ))
  
  
  HTML(paste0(
    "<div>",
    "<b>Group sample sizes in MDS:</b> ", counts_text, "<br/><br/>",
    transform_text, "<br/>",
    "</div>"
  ))
})


observeEvent(input$heatmapButton, { 
  countData <- wideData()
  metadataData <- metadata()
  
  metadataData <- metadataData %>% filter(SampleID %in% colnames(countData))
  countData <- countData[,metadataData$SampleID]
  metadataData$SampleID <- factor(metadataData$SampleID, levels =colnames(countData))
  metadataData <- metadataData[order(metadataData$SampleID), ]
  
  designFormula <- as.formula(paste("~", input$groupBy))
  dds <- DESeqDataSetFromMatrix(countData = countData, 
                                colData = metadataData, 
                                design = designFormula)
  
  dds<- DESeq(dds)
  norm_cts <- vst(dds, blind = TRUE, sum( rowMeans( counts(dds, normalized=TRUE)) > 5 )) %>% assay()
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:input$countThreshold]
  
  print(select)
  
  df <- as.data.frame(colData(dds)) %>% 
    dplyr::select(c(input$groupBy))
  
  if (input$rowScale && input$colScale) {
    showNotification("You can only check one scaling option (row or column).", type = "error")
    return()
  }
  
  pheatmap_params <- list(
    mat = norm_cts[select,],
    annotation_col = df,
    show_rownames = input$showRowLabels,
    show_colnames = input$showColLabels
  )
  
  if (input$rowCluster) {
    pheatmap_params$cluster_rows <- TRUE
  } else {
    pheatmap_params$cluster_rows <- FALSE
  }
  
  if (input$colCluster) {
    pheatmap_params$cluster_cols <- TRUE
  } else {
    pheatmap_params$cluster_cols <- FALSE
  }
  
  if (input$rowScale) {
    pheatmap_params$scale <- "row"
  } else if (input$colScale) {
    pheatmap_params$scale <- "column"
  } else {
    pheatmap_params$scale <- "none"
  }
  
  heatmap <- do.call(pheatmap, pheatmap_params)
  
  output$heatmap <- renderPlot({
    heatmap
  })
  
  # Download handler for heatmap
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      fmt <- input$plot_format
      paste0("exploratory_heatmap_", Sys.Date(),  ".", fmt)
    },
    content = function(file) {
      open_graphics_device(
        file,
        format    = input$plot_format,
        width_px  = input$plot_width,
        height_px = input$plot_height,
        res       = 96
      )
      print(heatmap)
      dev.off()
    }
  )
})

# Heatmap Help Text:
output$heatmapMessage <- renderUI({
  req(wideData(), metadata(), input$groupBy, input$countThreshold)
  
  countData <- wideData()
  metadataData <- metadata()
  
  # Align metadata and count matrix to shared samples
  metadataData <- metadataData %>% dplyr::filter(SampleID %in% colnames(countData))
  countData <- countData[, metadataData$SampleID, drop = FALSE]
  metadataData$SampleID <- factor(metadataData$SampleID, levels = colnames(countData))
  metadataData <- metadataData[order(metadataData$SampleID), ]
  
  # Build DESeq2 object and apply VST 
  designFormula <- stats::as.formula(paste("~", input$groupBy))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData = metadataData,
    design = designFormula
  )
  
  dds <- DESeq2::DESeq(dds)
  
  # Variance-stabilized matrix (blind = TRUE, consistent with your code)
  norm_cts <- vst(dds, blind = TRUE, sum( rowMeans( counts(dds, normalized=TRUE)) > 5 )) %>% SummarizedExperiment::assay()
  
  # Select top features by mean normalized counts
  select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)),
                  decreasing = TRUE)[seq_len(min(input$countThreshold, nrow(countData)))]
  
  # Group composition: counts per level of input$groupBy among included samples
  group_counts <- metadataData %>%
    dplyr::count(!!rlang::sym(input$groupBy), name = "n")
  
  group_text <- paste(
    paste0(group_counts[[input$groupBy]], " (n = ", group_counts$n, ")"),
    collapse = " | "
  )
  
  # Determine clustering and scaling options
  row_cluster <- if (isTRUE(input$rowCluster)) "Yes" else "No"
  col_cluster <- if (isTRUE(input$colCluster)) "Yes" else "No"
  scale_mode <- if (isTRUE(input$rowScale) && isTRUE(input$colScale)) {
    "Both (invalid) — only one can be selected"
  } else if (isTRUE(input$rowScale)) {
    "Row (each miRNA standardized across samples)"
  } else if (isTRUE(input$colScale)) {
    "Column (each sample standardized across miRNAs)"
  } else {
    "None (values in VST scale)"
  }
  
  # Compose transformation/analysis description in plain language
  transform_text <- HTML(paste0(
    "<b>Data processing:</b><br/>",
    "1) Counts are first adjusted by DESEQ2 size factors to correct for sequencing depth and sample composition.<br/>",
    "2) A variance stabilizing transformation (VST implemented in DESeq2) is then applied to make variability more uniform across expression levels, ",
    "VST prevents very highly expressed miRNAs from dominating the color scale, so patterns are easier to see and compare.<br/>",
    "3) Select the top ", length(select), " miRNA by mean normalized counts to focus the heatmap on the most expressed features.<br/>",
    "4) Optionally cluster rows/columns to reveal patterns, and apply scaling (row or column) if selected.<br/>"
  ))
  
  notes <- HTML(paste0(
    "<b>Notes:</b><br/>",
    "- VST values are on a transformed scale (not raw counts); color differences reflect relative expression after normalization.<br/>",
    "- Row scaling standardizes each gene (mean 0, unit variance), emphasizing within-gene contrasts across samples.<br/>",
    "- Column scaling standardizes each sample across genes, emphasizing each sample’s profile.<br/>"
  ))
  
  HTML(paste0(
    "<div>",
    "<b>Selected genes:</b> ", length(select), "<br/>",
    "<b>Group composition:</b> ", group_text, "<br/><br/>",
    transform_text, "<br/>",
    "<b>Clustering:</b> rows = ", row_cluster, ", cols = ", col_cluster, "<br/>",
    "<b>Scaling:</b> ", scale_mode, "<br/><br/>",
    notes,
    "</div>"
  ))
})

# Download handler for MDS data (outside observeEvent)
output$downloadData <- downloadHandler(
  filename = function() {
    paste("mds_xy_", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    req(mds_xy_reactive())
    write.csv(mds_xy_reactive(), file, row.names = FALSE)
  }
)