# All download handlers must be defined at the top level, not inside observeEvents


open_graphics_device <- function(file, format,
                                 width_px, height_px,
                                 res = 96) {
  format <- tolower(format)
  
  if (format == "png") {
    png(file,
        width  = width_px,
        height = height_px,
        res    = res)
    
  } else if (format == "tiff") {
    # Requires TIFF capability on the server: capabilities("tiff") should be TRUE
    tiff(file,
         width       = width_px,
         height      = height_px,
         res         = res,
         compression = "lzw")
    
  } else if (format == "pdf") {
    # pdf uses inches: convert pixels → inches with 'res'
    pdf(file,
        width  = width_px  / res,
        height = height_px / res)
    
  } else if (format == "eps") {
    # EPS via postscript(); also uses inches
    postscript(file,
               width      = width_px  / res,
               height     = height_px / res,
               horizontal = FALSE,
               onefile    = FALSE,
               paper      = "special")
    
  } else if (format == "svg") {
    # SVG (vector); also uses inches
    svg(file,
        width  = width_px  / res,
        height = height_px / res)
    
  } else {
    stop("Unsupported format: ", format)
  }
}



# Module 1: Exploratory Visualization Downloads
output$downloadStackedColumn <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Stacked_Column_Chart_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- stacked_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$downloadPlot <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("MDS_Plot_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- mds_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# Module 3: Differential miRNA Expression Downloads  
output$downloadVolcano <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Volcano_Plot_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- volcano_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$downloadSignificantHeatmap <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Significant_Heatmap_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- significant_heatmap_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# Module 5: Single miRNA Target Prediction Downloads
output$downloadDotplot <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Pathway_Dotplot_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- dotplot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$downloadBarplot <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Pathway_Barplot_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- barplot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# Module 4: Aggregated Pathway Downloads
output$downloadChord <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Chord_Plot_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- chord_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# output$downloadNetwork <- downloadHandler(
#   filename = function() {
#     fmt <- input$plot_format
#     paste0("Network_Plot_", Sys.Date(), ".", fmt)
#   },
#   content = function(file) {
#     open_graphics_device(
#       file,
#       format    = input$plot_format,
#       width_px  = input$plot_width,
#       height_px = input$plot_height,
#       res       = 300
#     )
#     plot_func <- network_plot_reactive()
#     if (is.function(plot_func)) {
#       print(plot_func())
#     } else {
#       print(plot_func)
#     }
#     dev.off()
#   }
# )

# Module 7: Correlation Analysis Downloads
output$downloadCorrelationPlot <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Correlation_Plot_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- correlation_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$corr_downloadChord <- downloadHandler(
  filename = function() {
    fmt <- input$plot_format
    paste0("Chord_Plot_", Sys.Date(), ".", fmt)
  },
  content = function(file) {
    open_graphics_device(
      file,
      format    = input$plot_format,
      width_px  = input$plot_width,
      height_px = input$plot_height,
      res       = 96
    )
    plot_func <- chord_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)
