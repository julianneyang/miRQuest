# Define UI for application

ui <- fluidPage(
  useShinyjs(),
  titlePanel("miRQuest: Interactive Analysis of MicroRNA Sequencing Data"),
  theme=bs_theme(version=5, bootswatch = "cosmo"),
  
  sidebarLayout(
    sidebarPanel(
      actionButton("useDefaultData", "Use Default Data", class = "btn-primary", style = "width: 100%; margin-bottom: 20px;"),
      tags$hr(),
      h6("Upload raw miRNA counts x sample table"),
      fileInput("countTable", "Upload miRNA Count Table", 
                accept = c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
      textOutput("mirna_status"),
      tags$hr(),
      h6("Upload sample metadata; assumes first column is named 'SampleID'"),
      fileInput("metadataFile", "Upload Metadata File", 
                accept = c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
      textOutput("metadata_status"),
      tags$hr(),
      h6("Optionally, add RNA-Seq data for the same samples"),
      fileInput("mrna_countTable", "Upload RNA-Seq Count Table", 
                accept = c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
      textOutput("mrna_status"),
      tags$hr(),
      h6("Select species"),
      checkboxInput("hsa", "Homo sapiens", FALSE),
      checkboxInput("mmu", "Mus musculus", FALSE),
      tags$hr(),
      
      # make user upload MTI file
      checkboxInput("validated", "Retrieve validated gene targets", FALSE),
      conditionalPanel(
        condition = "input.validated && input.hsa && !input.mmu",
        fileInput("file_hsa", "Upload hsa_MTI.csv", accept = c(".csv"))
      ),
      conditionalPanel(
        condition = "input.validated && input.mmu && !input.hsa",
        fileInput("file_mmu", "Upload mmu_MTI.csv", accept = c(".csv"))
      ),
      
      fluidRow(
        column(8, actionButton("submit", "Submit")),  # Submit button
        column(7, textOutput("error_message"))  # Error message next to the button
      ),
      tags$hr(),
      # Download options UI ----
      selectInput(
        "plot_format",
        "Download format (can change on per-plot basis):",
        choices = c(
          "PNG"       = "png",
          "TIFF" = "tiff",
          "PDF"       = "pdf",
          "EPS"       = "eps",
          "SVG"       = "svg"
        ),
        selected = "png"
      ),
      numericInput(
        "plot_width",
        "Plot width (pixels):",
        value = 800,
        min   = 400,
        max   = 4000,
        step  = 100
      ),
      numericInput(
        "plot_height",
        "Plot height (pixels):",
        value = 600,
        min   = 300,
        max   = 4000,
        step  = 100
      ),
      
      width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Quickstart",
            fluidPage(
            fluidRow(
              column(
                width = 12,
                div(
                  style = "margin-top: 10px; margin-bottom: 10px; color: #555;",
                  strong("Maintained as of: "),
                  span("2026-03-10")  
                  )
                )
            ),
            fluidRow(
              column(
                width = 12,
                h2("Welcome to MiRQuest Quickstart"),
                p("This section provides a brief overview to get you up and running quickly. 
           Use the diagram below to understand the suggested workflow, then explore 
           the rest of the app using the tabs above.")
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 12,
                h3("GitHub Repository"),
                p("Find the source code, documentation, example files, and issue tracker here:"),
                tags$a(
                  href = "https://github.com/MSDLLCpapers/miRQuest",  
                  target = "_blank",
                  "Open GitHub Repository"
                ),
                br(),
                p("If you encounter issues or have feature requests, please open an issue on GitHub.")
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 12,
                h3("Getting Started"),
                tags$ol(
                  tags$li("To begin, please go to the side panel and load your miRNA count data and metadata files. By default, demo data is loaded."),
                  tags$li("Please select your species and hit the Submit button. 
                          If you would like to retrieve validated gene targets from miRTarBase, please download and upload the species-specific *_MTI.csv file when prompted, from this link: "),
                  tags$a(
                    href = "https://awi.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/php/download.php",  
                    target = "_blank",
                    "Download miRTarBase file"
                  ),
                  tags$li("You can then either proceed to Exploratory Visualization or Differential miRNA Analysis."),
                  tags$li("For more help, please follow the arrows in the workflow diagram below.")
                )
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 12,
                h3("Suggested Use Diagram"),
                p("The diagram below illustrates the recommended workflow for using this application:"),
                tags$img(
                  src = "App_Schematic.png",  
                  alt = "Suggested Use Diagram",
                  style = "max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 6px;"
                )
              )
            ),
            hr(), 
            fluidRow(
              column(
                width = 12,
                h3("Walkthrough Tutorial Using the Demo Data"),
                p("Please see the following tutorial for more help"),
                tags$iframe(
                  src = "MiRQuest_tutorial.pdf",
                  style = "width:100%; height:600px; border:none;"
                )
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 12,
                h3("Manuscript"),
                p("To see an application of miRQuest and interpretation of results, please see the following paper:"),
                p("Yang JC, Sauter J, Adam G, Carr R. MiRQuest: An Interactive App for MiRNA Sequencing Analysis. in revision at BMC Bioinformatics."),
                p("If you use miRQuest in your work, please cite both the paper and the source code:"),
                p("Zenodo placeholder")
              )
            )
          )
        ),
        
        tabPanel("Exploratory Visualization",
                 selectInput("groupBy", "Select Metadata Column to Group Plots", choices = NULL),
                 helpText("You must select a column that encodes groups containing > 1 sample for the heatmap and volcano plots, e.g. not SampleID."),
                 tags$hr(),
                 fluidRow(
                   column(4,
                    checkboxInput("meanToggle", "Calculate Mean Relative Abundance", value = FALSE),  
                   ),
                   column(8, 
                          helpText("For the stacked column charts, if you are constrained by the space available, please try subsetting the data or shorten your SampleIDs.")
                          ),
                 ),
                 actionButton("plotButton", "Show miRnome as Stacked Column Chart"),
                 sliderInput("percentageFilter", 
                             "% miRNAs that comprise the sample", 
                             min = 0, 
                             max = 5, 
                             value = 1, 
                             step = 0.2),
    
                 fluidRow(
                   column(12,
                          uiOutput("facetMessage"),
                          plotOutput("miRNAPlot")
                   )
                 ),
                 downloadButton("downloadStackedColumn", "Download Plot"),
                 tags$hr(),
            
                 actionButton("MDSButton", "Show miRnome as a Single Point in 2D Space"),
                 checkboxInput("showSampleID", "Show Sample ID on MDS Plot", value = FALSE),
               
                 fluidRow(
                   column(12,
                          uiOutput("MDSMessage"),
                          plotOutput("MDSPlot")
                   )
                 ),
                 downloadButton("downloadData", "Download Dataframe"),
                 downloadButton("downloadPlot", "Download Plot"),
                 tags$hr(),

                 actionButton("heatmapButton", "Show miRnome as Heatmap"),
                 sliderInput("countThreshold", "Top n number of miRNA with the most counts", value = 20, min = 5, max=50, step = 5),
                 # Organizing checkboxes into three rows with two checkboxes per row
                 fluidRow(
                   column(6, checkboxInput("colCluster", "Cluster by Sample", value = TRUE)),
                   column(6, checkboxInput("rowCluster", "Cluster by miRNA", value = FALSE))
                 ),
                 fluidRow(
                   column(6, checkboxInput("rowScale", "Scale by miRNA", value = TRUE)),
                   column(6, checkboxInput("colScale", "Scale by Sample", value = FALSE))  
                 ),
                 fluidRow(
                   column(6, checkboxInput("showRowLabels", "Show Row Labels", value = FALSE)),  
                   column(6, checkboxInput("showColLabels", "Show Column Labels", value = FALSE))   
                 ),
                 fluidRow(
                   column(12,
                          uiOutput("heatmapMessage"),
                          plotOutput("heatmap")
                   )
                 ),
                 
                 downloadButton("downloadHeatmap", "Download Heatmap"),
                 tags$hr()
                 
                 
        ),
        tabPanel("Data Subsetting",
                 selectInput("metadataColumn", "Select Metadata Column", choices = NULL),
                 uiOutput("filterValueInput"),  # Dynamic UI for filter values
                 actionButton("subsetDataButton", "Generate Subset"),
                 downloadButton("downloadSubset", "Download Subset as CSV"),  # Download button
                 tableOutput("subsetDataTable")  # Table to display the subsetted data
                 
        ),
        tabPanel("Differential miRNA Expression Analysis",
                 #selectInput("designColumn", "Select Design Column", choices = NULL),
                 
                 ################################ NEW CHANGES ################################
                 # ADDED: select multiple covariates for the design
                 selectizeInput(
                   "designColumns",
                   "Design covariates (main effects):",
                   choices = NULL,  # populate server-side after metadata is available
                   multiple = TRUE
                 ),  # ADDED
                 
                 # ADDED: choose the variable of interest (one of the selected covariates)
                 selectInput(
                   "voi",
                   "Variable of interest (VOI):",
                   choices = NULL
                 ) , # ADDED
                 
                 # ADDED: choose reference level for the variable of interest (shown only if VOI is a factor)
                 uiOutput("voiRefLevelUI"),  # ADDED
                 textOutput("voi_sample_sizes"),       # CHANGED: VOI sample sizes text
                 uiOutput("covariate_sample_sizes") ,  # CHANGED: list of covariate sample sizes
                 #############################################################################
                 
                 actionButton("runDESeq", "Run DESeq2 Analysis"),
                 numericInput("DEM_padj_filter", "Adjusted P-value Threshold", 
                              value = 0.05, min = 0, max = 1, step = 0.01),
                 #actionButton("showDESeqTable",  "Preview Differentially Expressed miRNA"), 
                 textOutput("dem_count"), 
                 downloadButton("download_significant_miRNA", "Download Significantly Differentially Expressed miRNA"),
                 DTOutput("resSignificantTable"), 
                 tags$hr(),
                 selectInput("selectedMiRNA", "Select Differentially Expressed miRNA", choices = NULL),
                 selectInput("groupBy_Boxplot", "Select Metadata Column For Plot Grouping", choices = NULL),
                 # In your UI, near the boxplot output:
                 fluidRow(
                   column(12,
                          uiOutput("boxplotMessage"),
                          plotOutput("boxplot")
                   )
                 ),
                 downloadButton("downloadBoxplot", "Download Plot"),
                 downloadButton("download_plotCounts_selected", "Download Data"),
                 tags$hr(),
                 h6("The following Log2FoldChange filter affects both the volcano and heatmap plots. Click the button to redraw the plots."),
                 numericInput("DEM_log2FoldChange_filter", "Absolute Log2 Fold Change Threshold", 
                              min = 0, 
                              value = 1, 
                              step = 0.5),
                 # Multi-select: choices populated reactively (optional) or left empty for free typing

                 selectizeInput(
                   inputId = "volcano_select_labels",
                   label = "Label only these miRNAs (leave empty to label all)",
                   choices = NULL,         # will be populated from server
                   multiple = TRUE,
                   options = list(create = TRUE, placeholder = "Type to search or add miRNAs")
                 ),
                 
                 # Optional message area
                 textOutput("volcano_label_warning"),
                 sliderInput("volcano_plot_label_size", "Label Size",
                             min = 1, 
                             max = 6, 
                             value=3,
                             step=1),
                 # future_expansion_code: limma::voom analysis button
                 # actionButton("runVoom", "Run limma::voom Analysis"),  # New button for limma::voom
                 actionButton("Volcano_Plot_Button", "Show DE miRNA as Volcano Plot"),
                 plotOutput("volcanoPlot"),
                 downloadButton("downloadVolcano", "Download Volcano Plot"),
                 # future_expansion_code: limma volcano plot download
                 # downloadButton("downloadVolcanoLimma", "Download limma Volcano Plot"),
                 tags$hr(),
                 # Controls 
                 selectizeInput(
                   inputId = "heatmap_select_rows",
                   label = "Show only these miRNAs in heatmap (leave empty to use significant filter)",
                   choices = NULL,          # populated from server
                   multiple = TRUE,
                   options = list(create = TRUE, placeholder = "Type to search or add miRNAs")
                 ),
                 actionButton("significant_heatmap_button", "Show DE miRNA as Heatmap"),
                 checkboxInput("significant_heatmap_showRowLabels", "Show Row Labels", value = FALSE),
                 checkboxInput("significant_heatmap_showColLabels", "Show Column Labels", value = FALSE),
                 tags$div(
                   tags$strong("Data processing: "),
                   tags$br(),
                       "1) Counts are first adjusted by DESEQ2 size factors to correct for sequencing depth and sample composition.",
                   tags$br(),
                       "2) A variance stabilizing transformation (VST implemented in DESeq2) is then applied to make variability more uniform across expression levels.",
                   tags$br(),
                       "3) This heatmap shows only the significantly differentially expressed miRNA, scaled by row. Clustering is automatically done for both samples and miRNAs."
                 ),
                 plotOutput("significant_miRNA_heatmap"),
                 downloadButton("downloadSignificantHeatmap", "Download Plot"),
                 downloadButton("download_heatmap_matrix", "Download Data"),
                 
                 tags$hr(),
                 # future_expansion_code: DE miRNA data table for debugging/correlation analysis
                 # h5("DE miRNA Data for Correlation Analysis:"),
                 # tableOutput("mirna_DE_table"),
                 # downloadButton("download_mirna_DE_table", "Download DE miRNA Data"),
                 # tags$hr()
                 # Plot Outputs
                 
                 
                 
        ),
        
        tabPanel("Single DE miRNA Target Retrieval and Pathway Enrichment",
               
                 
                 
                 ###################################### NEW CHANGES ################################
                 fluidRow(
                   column(4, selectInput("DEMiRNA", "Select Differentially Expressed miRNA", choices = NULL)),
                   column(4,
                     radioButtons(
                     inputId = "result_scope_all",
                     label = "Show Predicted and/or Experimentally Validated Gene Targets",
                     choices = c(
                       "Both (predicted ∩ validated)" = "both",
                       "Predicted only" = "predicted",
                       "All (predicted ∪ validated)" = "all",
                       "Validated only" = "validated"
                     ),
                     selected = "both",
                     inline = FALSE
                    )
                 ),
                   column(4,
                     selectInput(
                       inputId = "sources",
                       label = "If predicted gene targets are desired, select source databases",
                       choices =  c("pictar", "diana", "targetscan", "miranda", "mirdb"),
                       selected =  c("pictar", "diana", "targetscan", "miranda", "mirdb"),   # default: all selected
                       multiple = TRUE
                     ),
                     uiOutput("single_min_src_ui")
                   )
                 ),
                 ############################################################################
                 
                 actionButton("runAnalysis", "Retrieve Gene Targets for a Single miRNA"),
                 

                 downloadButton("downloadgeneResults", "Download Predicted Genes"),
                 DTOutput("genetargets"), 
                 tags$hr(),
                 fluidRow(
                   column(6,
                      selectInput("single_miRNA_analysis_type", "Select Pathway Database",
                              choices = c("GO" = "GO", "Reactome" = "Reactome"),
                              selected = "Reactome"),  
                   ),
                   column(6, 
                # Background universe selector
                     selectInput("select_background_single", "Background universe for ORA",
                       choices = c(
                         "Organism_all",          # all ENTREZ IDs in OrgDb
                         "Custom_list",           # user-entered comma-separated ENTREZ IDs
                         "Validated_universe"    # union of miRTarBase targets (only when scope == validated_agg)
                       ),
                       selected = "Organism_all"
                     )
                   )
                ),
                 tags$div(
                   tags$p(strong("Background options and their biological questions:")),
                   tags$ul(
                     tags$li(
                       tags$strong("Organism_all: "),
                       "Are the target genes enriched in specific biology relative to all annotated genes in the organism?"
                     ),
                     tags$li(
                       tags$strong("Validated_universe: "),
                       "Are the target genes enriched compared to all miRTarBase-validated targets (same organism)? ",
                       em("Only applied when retrieving target genes from miRTarBase (validated).")
                     ),
                     tags$li(
                       tags$strong("Custom_list: "),
                       "Are the target genes enriched relative to a user-defined gene context (paste ENTREZ IDs below)?"
                     )
                   )
                 ),
                 
                 fluidRow(
                   column(12,
                     # Show custom input only when Custom_list is selected
                     conditionalPanel(
                       condition = "input.select_background == 'Custom_list'",
                       textInput(
                         "single_custom_universe_ids",
                         "Custom universe (comma-separated ENTREZ IDs)",
                         placeholder = "e.g., 7157, 1956, 7422"
                       ),
                       helpText("Paste ENTREZ gene IDs separated by commas. The analysis will test enrichment against this custom background.")
                     )
                   )
                 ),
                 
                 helpText("Note: 'Validated_universe' uses all miRTarBase targets and is only applied when targets retrieved are 'Validated' only. For other scopes, the app will fall back to 'Organism' unless you choose another option."),
                 numericInput("ORA_single_padj_Cutoff", "FDR-adjusted P-value < ", value = 0.05, min = 0, max = 1, step = 0.01),
                 actionButton("runORA_single_miRNA", "Perform Pathway Overrepresentation Analysis"), 
                 downloadButton("downloadFinalResults", "Download Enriched Pathways"),
                 DTOutput("resultsTable"),
                 tags$hr(),
                 actionButton("plot_pathways_single_miRNA", "Plot Enriched Pathways"), 
              layout_columns(
                card(card_header("Dotplot"),
                 plotOutput("ora_dotplot"),
                 downloadButton("downloadDotplot", "Download Plot")
                ),
                card(card_header("Barplot"),
                 plotOutput("ora_barplot"),
                 downloadButton("downloadBarplot","Download Plot")
                )
              )
        ),
        
        tabPanel("All DE miRNA Target Retrieval and Pathway Enrichment",
                 fluidRow(
                   column(4,
                     radioButtons(
                       inputId = "result_scope_all_agg",
                       label = "Show Predicted and/or Experimentally Validated Gene Targets",
                       choices = c(
                         "Both (predicted ∩ validated)" = "both_agg",
                         "Predicted only" = "predicted_agg",
                         "Validated only" = "validated_agg",
                         "All (predicted ∪ validated)" = "all_agg"
                       ),
                       selected = "both_agg",
                       inline = FALSE
                     )
                     
                   ),
                  column(4,helpText("If you have RNA-seq data, consider choosing a relaxed target retrieval method such as 'All' so you can use the 'miRNA-Gene Correlation Analysis' module to identify anti-correlated pairs with database support. 
                                    If you only have miRNA-seq data, consider choosing a more stringent method such as 'Both'.")
                  ),
                  column(4,
                     selectInput(
                       inputId = "sources_agg",
                       label = "If predicted gene targets are desired, select source databases",
                       choices =  c("pictar", "diana", "targetscan", "miranda", "mirdb"),
                       selected =  c("pictar", "diana", "targetscan", "miranda", "mirdb"),   # default: all selected
                       multiple = TRUE
                     ),
                     uiOutput("min_src_ui")
                     )
                 ),
                 actionButton("predict_all_genes", "Retrieve Targeted Genes for all DE miRNA"),  
                 
                 
                 downloadButton("all_downloadgeneResults", "Download Predicted Genes"),
                 DTOutput("all_genetargets"), 
                 tags$hr(),
                 
                 ########################### NEW CHANGES ##################################
                 fluidRow(
                  column(3,  h6("This will show you the top n nodes (miRNA or gene) with the most connections."),
                         numericInput("max_nodes","Max nodes:", value = 200, min = 10, step = 10)
                         
                         ),
                  column(3,  h6("This will enable you to show only nodes that have at least n number connections"),
                         numericInput("min_degree",
                                         "Minimum node degree:", value = 1, min = 0)
                         ),
                 
                 
                  column(3, radioButtons(
                   "miRNA_dir",
                   "miRNA direction to show",
                   choices = c("Both" = "both", "Upregulated only" = "up", "Downregulated only" = "down"),
                   selected = "both",
                   inline = TRUE
                  )),
                  # column(3, selectInput("miRNA_filter", "Filter miRNA", choices = NULL, multiple = TRUE)) # future expansion: filter by miRNA or gene
                 
                 ),
                 actionButton("showGeneNetwork", "Show miRNA–Gene Network"),
                 visNetworkOutput("miRNAGeneNetwork", height = "700px"),
                 downloadButton("dl_network_nodes", "Download nodes"),
                 downloadButton("dl_network_edges", "Download edges"),
                 ################################################################################
                 
                 tags$hr(),
                 
                 fluidRow(
                   column(4,
                     selectInput("select_up_or_down_miRNA", "Select Up- or Down-regulated miRNA",
                                 choices = c("Up" = "Up", "Down" = "Down"),
                                 selected = "Up")
                   ),
                   column(4,
                     # Background universe selector
                     selectInput(
                       "select_background",
                       "Background universe for ORA",
                       choices = c(
                         "Organism_all",          # all ENTREZ IDs in OrgDb
                         "Custom_list",           # user-entered comma-separated ENTREZ IDs
                         "Validated_universe"    # union of miRTarBase targets (only when scope == validated_agg)
                       ),
                       selected = "Organism_all"
                     )
                   ),
                    column(4,
                     selectInput("analysis_type", "Select Pathway Database",
                                 choices = c("GO" = "GO", "Reactome" = "Reactome"),
                                 selected = "Reactome")
                    )
                 ),
                 # CHANGE: Add helper text that explains the biological question each universe answers
                 tags$div(
                   tags$p(strong("Background options and their biological questions:")),
                   tags$ul(
                     tags$li(
                       tags$strong("Organism_all: "),
                       "Are the target genes enriched in specific biology relative to all annotated genes in the organism?"
                     ),
                     tags$li(
                       tags$strong("Validated_universe: "),
                       "Are the target genes enriched compared to all miRTarBase-validated targets (same organism)? ",
                       em("Only applied when retrieving target genes from miRTarBase (validated).")
                     ),
                     tags$li(
                       tags$strong("Custom_list: "),
                       "Are the target genes enriched relative to a user-defined gene context (paste ENTREZ IDs below)?"
                     )
                   )
                 ),
                 
                 # Show custom input only when Custom_list is selected
                 conditionalPanel(
                   condition = "input.select_background == 'Custom_list'",
                   textInput(
                     "custom_universe_ids",
                     "Custom universe (comma-separated ENTREZ IDs)",
                     placeholder = "e.g., 7157, 1956, 7422"
                   ),
                   helpText("Paste ENTREZ gene IDs separated by commas. The analysis will test enrichment against this custom background.")
                 ),
                 
                 # Optional guidance about Validated_universe applicability
                 helpText("Note: 'Validated_universe' uses all miRTarBase targets and is only applied when retrieved targets are 'Validated' only. For other scopes, the app will fall back to 'Organism' unless you choose another option."),
        
                 numericInput("ORA_all_padj_Cutoff", "FDR-adjusted P-value < ", value = 0.05, min = 0, max = 1, step = 0.01),
                 
                 actionButton("runORA", "Perform Pathway Overrepresentation Analysis"),
                 downloadButton("all_downloadFinalResults", "Download Enriched Pathways"),
                 DTOutput("all_resultsTable"),
                 tags$hr(),
                 
                 
                 # New inputs for user-defined parameters
                 h6("The following max. pathway and min. coverage controls affect both the chord and network plots."),
                 h6("The number of pathways shown will be affected by coverage selection."),
                 numericInput("num_pathways", "Maximum Number of Pathways to Show:", value = 10, min = 2),
                 h6("Coverage indicates the minimum fraction of genes in a pathway an miRNA must target to be shown."),
                 numericInput("min_coverage", "Minimum Coverage (0 to 1):", value = 0.1, min = 0, max = 1, step = 0.01),

                 tags$hr(),
                 actionButton("showChordPlot", "Show Chord Plot"),
                 tags$div(
                    style = "overflow: visible; width: 100%; height: 600px;",
                    plotOutput("chordPlot", width = "100%", height = "600px")
                 ),
                 downloadButton("downloadChord", "Download Plot"),
                 downloadButton("downloadChordData", "Download Data"),
                 numericInput("path_min_degree", "Minimum node degree: This will enable you to show only nodes that have at least n number connections", value = 1, min = 0),
                 
                 actionButton("showNetworkPlot", "Show Network Plot"),
                 visNetworkOutput("networkPlot", height = "700px"),
                 downloadButton("path_dl_network_nodes", "Download nodes"),
                 downloadButton("path_dl_network_edges", "Download edges"),

        ),
        tabPanel("Differential Gene Expression Analysis",
                 h6("Uses the same model design as is selected for miRNA"),
                 actionButton("mrna_runDESeq", "Run DESeq2 Analysis"),  
                 numericInput("padj_filter", "Adjusted P-value Threshold", 
                              value = 0.05, min = 0, max = 1, step = 0.01),
                 textOutput("deg_count"),  # Output for displaying the number of DEGs
                 downloadButton("download_DEG_Table", "Download Differentially Expressed Genes"),
                 tableOutput("mrna_resSignificantTable"), 
                 tags$hr()
        ),
        
        tabPanel("miRNA - Gene correlation analysis",
                 h6("Pairwise miRNA-gene Correlations"),
                 helpText("Computes pairwise correlations for all selected DE miRNAs and DE genes. Runtime increases with the number of features (miRNAs × genes). CAUTION: Please heed the numbers of miRNA and genes below (for #### x ### features, this can take 30 minutes)"),
                 
                 fluidRow(
                   column(4, numericInput("filter_DEG_correlation", "mRNA Log2 Fold Change Threshold", 
                                          value = 1, min = 0, step = 0.1)),
                   column(4, numericInput("filter_DEM_correlation", "miRNA Log2 Fold Change Threshold", 
                                          value = 1, min = 0, step = 0.1)),
                   column(4, numericInput("corr_padj_cutoff", "FDR-adjusted P-value < ", value = 0.05, min = 0, max = 1, step = 0.01))
                 ),
                 textOutput("DEG_DEM_filter_message"),
                 actionButton("runCorrelationAnalysis", "Run Correlation Analysis"), 
                 textOutput("num_neg_correlations"),
                 downloadButton("downloadNegCor", "Download Significant Negative Correlations"),
                 helpText("Please note that we will pass all miRNA-gene pairs with Spearman rho < 0 regardless of significance threshold to subsequent modules below, and we have filtered only the downloadable results to the specified threshold for your convenience."),
                 DTOutput("NegativeCorrelationTable"),
                 #tableOutput("NegativeCorrelationTable"),
                 tags$hr(),
                 ####################################################### NEW 2026 ###############################################
                 h6("Optional global permutation test of negative correlations"),
                 helpText("Tests whether the overall number of significantly negative miRNA–mRNA correlations reported above is higher than expected by chance. CAUTION: Runtime will be (number of permutations) x (length of time for the pairwise correlations above)!!! "),
                 
                 numericInput("n_perm_correlation", "Number of permutations", value = 100, min = 10, step = 50),
                 
                 actionButton("runGlobalPermutation", "Run global permutation test"),
                 
                 verbatimTextOutput("GlobalPermutationSummary"),
                 tags$hr(),
                 #################################################################################################################
                 actionButton("show_database_support", "See Correlations with Database Support"), 
                 helpText("This will identify the intersection of negatively correlated miRNA-gene pairs (whether significant or not) with targets retrieved for all DE miRNA."),
                 textOutput("degMessage"),
                 downloadButton("downloadsupportednegcor", "Download Supported Negative Correlations"),
                 DTOutput("SupportedNegativeCorrelationTable"),
                 #tableOutput("SupportedNegativeCorrelationTable"),
                 tags$hr(),
                 h6("The following controls affect both the correlation and network plots. To redraw the plots, click the button again"),
                 checkboxInput("select_supported_only", 
                               label = "Show Supported miRNA-Gene Pairs Only", 
                               value = FALSE),
                 fluidRow(
                   column(3, numericInput("correlationCutoff", "Correlation Cutoff:", 
                                          value = -0.5, min = -1, max = 0, step = 0.1)),  
                   column(3, numericInput("correlation_padj_Cutoff", "Correlation P-adj < ", 
                                          value = 0.05, min = 0, max = 1, step = 0.01)),  
                   column(3, numericInput("logratio_miRNA_Cutoff", "Absolute DE miRNA log2FoldChange Cutoff:", 
                                          value = 1, min = 0, step = 1)),  
                   column(3, numericInput("logratio_gene_Cutoff", "Absolute DE Gene log2FoldChange Cutoff", 
                                          value = 1, min = 0, step = 1))  
                 ),
                
                 tags$hr(),
                 actionButton("CorrelationPlotButton", "Show Correlation Plot"),
                 fluidRow(
                   column(12,
                          uiOutput("correlationMessage"),
                          plotlyOutput("correlationPlot")  
                   )
                 ),
                 downloadButton("downloadCorrelationPlot", "Download Plot"),
                 fluidRow(
                   column(4, h6("This will show you the top n nodes (miRNA or gene) with the most connections."),
                          numericInput("max_nodes_experimental", "Max nodes:", value = 200, min = 10, step = 10)),
                   column(4, h6("This will enable you to show only nodes that have at least n number connections"),
                          numericInput("min_degree_experimental", "Minimum node degree:", value = 10, min = 0)),
                   
                   
                   column(4, radioButtons(
                     "corr_miRNA_dir",
                     "miRNA direction to show",
                     choices = c("Both" = "both", "Upregulated only" = "up", "Downregulated only" = "down"),
                     selected = "both",
                     inline = TRUE
                   )),
                 ),
                 actionButton("showGeneNetwork_Experimental", "Show miRNA–Gene Network"),
                 visNetworkOutput("miRNAGeneNetwork_Experimental", height = "700px"),
                 downloadButton("corr_dl_network_nodes", "Download nodes"),
                 downloadButton("corr_dl_network_edges", "Download edges"),
                 tags$hr(),
                 
                 ###################################### NEW CHANGES #####################################
                 fluidRow(
                   column(
                     2,
                          selectInput("corr_select_up_or_down_miRNA", "Select Gene Targets of Up- or Down-regulated miRNA",
                                      choices = c("Up" = "Up", "Down" = "Down"),
                                      selected = "Up")
                   ),
                   column(
                     2,
                     radioButtons(
                       inputId = "corr_select_supported",
                       label = "Select scope of miRNA–gene pairs for Pathway ORA",
                       choices = c(
                         "DB-supported pairs" = "Restrict",
                         "All negative correlations" = "Relax",
                         "Anticorrelated genes with padj cutoff" = "RestrictByPadj"
                       ),
                       selected = "Relax",
                       inline = FALSE
                     ),
                   ),
                   column(
                     5,
                     tags$div(
                       tags$strong("Data processing: "),
                       tags$br(),
                       "DB-supported pairs: Only miRNA–mRNA pairs that are anti-correlated AND supported by target prediction/validation databases.",
                       tags$br(),
                       "All negative correlations: All genes from anti-correlated miRNA–mRNA pairs, regardless of database support or p cutoff.",
                       tags$br(),
                       "Anticorrelated genes with padj cutoff: All genes from anti-correlated miRNA–mRNA pairs, regardless of database support, that meet a user-specified adjusted p-value threshold for correlation."
                     )
                   ),
                   column(
                    3,
                     conditionalPanel(
                       condition = "input.corr_select_supported == 'RestrictByPadj'",
                       numericInput(
                         inputId = "corr_padj_threshold",
                         label = "Adjusted p-value threshold for correlation",
                         value = 0.05,
                         min = 0,
                         max = 1,
                         step = 0.005
                       )
                     )
                   )
                 ),
                 fluidRow(
                   column(2,
                          # Background universe selector
                          selectInput(
                            "corr_select_background",
                            "Background universe for ORA",
                            choices = c(
                              "Organism_all",          # all ENTREZ IDs in OrgDb
                              "Custom_list",           # user-entered comma-separated ENTREZ IDs
                              "Experimental_list"    # union of miRTarBase targets (only when scope == validated_agg)
                            ),
                            selected = "Experimental_list"
                          )
                   ),
                   
                   column(6,
                          tags$div(
                            tags$p(strong("Background options and their biological questions:")),
                            tags$ul(
                              tags$li(
                                tags$strong("Organism_all: "),
                                "Are the target genes in anticorrelated miRNA-gene pairs enriched in specific biology relative to all annotated genes in the organism?"
                              ),
                              tags$li(
                                tags$strong("Custom_list: "),
                                "Are the target genes in anticorrelated miRNA-gene pairs enriched relative to a user-defined gene context (paste ENTREZ IDs below)?"
                              ),
                              tags$li(
                                tags$strong("Experimental_list: "),
                                "Are the target genes in anticorrelated miRNA-gene pairs enriched compared to all genes that were eligible for the correlation analysis? (recommended)",
                                em("We define these as the significant DE genes with absolute log2foldChange of >= 0.5.")
                              ),
                            )
                          )
                     
                   ),
                   column(2,
                          selectInput("corr_analysis_type", "Select Pathway Database",
                                      choices = c("GO" = "GO", "Reactome" = "Reactome"),
                                      selected = "Reactome")
                   ),
                   column(2, numericInput("corr_pathway_padj_cutoff", "FDR-adjusted P-value < ", value = 0.05, min = 0, max = 1, step = 0.01))
                 ),
                 fluidRow(
                   column(12,
                          # Show custom input only when Custom_list is selected
                          conditionalPanel(
                            condition = "input.corr_select_background == 'Custom_list'",
                            textInput(
                              "custom_universe_ids",
                              "Custom universe (comma-separated ENTREZ IDs)",
                              placeholder = "e.g., 7157, 1956, 7422"
                            ),
                            helpText("Paste ENTREZ gene IDs separated by commas. The analysis will test enrichment against this custom background.")
                          )
                   )
                 ),
                 actionButton("Pathway_ORA_Experimental", "Run Pathway Overrepresentation Analysis"),
                 DTOutput("corr_pathway_ORA_table"),
                 #tableOutput("corr_pathway_ORA_table"),
                 # New inputs for user-defined parameters
                 tags$hr(),
                 h6("The following max. pathway and min. coverage controls affect both the chord and network plots."),
                 h6("The number of pathways shown will be affected by coverage selection."),
                 numericInput("corr_num_pathways", "Maximum Number of Pathways to Show:", value = 10, min = 2),
                 h6("Coverage indicates the minimum fraction of genes in a pathway an miRNA must target to be shown."),
                 numericInput("corr_min_coverage", "Minimum Coverage (0 to 1):", value = 0.1, min = 0, max = 1, step = 0.01),
                 tags$hr(),
                 
                  actionButton("corr_showChordPlot", "Show Chord Plot"),
                 tags$div(
                   style = "overflow: visible; width: 100%; height: 600px;",
                   plotOutput("corr_chordPlot", width = "100%", height = "600px")),
                 downloadButton("corr_downloadChord", "Download Plot"),
                 downloadButton("corr_downloadChordData", "Download Data"),
                 h6("This will enable you to show only nodes that have at least n number connections"),
                 numericInput("corr_min_degree", "Minimum node degree:", value = 1, min = 0),
                 actionButton("corr_showNetworkPlot", "Show Network Plot"),
                 visNetworkOutput("corr_networkPlot", height = "700px"),
                 downloadButton("path_corr_dl_network_nodes", "Download nodes"),
                 downloadButton("path_corr_dl_network_edges", "Download edges")
                 ###############################################################################
        )
      )
    )
  )
)


