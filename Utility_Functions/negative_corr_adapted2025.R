#' This function has been minimally modified from https://rdrr.io/bioc/anamiR/src/R/negative_cor.R (GPL-2 license) and made to be compatible with miRQuest.
#' Please cite: Wang TT, Lee CY, Lai LC, Tsai MH, Lu TP, Chuang EY. anamiR: integrated analysis of MicroRNA and gene expression profiling. BMC Bioinformatics. 2019 May 14;20(1):239. doi: 10.1186/s12859-019-2870-x. PMID: 31088348; PMCID: PMC6518761.
#' @import stats
#' @export
negative_cor <- function(
  mrna_data,
  mirna_data,
  method = c("pearson",
             "kendall",
             "spearman"),
  cut.off = -0.5
) {

  method <- match.arg(method)

  common_column <- intersect(colnames(mrna_data),
                             colnames(mirna_data))

  mrna_data <- as.data.frame(mrna_data) %>% dplyr::select(all_of(common_column))

  mirna_data <- as.data.frame(mirna_data) %>% dplyr::select(all_of(common_column))

  cal_cor <- function(data_1, data_2, cor_cut) {
    n <- 1
    corr <- list()
    
    withProgress(message = 'Calculating correlations...', value = 0, {
      for (i in seq_len(nrow(data_1))) {
        for (j in seq_len(nrow(data_2))) {
          mrna <- as.numeric(data_1[i, 1:(ncol(mrna_data) - 5)])
          mirna <- as.numeric(data_2[j, 1:(ncol(mirna_data) - 5)])
           # tmp <- stats::cor(mrna, mirna, method = method)
          
          valid <- !(is.na(mrna) | is.na(mirna))    
          x <- mrna[valid]                          
          y <- mirna[valid]                         
          
          # Ensure enough data points for cor.test and use cor.test
          if (length(x) >= 3) {                  
            
            test <- try(
              suppressWarnings(                    # NEW: suppress warnings from cor.test
                stats::cor.test(x, y, method = method)
              ),
              silent = TRUE
            ) 
            
            if (inherits(test, "try-error")) {                                     
              next                                                                  
            }                                                                       
            tmp  <- unname(test$estimate)           
            pval <- test$p.value                     
            
            #print(paste0("The row indices are ", i, " for mrna, and ", j, " for mirna"))
            #print(tmp)
            
            if (tmp < cor_cut) {
              # MODIFIED: Insert pval as a new element (4th position) in the stored result
              corr[[n]] <- c(
                row.names(data_2)[j],               # miRNA
                row.names(data_1)[i],               # Gene
                tmp,                                # correlation
                pval,                               # NEW: correlation p-value
                data_2[["log_ratio"]][j],
                data_2[["P-adjust"]][j],
                data_2[["mean_case"]][j],
                data_2[["mean_control"]][j],
                data_1[["log_ratio"]][i],
                data_1[["P-adjust"]][i],
                data_1[["mean_case"]][i],
                data_1[["mean_control"]][i]
              )
              n <- n + 1
            }
          }  
        
        }
        # Update progress based on the mrna_data row
        incProgress(1/ nrow(data_1), detail = paste("Correlating mRNA", i, "of", nrow(data_1), "with", nrow(data_2), "miRNAs"))
        
      }
    })
    return(corr)
  }

  corr <- cal_cor(mrna_data, mirna_data, cut.off)

  corr <- if (length(corr)) do.call(rbind, corr) else NULL    

  if (is.null(corr)) {
    cut.off <- cut.off + 0.2
    corr <- cal_cor(mrna_data, mirna_data, cut.off)
    corr <- if (length(corr)) do.call(rbind, corr) else NULL
  }
  
  # If still NULL, return empty data frame with expected columns
  if (is.null(corr)) {                                        # NEW
    last_column <- "Correlation"                              # NEW
    colnames_out <- c(                                        # NEW
      "miRNA", "Gene", last_column, "P-value(correlation)",   # NEW: added p-value column name
      "logratio_miRNA", "P-adjust(miRNA)",
      "mean_case(miRNA)", "mean_control(miRNA)",
      "logratio_gene", "P-adjust(gene)",
      "mean_case(gene)", "mean_control(gene)"
    )
    return(data.frame(matrix(ncol = length(colnames_out), nrow = 0,
                             dimnames = list(NULL, colnames_out))))   # NEW
  }
  
  # Convert to data.frame and set column names
  corr <- as.data.frame(corr, stringsAsFactors = FALSE)      
  
  
  last_column <- paste("Correlation")
  colnames(corr) <- c("miRNA", "Gene", last_column, "P-value(correlation)", 
                      "logratio_miRNA", "P-adjust(miRNA)",
                      "mean_case(miRNA)", "mean_control(miRNA)",
                      "logratio_gene", "P-adjust(gene)",
                      "mean_case(gene)", "mean_control(gene)")
  
  # coerce numeric columns appropriately (since we built with c())
  numeric_cols <- c(last_column, "P-value(correlation)",
                    "logratio_miRNA", "P-adjust(miRNA)",
                    "mean_case(miRNA)", "mean_control(miRNA)",
                    "logratio_gene", "P-adjust(gene)",
                    "mean_case(gene)", "mean_control(gene)")
  corr[numeric_cols] <- lapply(corr[numeric_cols], function(z) suppressWarnings(as.numeric(z)))  # MODIFIED
  
  return(corr)
}


