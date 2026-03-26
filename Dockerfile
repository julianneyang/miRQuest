# Use Bioconductor Docker image with R 4.4.2 and Bioconductor 3.20
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Metadata labels
LABEL org.opencontainers.image.title="miRQuest"
LABEL org.opencontainers.image.description="miRQuest - Interactive microRNA Analysis Shiny Application"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.authors="Julianne Yang <juliannecyang@gmail.com>, Jake Sauter <jake.sauter3@gmail.com>"
LABEL org.opencontainers.image.source="https://github.com/MSDLLCpapers/miRQuest"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.ref.name="mirquest_v1.0.0"
LABEL org.opencontainers.image.documentation="Associated with: Yang & Sauter et al., 'MiRQuest: A user-friendly web app for the interactive analysis and visualization of microRNA sequencing data'"

# Set working directory
WORKDIR /app

# Install system dependencies required for R packages
# These are commonly needed for Bioconductor and other packages
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libgsl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libicu-dev \
    libmariadb-dev \
    libpq-dev \
    libsqlite3-dev \
    libv8-dev \
    libcairo2-dev \
    libxt-dev \
    pandoc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables for R and SSL certificates
ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt

# Install required Bioconductor and CRAN packages using BiocManager
# This ensures proper dependency resolution for Bioconductor packages
RUN R -e "options(repos = BiocManager::repositories()); \
    BiocManager::install(c( \
    'DESeq2', 'edgeR', 'limma', \
    'miRNAtap', 'miRNAtap.db', \
    'clusterProfiler', 'ReactomePA', 'topGO', 'fgsea', 'enrichplot', 'DOSE', 'GOSemSim', \
    'biomaRt', 'org.Hs.eg.db', 'org.Mm.eg.db', 'GO.db', 'reactome.db', \
    'graph', 'graphite', 'EnhancedVolcano', \
    'shiny', 'tidyverse', 'dplyr', 'pheatmap', 'ggplot2', 'ggrepel', 'ggtree', 'ggraph', \
    'data.table', 'matrixStats', 'sqldf', 'conflicted', 'paletteer', \
    'circlize', 'patchwork', 'cowplot', 'scatterpie', \
    'renv', 'vroom', 'here', 'shinyjs', 'bslib', 'igraph', \
    'DT', 'plotly', 'visNetwork' \
    ), ask = FALSE, update = FALSE)"

# Copy the application files
COPY global.R global.R
COPY server.R server.R
COPY ui.R ui.R
COPY Server_Code/ Server_Code/
COPY Utility_Functions/ Utility_Functions/
COPY inst/ inst/
COPY www/ www/

# Expose port for Shiny app (default Shiny port is 3838)
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp(host='0.0.0.0', port=3838)"]
