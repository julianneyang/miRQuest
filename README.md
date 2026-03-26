# miRQuest: Interactive Analysis and Visualization of miRNA Sequencing Data

## Purpose

MiRQuest is an RShiny application for the analysis of microRNA-sequencing data. MicroRNA (miRNA) are small noncoding RNA that repress translation or induce degradation of target messenger RNA (mRNA). The app features common visualization methods for high-throughput sequencing data, such as multidimensional scaling, stacked column charts, heatmaps, and boxplots to compare expression across groups for a user-specified miRNA of interest. MiRQuest further provides support for differential miRNA and gene expression analysis. Unique to miRNome sequencing data analysis, users may retrieve predicted gene targets of differentially expressed miRNA and follow up with pathway overrepresentation analysis of the gene targets. Finally, if users upload paired bulk mRNA sequencing data, they may identify differentially expressed genes and further retrieve negatively correlated miRNA-gene pairs.

## Workflow Overview

Please view the Quickstart tab in the app or read the tutorial in the `www` folder.

## Reproducibility

MiRQuest provides **two complementary approaches** for ensuring computational reproducibility:

### Option 1: renv (R Package Management)
Use `renv` to restore the exact R package versions used in this project. This approach works well if you have R 4.4.2 already installed and want to work within your existing R environment.

**Best for:**
- Users with R 4.4.2 installed
- Active development and modification of the code
- When you need to integrate with other R projects

### Option 2: Docker (Complete Environment Isolation)
Use Docker to recreate the **entire computational environment**, including R version, system libraries, and all dependencies in a self-contained container.

**Best for:**
- **Publications and long-term reproducibility** (strongly recommended)
- Users with R versions < 4.4.2
- Avoiding conflicts with existing R installations
- Cross-platform consistency (Windows/Mac/Linux)
- Archiving computational environments with publications

**We strongly recommend Docker (Method 3 below) for publication-quality reproducibility**, as it captures not just R package versions but the complete system environment. However, renv (Methods 1-2 below) is suitable for development and quick testing if you already have R 4.4.2 installed.

### Using renv Lock Files for Different Operating Systems

To ensure proper package installation on your operating system, use the appropriate renv lock file:

- **Windows**: Use `renv_for_Windows.lock`
  ```R
  renv::restore(lockfile = "renv_for_Windows.lock")
  ```

- **macOS**: Use `renv_for_Mac.lock`
  ```R
  renv::restore(lockfile = "renv_for_Mac.lock")
  ```

- **Linux**: Use `renv.lock` (default) or `renv_lock_for_Linux.lock`
  ```R
  renv::restore()  # Uses renv.lock by default
  # OR
  renv::restore(lockfile = "renv_lock_for_Linux.lock")
  ```

If you encounter issues with package installation, you may need to manually specify the lock file for your operating system.

```R
> R.version
               _                                
platform       x86_64-w64-mingw32               
arch           x86_64                           
os             mingw32                          
crt            ucrt                             
system         x86_64, mingw32                  
status                                          
major          4                                
minor          4.4                              
year           2024                             
month          10                               
day            31                               
svn rev        87279                            
language       R                                
version.string R version 4.4.2 (2025-06-03 ucrt)
nickname       Pile of Leaves
```

## Installation Instructions

### Method 1: Clone github repository and double click on .bat file

*[Recommended for Windows without computational experience]*
- To clone without Git: Click the green `Code` button at the upper right of the repository, then click `Download ZIP`.
- Extract the files from the zipped files.
- Double click on the `Create_Exe_File.R` script. This should automatically launch RStudio; if prompted, select "Open With RStudio". Click "Run".
- Close RStudio. Double click on the newly created "miRQuest.bat" file.
- To stop the app, close the terminal.
- Assumptions
  - Must be a Windows User
  - Have a preexisting installation of R, Rstudio, and Rtools (matching version of R; e.g. if running R 4.4.2, install Rtools 44)

### Method 2: Clone github repository, open .Rproj file, and click Run App

*[Recommended for non-Windows users without computational experience]*
- To clone without Git: Click the green `Code` button at the upper right of the repository, then click `Download ZIP`.
- Extract the files from the zipped files.
- Double click on the `miRnome-Rshiny.Rproj` file. This will launch RStudio.
- In the upper left, go to File > Open > `global.R`
- Click the `Run App` button.
- To stop the app, close the app window or click the "STOP" sign in RStudio.
- Assumptions
  - Have a preexisting installation of R, RStudio, and Rtools

### Method 3: Docker Deployment (Recommended for Publications)

*[**STRONGLY RECOMMENDED** for publications, reproducibility, and archival purposes]*

Docker ensures exact reproducibility by encapsulating R version, all packages, system libraries, and dependencies in a container that works identically across Windows, macOS, and Linux. Images can be archived alongside publications and run identically years later.

#### Prerequisites

- Docker Engine 20.10+ ([Docker Desktop](https://www.docker.com/products/docker-desktop/) for Windows/Mac)
- 4GB+ RAM (8GB recommended)
- 10GB+ free disk space
- Verify installation: `docker --version && docker run hello-world`

#### Quick Start

**Option A: Docker Compose (Easiest)**
```bash
git clone https://github.com/MSDLLCpapers/miRQuest.git
cd miRQuest
docker-compose up -d          # Build and start (15-30 min first time)
# Access at http://localhost:3838
docker-compose down           # Stop
```

**Option B: Standard Docker**
```bash
git clone https://github.com/MSDLLCpapers/miRQuest.git
cd miRQuest
docker build -t mirquest:v1.0.0 .                              # Build image (15-30 min)
docker run -d -p 3838:3838 --name mirquest-app mirquest:v1.0.0 # Start
# Access at http://localhost:3838
docker logs -f mirquest-app                          # View logs
docker stop mirquest-app                             # Stop
docker start mirquest-app                            # Restart
docker rm mirquest-app                               # Remove
```

**Build Process:** Base image download (2-5 min) → System dependencies (3-5 min) → R packages (10-20 min) → App files (1 min)

#### Configuration

**Port Configuration:**
By default, the app runs on port 3838. To use a different port:
```bash
# Docker Compose: Edit docker-compose.yml
ports:
  - "8080:3838"  # Access at http://localhost:8080

# Docker CLI:
docker run -d -p 8080:3838 --name mirquest-app mirquest:v1.0.0
```

**Persistent Data Storage:**
To save analysis outputs or uploads to your host machine:
```bash
# Docker Compose: Uncomment volumes in docker-compose.yml
volumes:
  - ./output:/app/output

# Docker CLI:
docker run -d -p 3838:3838 -v $(pwd)/output:/app/output mirquest:v1.0.0
```

**Development Mode:**
For development with live code reloading:
```bash
docker run -d -p 3838:3838 -v $(pwd):/app mirquest:v1.0.0
```

#### Troubleshooting

**Build Issues:**
```bash
# Memory issues: Increase Docker memory allocation (recommended: 4GB+)
# Docker Desktop: Settings > Resources > Memory

# Build fails: Clear cache and rebuild
docker builder prune -a && docker build --no-cache -t mirquest:v1.0.0 .

# Network timeouts: Retry the build
docker-compose build --no-cache

# Timeout errors: Disable BuildKit
DOCKER_BUILDKIT=0 docker build -t mirquest:v1.0.0 .

# Check build logs
docker-compose build 2>&1 | tee build.log
```

**Runtime Issues:**
```bash
# Container starts but app doesn't load:
docker logs mirquest-app  # Check logs
lsof -i :3838  # Verify port not in use (macOS/Linux)
netstat -ano | findstr :3838  # Windows

# Port conflict: Use different port
docker run -d -p 3839:3838 --name mirquest-app mirquest:v1.0.0

# Out of memory errors: Increase Docker memory limit to at least 4GB

# Biomart connection errors: App uses "useast" Ensembl mirror
# Usually transient; refresh the app or restart the container
```

#### Archiving for Publications

The Docker image associated with the paper *"MiRQuest: A user-friendly web app for the interactive analysis and visualization of microRNA sequencing data"* is tagged as **`mirquest:v1.0.0`** (label: `mirquest_v1.0.0`). This explicit tag ensures that the exact computational environment described in the paper can be rebuilt and archived independently of the `:latest` tag, which may change over time.

**Build the publication image:**
```bash
git checkout v1.0.0                                  # Use the release tag
docker build -t mirquest:v1.0.0 .                    # Build with explicit version tag
```

**Save Docker image for long-term archival:**
```bash
docker save mirquest:v1.0.0 | gzip > mirquest-v1.0.0.tar.gz  # ~1.5-2.5 GB
sha256sum mirquest-v1.0.0.tar.gz > checksums.txt
```

**Load archived image to reproduce results:**
```bash
docker load < mirquest-v1.0.0.tar.gz
docker run -d -p 3838:3838 --name mirquest-app mirquest:v1.0.0
```

**Long-Term Archiving Guidance:**

For publication-grade reproducibility, we recommend archiving the Docker image in a persistent repository that provides a DOI:

1. **[Zenodo](https://zenodo.org)** (recommended) - Free, DOI-minting, supports files up to 50 GB, backed by CERN. Upload `mirquest-v1.0.0.tar.gz` and `checksums.txt`.
2. **[Figshare](https://figshare.com)** - Free, DOI-minting, supports files up to 20 GB per file.
3. **Institutional data repository** - Many universities and research institutions host long-term data archives.

When archiving, include the following alongside the Docker image tarball:
- `checksums.txt` (SHA-256 checksum for image integrity verification)
- This `README.md`
- The `Dockerfile`
- The `renv.lock` file (R package versions)
- Git release tag or commit SHA (e.g., `v1.0.0`)
- Package version report (see Verifying Environment below)

**Example Methods Section Citation:**
> "Data analysis was performed using miRQuest v1.0.0 (Yang & Sauter et al.; GitHub: MSDLLCpapers/miRQuest, tag: v1.0.0) in a Docker container (archived at Zenodo: [DOI]). The image is based on Bioconductor RELEASE_3_20, containing R 4.4.2, Bioconductor 3.20, and exact package versions specified in the Dockerfile and renv.lock."

#### Verifying Environment

Generate a package version report for publication supplementary materials:

```bash
docker run -d -p 3838:3838 --name mirquest-app mirquest:v1.0.0
docker exec mirquest-app Rscript -e "write.csv(installed.packages()[,c('Package','Version')], '/tmp/packages.csv')"
docker cp mirquest-app:/tmp/packages.csv ./mirquest-package-versions.csv
```

**Expected versions:** R 4.4.2, Bioconductor 3.20, DESeq2 1.44.0, miRNAtap 1.38.0, clusterProfiler 4.12.0

#### Advantages of Docker Deployment

- ✅ **Complete isolation** - No conflicts with existing R installations
- ✅ **Reproducibility** - Exact same environment every time
- ✅ **Version control** - All package versions locked to R 4.4.2 and Bioconductor 3.20
- ✅ **Cross-platform** - Works identically on Windows, macOS, and Linux
- ✅ **Publication-ready** - Docker images can be archived with publications
- ✅ **No R installation required** - Only Docker needed

**Image Size:** The final Docker image is approximately 3-4 GB (R base ~700 MB, system libraries ~500 MB, R packages ~2-3 GB).

---

## Sample Data

The demo dataset included in `inst/data` was derived from TCGA-COAD and TCGA-READ colorectal cancer data. The script used to create it is available at `Data/make_demo_dataset.R`. Note that this script requires additional R packages **not included** in the app's `renv.lock` or Docker image, as they are only needed for data preparation and not for running miRQuest itself:

- `TCGAbiolinks` - TCGA data retrieval
- `SummarizedExperiment` - Bioconductor data structures
- `miRBaseConverter` - miRNA ID conversion to mature miRBase v22
- `caret` - classification utilities
- `scutr` - Tomek link undersampling

To install these packages for recreating the demo dataset:
```R
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "miRBaseConverter"))
install.packages(c("caret", "scutr"))
```

See `inst/data` to view how sample input files should look like.

**COREAD_Downsampled_Tumor_vs_Normal_Mature_miRNA.csv**

- This is the miRNA input file
- Requires: raw counts, formatted as a .csv, with features as mature miRNA IDs from miRBase v.22

**COREAD_Tumor_vs_Normal_Metadata.csv**

- This is the metadata input file
- Requires: .csv, with SampleID as the first column header

**COREAD_Downsampled_Tumor_vs_Normal_mRNA.csv**

- This is the RNA-seq input file
- Requires: raw counts, formatted as a .csv, with features as ENSEMBL IDs
---

## Development: Modular Server Architecture

This application uses a modular server architecture to improve code maintainability and development workflow. The server logic has been separated into focused, functional modules while maintaining all reactive relationships.

### File Structure

#### Core Files
- `server.R` - Main server file with reactive variable definitions and module sourcing
- `global.R` - Libraries, settings, and helper function sourcing
- `ui.R` - User interface definition

#### Server Code Modules (`Server_Code/`)

The server logic is organized into the following modules:

1. **`data_input_handlers.R`**
  - Module 0A: miRNA counts data reading and processing
  - Module 0B: Metadata file handling and UI updates
  - Species validation logic

2. **`data_subsetting.R`**
  - Module 2: Data filtering and subsetting functionality
  - Dynamic UI generation for filter values
  - Download handlers for subsetted data

3. **`exploratory_visualization.R`**
  - Module 1: Exploratory data visualization
  - Stacked column charts, MDS plots, heatmaps
  - Data normalization and transformation

4. **`differential_mirna_analysis.R`**
  - Module 3: DESeq2 differential miRNA expression analysis
  - Volcano plots, significance testing
  - Alternative limma::voom analysis
  - Boxplot generation for individual miRNAs

5. **`differential_gene_analysis.R`**
  - Module 6: DESeq2 differential gene expression analysis
  - mRNA data processing and analysis
  - Results filtering and visualization

6. **`correlation_analysis.R`**
  - Module 7: miRNA-mRNA correlation analysis
  - Negative correlation identification
  - Database support validation
  - Correlation plotting and filtering

7. **`target_prediction_all.R`**
  - Module 4: Bulk target gene prediction for all DE miRNAs
  - Pathway overrepresentation analysis
  - Up/down-regulated miRNA separation
  - Progress tracking for large analyses

8. **`target_prediction_single.R`**
  - Module 5: Single miRNA target prediction
  - Individual pathway analysis
  - Dotplot and barplot generation

9. **`pathway_visualization.R`**
  - Chord plot generation
  - Network plot visualization
  - Interactive pathway-miRNA relationships

10. **`future_expansion_code.R`**
  - Commented-out code for future development
  - Alternative implementations
  - Experimental features

### Reactive Variables

All reactive variables are defined in the main `server.R` file and are accessible across all modules:

#### Species and Database Reactives
- `mart()` - Biomart connection
- `org_string()` - Organism-specific gene symbols
- `org_abbrev()` - Species abbreviation
- `org_db()` - Organism database
- `organism_reactive()` - Full organism name

#### Data Storage Reactives
- `longData()`, `wideData()` - miRNA count data
- `metadata()` - Sample metadata
- `DESEQ_obj()` - DESeq2 object
- `res_significant()` - Significant DE miRNAs
- `mrna_res_significant()` - Significant DE genes

#### Analysis Results Reactives
- `predicted_target_reactive()` - All predicted targets
- `neg_cor_reactive()` - Negative correlations
- `All_miRNA_Pathways_reactive()` - Pathway analysis results
- And many more...

### Key Architectural Features

#### Tidyverse Code Style
- Proper indentation and spacing
- Function arguments on separate lines for clarity
- Consistent naming conventions

#### Error Handling
- Notification systems for long-running processes
- Validation checks for required inputs
- Graceful handling of missing data

#### Download Capabilities
- CSV downloads for all major results
- Plot downloads in multiple formats
- Consistent file naming with timestamps

#### Progress Tracking
- Progress bars for intensive computations
- User feedback during database queries
- Notification management

### Module Loading

The modular structure is automatically sourced in `server.R`:

```r
files_to_source <- list.files("Server_Code", pattern = "\\.R$", ignore.case = TRUE)

for (file in files_to_source) {
    source(paste0("Server_Code/", file), local = TRUE)
}
```

This approach maintains all reactive relationships while providing clean separation of concerns. Each module focuses on specific functionality while having access to all shared reactive variables. The pattern matching ensures only R files are sourced, avoiding non-code files like README documents.

### Why Separate miRNA and mRNA Analysis?

The application includes two distinct DESeq2 workflows:

#### miRNA Analysis Pipeline
- **Input**: `input$countTable` (miRNA counts)
- **Purpose**: Find differentially expressed miRNAs
- **Outputs**: Used for target prediction, pathway analysis, visualization
- **Downstream**: Feeds into correlation analysis, single miRNA analysis

#### mRNA Analysis Pipeline
- **Input**: `input$mrna_countTable` (gene expression counts)
- **Purpose**: Find differentially expressed genes
- **Outputs**: Used specifically for correlation with miRNA results
- **Downstream**: Only used in Module 7 correlation analysis

This separation allows analysis of two different molecular data types from the same samples to understand their regulatory relationships.

### Benefits of Modular Architecture

1. **Maintainability** - Each module handles specific functionality
2. **Readability** - Code is organized by biological purpose
3. **Debugging** - Issues can be isolated to specific modules
4. **Collaboration** - Multiple developers can work on different modules
5. **Testing** - Individual modules can be tested independently
6. **Scalability** - New features can be added as new modules

### Dependencies

Each module assumes access to:
- All reactive variables defined in `server.R`
- All libraries loaded in `global.R`
- Helper functions sourced from the `scripts/` directory
- Standard Shiny server environment (`input`, `output`, `session`)
---

## Contributing

When contributing to this project, please:

1. Follow the established modular structure
2. Add new functionality to the appropriate module
3. Maintain tidyverse coding standards
4. Test individual modules before integration
5. Update documentation as needed

## Contributors

- Julianne Yang - <juliannecyang@gmail.com>, <julianne.yang@msd.com>
- Jake Sauter - <jake.sauter3@gmail.com>, <jake.sauter@msd.com>

## Contact

For questions, issues, or suggestions, please reach out to the contributors via email or create an issue on GitHub.
