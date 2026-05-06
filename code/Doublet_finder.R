library(Seurat)
library(ggplot2)
library(tidyr)
library(presto)
library(DropletUtils)
library(DoubletFinder)

# -- Directories ---------------------------------------------------------------
# Resolve repo root: REPO_DIR env var > here::here() > getwd().
# Run this script from the repo root, or set REPO_DIR explicitly.
base_dir <- Sys.getenv("REPO_DIR", unset = NA)
if (is.na(base_dir)) {
  base_dir <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
}
out_dir  <- file.path(base_dir, "data")
plot_dir <- file.path(base_dir, "plot", "ctrl")
cb_dir   <- file.path(base_dir, "data", "cellbender")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# ==============================================================================
# 1. Load CellBender-filtered data + DoubletFinder per sample
# ==============================================================================

# Picks up both layouts:
#   data/cellbender/<name>_cellbender_filtered.h5                       (Zenodo tar.gz, flat)
#   data/cellbender/<name>_cellbender/<name>_cellbender_filtered.h5     (per-sample run)
h5_files     <- list.files(cb_dir, pattern = "_filtered\\.h5$",
                           full.names = TRUE, recursive = TRUE)
sample_names <- gsub("_cellbender_filtered\\.h5$", "", basename(h5_files))
names(h5_files) <- sample_names

sobj_list <- lapply(sample_names, function(s) {
  cat(sprintf("=== %s ===\n", s))
  
  # Load
  sce <- read10xCounts(h5_files[s])
  dat <- as(counts(sce), "dgCMatrix")
  colnames(dat) <- colData(sce)$Barcode
  rownames(dat) <- make.unique(rowData(sce)$Symbol)
  
  obj <- CreateSeuratObject(
    counts       = dat,
    project      = s,
    min.cells    = 5,
    min.features = 200
  )
  obj$orig.ident    <- s
  obj[["ensembl_id"]] <- rowData(sce)$ID[match(rownames(obj),
                                               make.unique(rowData(sce)$Symbol))]
  
  # Light QC for DoubletFinder (permissive)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- subset(obj, nCount_RNA > 100  & nFeature_RNA > 100  & percent.mt < 15)
  
  # Preprocessing
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, seed.use = 42, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:10, seed.use = 42, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.3, random.seed = 42, verbose = FALSE)
  
  # pK sweep
  sweep_res   <- paramSweep(obj, PCs = 1:10, sct = FALSE, num.cores=1)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  pk_res      <- find.pK(sweep_stats)
  pk_best     <- as.numeric(as.character(
    pk_res$pK[which.max(pk_res$BCmetric)]
  ))
  
  # Homotypic adjustment (5k loaded = ~4% rate)
  homotypic_prop <- modelHomotypic(obj$seurat_clusters)
  nExp           <- round(0.04 * ncol(obj))
  nExp_adj       <- round(nExp * (1 - homotypic_prop))
  
  cat(sprintf("  Cells: %d | pK: %.4f | nExp_adj: %d\n",
              ncol(obj), pk_best, nExp_adj))
  
  # Run DoubletFinder
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25,
                       pK = pk_best, nExp = nExp_adj, sct = FALSE)
  
  # Clean up metadata
  df_class_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  df_pann_col  <- grep("^pANN", colnames(obj@meta.data), value = TRUE)
  
  obj$doublet_class <- obj@meta.data[[df_class_col]]
  obj$doublet_score <- obj@meta.data[[df_pann_col]]
  obj@meta.data[[df_class_col]] <- NULL
  obj@meta.data[[df_pann_col]]  <- NULL
  
  cat(sprintf("  Doublets: %d (%.1f%%)\n",
              sum(obj$doublet_class == "Doublet"),
              100 * mean(obj$doublet_class == "Doublet")))
  
  # Strip DoubletFinder preprocessing (keep only counts + metadata)
  # so it doesn't interfere with the real integration later
  DefaultAssay(obj) <- "RNA"
  obj[["pca"]]  <- NULL
  obj[["umap"]] <- NULL
  
  obj
})
names(sobj_list) <- sample_names

# Save per-sample doublet stats
doublet_summary <- do.call(rbind, lapply(sobj_list, function(obj) {
  data.frame(
    sample    = obj$orig.ident[1],
    n_cells   = ncol(obj),
    n_singlet = sum(obj$doublet_class == "Singlet"),
    n_doublet = sum(obj$doublet_class == "Doublet"),
    pct_doublet = 100 * mean(obj$doublet_class == "Doublet")
  )
}))
write.csv(doublet_summary, file.path(out_dir, "doublet_summary.csv"), row.names = FALSE)

# -- Merge all samples ---------------------------------------------------------
dat_big <- merge(
  sobj_list[[1]],
  y            = sobj_list[-1],
  add.cell.ids = sample_names,
  project      = "merged"
)

# Check doublet distribution
table(dat_big$doublet_class)

saveRDS(dat_big, file = file.path(out_dir, "dat_big_filtered_cellbender_doubletfinder.rds"))
 