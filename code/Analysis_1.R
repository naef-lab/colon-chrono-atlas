
# ==============================================================================
# 0. Setup
# ==============================================================================

library(Seurat)
library(ggplot2)
library(tidyr)
library(presto)
library(DropletUtils)
library(patchwork)

make_contour <- function(seurat_obj, xlim, ylim, color = "grey60", n_smooth = 300) {
  umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_df) <- c("x", "y")
  umap_df <- umap_df[umap_df$x >= xlim[1] & umap_df$x <= xlim[2] &
                       umap_df$y >= ylim[1] & umap_df$y <= ylim[2], ]
  
  hull_idx <- chull(umap_df$x, umap_df$y)
  hull_df <- umap_df[c(hull_idx, hull_idx[1]), ]
  
  # Parametric spline smoothing
  t_orig <- seq_len(nrow(hull_df))
  t_smooth <- seq(1, nrow(hull_df), length.out = n_smooth)
  smooth_df <- data.frame(
    x = spline(t_orig, hull_df$x, xout = t_smooth, method = "periodic")$y,
    y = spline(t_orig, hull_df$y, xout = t_smooth, method = "periodic")$y
  )
  
  geom_path(
    data        = smooth_df,
    aes(x = x, y = y),
    color       = color,
    linetype    = "dashed",
    linewidth   = 0.3,
    inherit.aes = FALSE
  )
}

common_theme <-
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        aspect.ratio = 1,
        plot.title   = element_text(size = 12, face = "italic", hjust = 0.5),
        legend.title = element_blank(),
        legend.text  = element_text(size = 8),
        legend.key.height = unit(0.6, "cm"),
        legend.key.width  = unit(0.3, "cm")
  )


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
# 1. Load CellBender-filtered, doubletfinder-ran data
# ==============================================================================

dat_big <- readRDS(file.path(out_dir, "dat_big_filtered_cellbender_doubletfinder.rds"))
dat_big <- dat_big[,dat_big$doublet_class == "Singlet"]

# ==============================================================================
# 2. Sample name mapping and metadata
# ==============================================================================

# -- Standardized sample names -------------------------------------------------
name_map <- c(
  # CTRL (wildtype, distal colon only)
  "M1_ZT1"  = "M1_ZT01_DC_CTRL",
  "M2_ZT7"  = "M2_ZT07_DC_CTRL",
  "M3_ZT13" = "M3_ZT13_DC_CTRL",
  "M4_ZT19" = "M4_ZT19_DC_CTRL",
  "F1_ZT1"  = "F1_ZT01_DC_CTRL",
  "F2_ZT7"  = "F2_ZT07_DC_CTRL",
  "F3_ZT13" = "F3_ZT13_DC_CTRL",
  "F4_ZT19" = "F4_ZT19_DC_CTRL",
  # DSS (distal colon)
  "P47F4_10R_ZT07_CD" = "F5_ZT07_DC_DSS",
  "P47F5_10R_ZT01_CD" = "F6_ZT01_DC_DSS",
  "P48M2_10R_ZT19_CD" = "M5_ZT19_DC_DSS",
  "P48M3_10R_ZT01_CD" = "M6_ZT01_DC_DSS",
  "P48M5_10R_ZT07_CD" = "M7_ZT07_DC_DSS",
  "P49F1_10D_ZT13_CD" = "F7_ZT13_DC_DSS",
  "P49F2_10R_ZT19_CD" = "F8_ZT19_DC_DSS",
  "P50M4_10D_ZT13_CD" = "M8_ZT13_DC_DSS",
  # DSS (proximal colon)
  "P47F4_10R_ZT07_CP" = "F5_ZT07_PC_DSS",
  "P47F5_10R_ZT01_CP" = "F6_ZT01_PC_DSS",
  "P48M2_10R_ZT19_CP" = "M5_ZT19_PC_DSS",
  "P48M3_10R_ZT01_CP" = "M6_ZT01_PC_DSS",
  "P48M5_10R_ZT07_CP" = "M7_ZT07_PC_DSS",
  "P49F1_10D_ZT13_CP" = "F7_ZT13_PC_DSS",
  "P49F2_10R_ZT19_CP" = "F8_ZT19_PC_DSS",
  "P50M4_10D_ZT13_CP" = "M8_ZT13_PC_DSS"
)

# Clean raw suffixes from orig.ident
dat_big$orig.ident <- gsub("_(sample_)?raw_feature_bc_matrix", "",
                           dat_big$orig.ident)

# Map to standardized names
dat_big@meta.data$sample_name <- name_map[as.character(dat_big$orig.ident)]

# -- Factor ordering -----------------------------------------------------------
sample_order <- c(
  # CTRL DC (by ZT)
  "F1_ZT01_DC_CTRL", "M1_ZT01_DC_CTRL",
  "F2_ZT07_DC_CTRL", "M2_ZT07_DC_CTRL",
  "F3_ZT13_DC_CTRL", "M3_ZT13_DC_CTRL",
  "F4_ZT19_DC_CTRL", "M4_ZT19_DC_CTRL",
  # DSS DC (by ZT)
  "F6_ZT01_DC_DSS",  "M6_ZT01_DC_DSS",
  "F5_ZT07_DC_DSS",  "M7_ZT07_DC_DSS",
  "F7_ZT13_DC_DSS",  "M8_ZT13_DC_DSS",
  "F8_ZT19_DC_DSS",  "M5_ZT19_DC_DSS",
  # DSS PC (by ZT)
  "F6_ZT01_PC_DSS",  "M6_ZT01_PC_DSS",
  "F5_ZT07_PC_DSS",  "M7_ZT07_PC_DSS",
  "F7_ZT13_PC_DSS",  "M8_ZT13_PC_DSS",
  "F8_ZT19_PC_DSS",  "M5_ZT19_PC_DSS"
)

dat_big$sample_name <- factor(dat_big$sample_name, levels = sample_order)

# -- Parse metadata from structured sample names ------------------------------
dat_big$sex       <- gsub("^([MF]).*", "\\1", dat_big$sample_name)
dat_big$mouse     <- gsub("^([MF][0-9]+)_.*", "\\1", dat_big$sample_name)
dat_big$ZT        <- as.numeric(gsub(".*ZT([0-9]+).*", "\\1", dat_big$sample_name))
dat_big$region    <- ifelse(grepl("_DC_", dat_big$sample_name), "DC", "PC")
dat_big$condition <- ifelse(grepl("CTRL", dat_big$sample_name), "CTRL", "DSS")


# ==============================================================================
# 3. Quality control
# ==============================================================================

# -- Compute QC metrics --------------------------------------------------------
dat_big[["percent.mt"]] <- PercentageFeatureSet(dat_big, pattern = "^mt-")
dat_big[["percent.RP"]] <- PercentageFeatureSet(dat_big, pattern = "^Rp[sl]")

# -- QC thresholds -------------------------------------------------------------
qc_thresholds <- list(
  nCount_RNA   = c(lower = 300,  upper = 15000),
  nFeature_RNA = c(lower = 300,  upper = 5000),
  percent.mt   = c(lower = 0,    upper = 10)
)

# -- QC violin plots (pre-filtering) ------------------------------------------
region_colors <- c("DC" = "#E74C3C", "PC" = "#3498DB")

pdf(file.path(plot_dir, "violplot_QC_prefiltering_thresholds_cellbender.pdf"),
    width = 16, height = 6)

for (feat in names(qc_thresholds)) {
  th <- qc_thresholds[[feat]]
  
  p <- ggplot(dat_big@meta.data,
              aes(x = sample_name, y = .data[[feat]] + 0.1, fill = region)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_jitter(size = 0.1, alpha = 0.05, width = 0.3) +
    scale_fill_manual(values = region_colors) +
    scale_y_log10() +
    geom_vline(xintercept = c(8.5, 16.5), linetype = "dashed", color = "grey50") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 8)
    ) +
    ggtitle(feat)
  
  if (feat != "percent.mt") {
    p <- p +
      geom_hline(yintercept = th["lower"], linetype = "dashed", color = "black") +
      geom_hline(yintercept = th["upper"], linetype = "dashed", color = "black")
  } else {
    p <- p +
      geom_hline(yintercept = th["upper"], linetype = "dashed", color = "black")
  }
  
  print(p)
}

dev.off()

# -- Apply filters -------------------------------------------------------------
keep_cells <- dat_big$nCount_RNA   > qc_thresholds$nCount_RNA["lower"]   &
  dat_big$nCount_RNA   < qc_thresholds$nCount_RNA["upper"]   &
  dat_big$nFeature_RNA > qc_thresholds$nFeature_RNA["lower"] &
  dat_big$nFeature_RNA < qc_thresholds$nFeature_RNA["upper"] &
  dat_big$percent.mt   > qc_thresholds$percent.mt["lower"]   &
  dat_big$percent.mt   < qc_thresholds$percent.mt["upper"]

dat_big_filtered <- dat_big[, keep_cells]

cat(sprintf("Cells before filtering: %d\n", ncol(dat_big)))
cat(sprintf("Cells after filtering:  %d\n", ncol(dat_big_filtered)))

# -- Post-filtering cell count barplot -----------------------------------------
cell_counts <- as.data.frame(table(dat_big_filtered$sample_name))
colnames(cell_counts) <- c("sample_name", "n_cells")
cell_counts$region <- ifelse(grepl("_DC_", cell_counts$sample_name), "DC", "PC")
cell_counts$sample_name <- factor(cell_counts$sample_name, levels = sample_order)

pdf(file.path(plot_dir, "cell_counts_postfiltering.pdf"), width = 16, height = 6)
ggplot(cell_counts, aes(x = sample_name, y = n_cells, fill = region)) +
  geom_col() +
  scale_fill_manual(values = region_colors) +
  geom_text(aes(label = n_cells), vjust = -0.3, size = 2.5) +
  geom_vline(xintercept = c(8.5, 16.5), linetype = "dashed", color = "grey50") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  ylab("Number of cells") +
  ggtitle("Cells per sample (post-filtering)")
dev.off()


# ==============================================================================
# 4. Normalization, RPCA integration, and clustering
# ==============================================================================

dat_big_filtered <- NormalizeData(dat_big_filtered)
dat_big_filtered <- FindVariableFeatures(dat_big_filtered,
                                         selection.method = "vst",
                                         nfeatures = 2000)
dat_big_filtered <- ScaleData(dat_big_filtered)
dat_big_filtered <- RunPCA(dat_big_filtered)

# -- RPCA integration ---------------------------------------------------------
options(future.globals.maxSize = 32 * 1024^3)

data <- IntegrateLayers(
  object         = dat_big_filtered,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  verbose        = FALSE
)

# -- Clustering and UMAP ------------------------------------------------------
data <- FindNeighbors(data, reduction = "integrated.rpca", dims = 1:30)
data <- FindClusters(data, resolution = 0.3, cluster.name = "rpca_clusters", random.seed = 42)
data <- RunUMAP(data, reduction = "integrated.rpca", dims = 1:30, random.seed = 42)
DimPlot(data, reduction = "umap", group.by = c("rpca_clusters", "condition"))

saveRDS(data, file.path(out_dir, "data_integrated_cellbender_df.rds"))

# ==============================================================================
# 5. Cluster marker detection
# ==============================================================================

data <- readRDS(file.path(out_dir, "data_integrated_cellbender_df.rds"))

# Proceed with JoinLayers
DefaultAssay(data) <- "RNA"
data <- JoinLayers(data)

# ==============================================================================
# 6. Cell type annotation
# ==============================================================================

# -- Manual annotation based on marker genes ----------------------------------
cluster_to_celltype <- c(
  "0"  = "B",
  "1"  = "Colonocytes",
  "2"  = "ISC_TA",
  "3"  = "Plasma",
  "4"  = "Trans_Goblet",
  "5"  = "Macrophages",
  "6"  = "T",
  "7"  = "T",
  "8"  = "Goblet",
  "9"  = "Fibroblasts", #"Interstitial",
  "10" = "SMC",
  "11" = "ISC_TA",
  "12" = "Fibroblasts", #"Telocytes",
  "13" = "Dendritic",
  "14" = "Lymphatic",
  "15" = "B_Proliferative",
  "16" = "Enteroendocrine",
  "17" = "Fibroblasts", #"Trophocytes",
  "18" = "Endothelial",
  "19" = "Goblet",
  "20" = "ICC",
  "21" = "Mesothelial",
  "22" = "Enteroendocrine",
  "23" = "T_Proliferative",
  "24" = "ILC3",
  "25" = "T_Proliferative",
  "26" = "Glia",
  "27" = "Dendritic",
  "28" = "Neurons",
  "29" = "Plasma",
  "30" = "ICC",
  "31" = "Tuft"
)

data@meta.data$celltype <- cluster_to_celltype[as.character(data$rpca_clusters)]
Idents(data) <- "celltype"


# ==============================================================================
# Subclustering ISC / TA
# ==============================================================================

data <- FindSubCluster(data, "ISC_TA", "RNA_snn", "CTA_cluster", resolution = 0.2)

TA_clusters <- grep("^ISC_TA_", unique(data$CTA_cluster), value = TRUE)

# -- Marker detection on subclusters only --------------------------------------
Markers <- wilcoxauc(
  data,
  group_by     = "CTA_cluster",
  groups_use   = TA_clusters,
  seurat_assay = "RNA"
)

# -- Score stemness and proliferation ------------------------------------------
stemness_genes <- c(
  "Lgr5", "Ascl2", "Olfm4", "Axin2", "Rnf43", "Ephb2", "Lrig1", "Sox9", "Smoc2"
)
proliferative_genes <- c(
  "Mki67", "Top2a", "Cdk1", "Ube2c", "Birc5", "Pcna", "Hmgb2", "Tubb5"
)

stemness_genes      <- intersect(stemness_genes, rownames(data))
proliferative_genes <- intersect(proliferative_genes, rownames(data))

data <- AddModuleScore(data, features = list(stemness_genes), name = "Stemness")
data <- AddModuleScore(data, features = list(proliferative_genes), name = "Proliferative")

FeatureScatter(subset(data, CTA_cluster %in% TA_clusters),
               feature1 = "Stemness1", feature2 = "Proliferative1",
               group.by = "CTA_cluster",split.by ='CTA_cluster') +
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") & theme(aspect.ratio=1, legend.position ='none') 

# -- Assign states (only for ISC_TA cells) ------------------------------------
is_iscta <- data$CTA_cluster %in% TA_clusters

data$state <- as.character(data$CTA_cluster) 

data$state[is_iscta &
             data$Stemness1 > 0.25 &
             data$Proliferative1 < 0.25] <- "ISC"

data$state[is_iscta &
             data$Stemness1 > 0.25 &
             data$Proliferative1 > 0.25] <- "TA/ISC_Proliferative"

data$state[is_iscta &
             data$Stemness1 < 0.25 &
             data$Proliferative1 > 0.25] <- "TA/ISC_Proliferative"

data$state[is_iscta &
             data$Stemness1 < 0.25 &
             data$Proliferative1 < 0.25] <- "Trans_Colonocytes"

# -- Remove contamination / low-quality subclusters ----------------------------
data <- subset(data, CTA_cluster != "ISC_TA_6" & CTA_cluster != "ISC_TA_8")

# -- Final celltype annotation ------------------------------------------------
# Combine: non-ISC_TA cells keep their original celltype,
# ISC_TA cells get the score-based state
data$celltype <- ifelse(data$CTA_cluster %in% TA_clusters,
                        data$state,
                        data$celltype)
Idents(data) ='celltype'
iscta_cells <- subset(data, CTA_cluster %in% TA_clusters)

state_colors <- c(
  "ISC"               = "darkgreen",
  "TA/ISC_Proliferative"  = "darkred",
  "Trans_Colonocytes" = "darkblue"
)


umap_a = DimPlot(iscta_cells, group.by = 'celltype',alpha=0.1, pt.size=0.05) +
  scale_color_manual(values=state_colors) + 
  theme(aspect.ratio=1, legend.position ='none') + 
  ggtitle('ISC/TA') +
  common_theme + 
  theme(plot.title  = element_text(face = "plain", hjust = 0.5)) +
  xlim(-10,6) + ylim(5,12)  

feature_b = FeaturePlot(iscta_cells,
                 features = c('Axin2',"Top2a",'Krt8','Stemness1','Proliferative1'),
                 reduction = "umap",
                 dims      = c(1, 2),
                 label     = FALSE,
                 raster    = TRUE,
                 cols      = c("gray100", "darkred"),
                 pt.size   = 4,
                 order     = TRUE,
                 alpha     = 0.5,
                 min.cutoff='q10',
                 max.cutoff='q95',
                 raster.dpi = c(1024, 1024),
                 ncol=3) & 
  xlim(-10,6) &
  ylim(5,12)   & 
  common_theme & 
  scale_color_gradient(low = "gray100", high = "darkred", breaks = scales::breaks_pretty(n = 4))
feature_b[[4]] <- feature_b[[4]] + ggtitle("Stemness score")
feature_b[[5]] <- feature_b[[5]] + ggtitle("Proliferative score") 


# ==============================================================================
# Subclustering Stromal/Fibroblasts
# ==============================================================================


data =FindSubCluster(data, "Fibroblasts","RNA_snn" ,"FBA_cluster",resolution=0.05)

fibro_clusters <- grep("^Fibroblasts_", unique(data$FBA_cluster), value = TRUE)
 
umap_c = DimPlot(subset(data,FBA_cluster %in% fibro_clusters), group.by = 'FBA_cluster',alpha=0.1, pt.size=0.05) +
  scale_color_manual(values=c('orange','darkgreen','darkred','black')) + 
  common_theme + 
  theme(plot.title  = element_text(face = "plain", hjust = 0.5),legend.position = 'None') +
  ggtitle("Stromal")  + xlim(c(-16,-8)) + ylim(c(-5,4))

feature_d = FeaturePlot(subset(data,FBA_cluster %in% fibro_clusters),
                 features = c("Pdgfra","Adamdec1","Wnt5a","Cd81","Grem1","Clu"),
                 reduction = "umap",
                 dims      = c(1, 2),
                 label     = FALSE,
                 raster    = TRUE,
                 cols      = c("gray100", "darkred"),
                 pt.size   = 4,
                 order     = TRUE,
                 alpha     = 0.5,
                 min.cutoff='q10',
                 max.cutoff='q90',
                 raster.dpi = c(1024, 1024),
                 ncol=3) & 
  common_theme & 
  xlim(c(-16,-8)) & 
  ylim(c(-5,4)) & 
  scale_color_gradient(low = "gray100", high = "darkred", breaks = scales::breaks_pretty(n = 4))


## give them new identity
data$FBA_cluster[data$FBA_cluster == 'Fibroblasts_0'] = 'Interstitial'
data$FBA_cluster[data$FBA_cluster == 'Fibroblasts_1'] = 'Telocytes'
data$FBA_cluster[data$FBA_cluster == 'Fibroblasts_2'] = 'Trophocytes'
data$FBA_cluster[data$FBA_cluster == 'Fibroblasts_3'] = 'Reticular'
Idents(data) = 'FBA_cluster'
data$celltype = as.character(Idents(data))

# ==============================================================================
# Subclustering Stromal/Fibroblasts
# ==============================================================================

data=FindSubCluster(data, "SMC","RNA_snn" ,"SMC_cluster",resolution=0.05)
muscle_clusters <- grep("^SMC_", unique(data$SMC_cluster), value = TRUE)
data$SMC_cluster[data$SMC_cluster == 'SMC_0'] = 'MP'
data$SMC_cluster[data$SMC_cluster == 'SMC_1'] = 'MM'
data$SMC_cluster[data$SMC_cluster == 'SMC_2'] = 'LPM'
data = subset(data,  SMC_cluster != 'SMC_3')
Idents(data) = 'SMC_cluster'
data$celltype = as.character(Idents(data))


umap_e = DimPlot(subset(data,SMC_cluster %in% c('LPM', 'MP', 'MM')), alpha=0.1, pt.size=0.05) +
  scale_color_manual(values=c('darkgreen','darkblue','darkred'))  + xlim(c(-15,-2)) + ylim(c(-2.5,5)) + 
  common_theme + 
  theme(plot.title  = element_text(face = "plain", hjust = 0.5),legend.position = 'None') +
  ggtitle("SMC")  



feature_f = FeaturePlot(subset(data,SMC_cluster %in% c('LPM', 'MP', 'MM')),
                features = c("Adamdec1","Hhip",'Grem2','Rspo3','Des','Actg2'),
                reduction = "umap",
                dims      = c(1, 2),
                label     = FALSE,
                raster    = TRUE,
                cols      = c("white", "darkred"),
                pt.size   = 4,
                order     = TRUE,
                alpha     = 0.5,
                raster.dpi = c(1024, 1024),
                min.cutoff='q10',
                max.cutoff='q99', ncol=3) &  
  xlim(c(-15,-2)) & 
  ylim(c(-2.5,5)) & 
  common_theme & 
  scale_color_gradient(low = "gray100", high = "darkred", breaks = scales::breaks_pretty(n = 4))

# -- Fig S4A/B---------------------------

# Build contours
contour_iscta <- make_contour(
  subset(data, celltype %in% c("ISC", "TA", "Trans_Colonocytes")),
  xlim = c(-10, 6), ylim = c(5, 12)
)

contour_fibro <- make_contour(
  subset(data, celltype %in% c("Interstitial", "Telocytes", "Trophocytes", "Reticular")),
  xlim = c(-16, -8), ylim = c(-5, 4)
)

contour_smc <- make_contour(
  subset(data, celltype %in% c("MP", "MM", "LPM")),
  xlim = c(-15, -2), ylim = c(-2.5, 5)
)

# Add to ISC/TA feature plots
for (i in seq_along(feature_b)) {
  feature_b[[i]] <- feature_b[[i]] + contour_iscta
}

# Add to fibroblast feature plots
for (i in seq_along(feature_d)) {
  feature_d[[i]] <- feature_d[[i]] + contour_fibro
}

# Add to SMC feature plots
for (i in seq_along(feature_f)) {
  feature_f[[i]] <- feature_f[[i]] + contour_smc
}

# Also add to UMAPs
umap_a <- umap_a + contour_iscta
umap_c <- umap_c + contour_fibro
umap_e <- umap_e + contour_smc

# Reassemble rows
p_row_1 <- (wrap_elements(umap_a) | wrap_elements(feature_b)) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(tag_levels = list(c("A", "B")))

p_row_2 <- (wrap_elements(umap_c) | wrap_elements(feature_d)) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(tag_levels = list(c("C", "D")))

p_row_3 <- (wrap_elements(umap_e) | wrap_elements(feature_f)) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(tag_levels = list(c("E", "F")))

p_full <- (p_row_1 / p_row_2 / p_row_3) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F")))

ggsave(file.path(plot_dir, "FigS4.pdf"),
       p_full, width = 210, height = 297, units = "mm", dpi = 300)

# -- Fig S3A/B/C/D---------------------------
celltype_to_broad <- c(
  # Immune
  "T"                    = "Immune",
  "B"                    = "Immune",
  "Plasma"               = "Immune",
  "B_Proliferative"      = "Immune",
  "T_Proliferative"      = "Immune",
  "Macrophages"          = "Immune",
  "Dendritic"            = "Immune",
  "ILC3"                 = "Immune",
  # Epithelial
  "Colonocytes"          = "Epithelial",
  "Trans_Colonocytes"    = "Epithelial",
  "Trans_Goblet"         = "Epithelial",
  "Goblet"               = "Epithelial",
  "ISC"                  = "Epithelial",
  "TA/ISC_Proliferative" = "Epithelial",
  "Enteroendocrine"      = "Epithelial",
  "Tuft"                 = "Epithelial",
  # SMC
  "MP"                   = "SMC",
  "MM"                   = "SMC",
  "LPM"                  = "SMC",
  # Fibroblast
  "Interstitial"         = "Fibroblast",
  "Telocytes"            = "Fibroblast",
  "Trophocytes"          = "Fibroblast",
  "Reticular"            = "Fibroblast",
  # Other
  "Endothelial"          = "Endothelial",
  "Lymphatic"            = "Lymphatic",
  "ICC"                  = "ICC",
  "Glia"                 = "Enteric_Nervous",
  "Neurons"              = "Enteric_Nervous",
  "Mesothelial"          = "Mesothelial"
)

# Verify before applying
missing <- setdiff(unique(data$celltype), names(celltype_to_broad))
if (length(missing) > 0) cat("Missing:", paste(missing, collapse = ", "), "\n")

data@meta.data$celltype.broad <- celltype_to_broad[data$celltype]

# Define the panel groupings explicitly
broad_groups <- list(
  "Immune"   = c("T", "B", "Plasma", "Macrophages", "Dendritic",
                 "ILC3", "T_Proliferative", "B_Proliferative"),
  "Epithelial" = c("ISC", "TA/ISC_Proliferative", "Trans_Colonocytes",
                   "Colonocytes", "Trans_Goblet", "Goblet",
                   "Enteroendocrine", "Tuft"),
  "Stromal & SMC" = c("Interstitial", "Telocytes", "Trophocytes",
                      "Reticular", "MP", "MM", "LPM"),
  "Others"   = c("Endothelial", "Lymphatic", "ICC",
                 "Glia", "Neurons", "Mesothelial")
)

panel_titles <- names(broad_groups)
panel_colors <- c("Immune" = "#E67E22", "Epithelial" = "#3498DB",
                  "Stromal & SMC" = "darkred", "Others" = "black")

plot_list <- list()

for (i in seq_along(broad_groups)) {
  ct <- intersect(broad_groups[[i]], unique(data$celltype))  # safety check
  cells_subset <- subset(data, celltype %in% ct)
  cells_subset$celltype <- factor(cells_subset$celltype, levels = ct)
  Idents(cells_subset) <- "celltype"
  
  Markers <- wilcoxauc(cells_subset, "celltype", seurat_assay = "RNA")
  Markers_df <- as.data.frame(
    top_markers(Markers, n = 5, auc_min = 0.5, padj_max = 0.01)
  )
  selected_markers <- unique(unlist(as.vector(Markers_df[, -1])))
  
  p <- DotPlot(cells_subset, features = selected_markers, dot.scale = 3,
               group.by = "celltype", scale = TRUE, cluster.idents = FALSE)
  
  # Extract the data from the plot and rebuild the point layer
  p <- p +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled),
               shape = 21, stroke = 0.3, color = "black") +
    scale_fill_gradient2(low = "steelblue", mid = "lightgrey",
                         high = "darkgoldenrod1",
                         limits = c(-2, 2), oob = scales::squish, name = "Average Expression") +
    guides(color = "none") +  # hide the original color legend
    RotatedAxis() +
    theme(axis.text.x  = element_text(size = 5, face = "italic"),
          axis.text.y  = element_text(size = 9),
          axis.title   = element_blank(),
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.3, "cm")) +
    ggtitle(panel_titles[i]) +
    theme(plot.title = element_text(color = panel_colors[panel_titles[i]],
                                    size = 13, face = "bold", hjust = 0.5))
  plot_list[[i]] <- p
}

#remove legend in the second and third
plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "none")
plot_list[[3]] <- plot_list[[3]] + theme(legend.position = "none")

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(plot_dir, "FigS3_dotplot_markers.pdf"),
       combined_plot, width = 297, height = 210, units = "mm", dpi = 300)


# ==============================================================================
# 7. Save final object
# ==============================================================================

saveRDS(data, file.path(out_dir, "data_integrated_celltype_cellbender_df.rds"))

cat("Done. Final object saved.\n")
