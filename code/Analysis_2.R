
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggstance)
library(enrichR)
library(ggtext)
library(dplyr)
library(reshape2)
library(parallel)
library(edgeR)
library(openxlsx)
library(ggrepel)
library(colorspace)

# -- Function definition --------------------------------------------------------------

f24_R2_cycling=function(x, t=2*(0:(length(x)-1)), period=24, group, single, offset=0)
{
  sel =  group %in% single
  if(!any(sel)){warning("Your single condition is not included in 'group'")}
  t  = t[sel]
  group  = group[sel]
  x= x[sel]
  n = length(x)
  #mu=mean(x)
  nb.timepoints=length(x)
  if(n<4)
  {
    if(n==0) c(nb.timepoints=nb.timepoints, mean=NA, amp=NA, relamp=NA,phase=NA,pval=NA)
    else
    {
      c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)
    }
  }
  else
  {
    sig2=var(x)
    c=cos(2*pi*t/period)
    s=sin(2*pi*t/period)
    A = mean(x*c)-mean(x)*mean(c)
    B = mean(x*s)-mean(x)*mean(s)
    c1 = mean(c^2)-mean(c)^2
    c2 = mean(c*s)-mean(c)*mean(s)
    c3 = mean(s^2)-mean(s)^2
    b = (A*c2-B*c1)/(c2^2-c1*c3)
    a = (A-b*c2)/c1
    mu = mean(x)-a*mean(c)-b*mean(s)
    #	b=2*mean(x*s)
    x.hat=mu+a*c+b*s
    sig2.1=var(x-x.hat)
    if(is.na(a)||is.na(b)) {c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)}
    else
    {
      p=3
      R2=0
      if(sig2>0) R2=1-sig2.1/sig2
      # http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
      # I checked that it works
      amp=max(x)-min(x)
      phase=period/(2*pi)*atan2(b, a)
      if(phase<0) phase=phase+period
      if(phase>period) phase=phase-period
      phase=(phase+offset)%%period
      pval = pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)
      c(nb.timepoints=nb.timepoints, mean =mean(x), amp=2*sqrt(a^2+b^2),relamp=sqrt(a^2+b^2)/(mu),phase=phase, pval=pval,tot_err=sum((x-x.hat)^2),a=a,b=b)
    }
  }
}
plot_data <- function(gene, subgroup = NULL, ss = 8) {
  
  DF <- NULL
  for (k in gene) {
    DF <- rbind(DF, data.frame(CN = as.numeric(fetched_data_bn[k, ]),
                               fetched_data_meta,
                               gene = k))
  }
  
  DF$ZT <- as.numeric(DF$ZT)
  DF <- DF[!DF$celltype %in% k_rem, ]
  
  DF <- rbind(DF, DF)
  DF[((nrow(DF) / 2) + 1):nrow(DF), "ZT"] <- DF[((nrow(DF) / 2) + 1):nrow(DF), "ZT"] + 24
  
  if (!is.null(subgroup)) {
    DF <- DF[DF$celltype %in% subgroup, ]
    DF$celltype <- factor(DF$celltype, levels = subgroup)
  }
  
  g1 <- ggplot(DF, aes(x = ZT, y = CN)) +
    geom_point(aes(col = cond), alpha = 0.7, size = 1) +
    stat_summary(aes(y = CN, group = cond, col = cond),
                 fun = mean, geom = "line", linewidth = 0.6, alpha = 0.7) +
    scale_color_manual(values = c("#EC6166", "#b8070d", "#6C9FCA")) +
    geom_vline(xintercept = 24, linetype = "dashed", color = "grey50") +
    theme_classic() +
    theme(aspect.ratio     = 1,
          axis.text.x      = element_text(size = ss),
          axis.text.y      = element_text(size = ss),
          strip.background = element_blank(),
          strip.text       = element_text(size = ss),
          plot.title       = element_text(face = "italic", hjust = 0.5),
          legend.position  = "none") +
    ylab("Log2 CPM") + xlab("ZT")
  
  if (length(gene) == 1) {
    g1 <- g1 +
      facet_wrap(~ celltype, scales = "free_y") +
      ggtitle(gene)
  } else {
    g1 <- g1 +
      facet_grid(gene ~ celltype, scales = "free_y") +
      theme(strip.text.y = element_text(face = "italic", size = ss))
  }
  
  return(g1)
}

# ==============================================================================
# Complex SVD function
# ==============================================================================
run_csvd <- function(genes, DS_list, k_rem, color_variations,
                     pval_threshold = 0.2, mode = 1,l.pos='none',g.s=3) {
  
  # Build complex matrix from cosinor coefficients
  DSS <- t(do.call(rbind, lapply(DS_list, function(x) {
    complex(real = x[genes, "a"], imaginary = x[genes, "b"])
  })))
  rownames(DSS) <- genes
  
  # Set non-significant or missing to zero
  DSS_pv <- t(do.call(rbind, lapply(DS_list, function(x) x[genes, "pval"])))
  DSS[is.na(DSS)] <- 0
  DSS[is.na(DSS_pv) | DSS_pv > pval_threshold] <- 0

  colnames(DSS) <- names(DS_list)
  # Remove low-cell-count celltypes
  keep <- !gsub(".+,", "", colnames(DSS)) %in% k_rem
  DSS <- DSS[, keep]
  genes=rownames(DSS)
  # SVD
  SVD <- svd(DSS)
  
  # Rotate modes for interpretability
  for (k in seq_along(SVD$d)) {
    mn  <- sum(SVD$v[, k])
    rot <- Conj(mn) / Mod(mn)
    SVD$u[, k] <- SVD$u[, k] * rot * max(Mod(SVD$v[, k])) * SVD$d[k]
    SVD$v[, k] <- Conj(SVD$v[, k] * rot / max(Mod(SVD$v[, k])))
  }
  
  # Extract mode i
  i <- mode
  varexp <- round((SVD$d[i]^2) / sum(SVD$d^2) * 100, 0)
  
  # Gene space (U)
  SVD.u <- data.frame(
    labs = genes,
    x    = SVD$u[, i]
  )
  
  p_gene <- ggplot(SVD.u, aes(x = Arg(x) %% (2 * pi) * 12 / pi, y = Mod(x))) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = labs), size = 4, max.overlaps = 50,
                    fontface = "italic") +
    coord_polar() +
    scale_x_continuous(breaks = seq(0, 24, by = 4), expand = c(0, 0),
                       limits = c(0, 24)) +
    ylim(c(0, NA)) +
    theme_minimal() +
    theme(text = element_text(size = 12)) +
    labs(title = paste0("Gene space"),
         x = "ZT", y = expression(Log[2]~fold-change))
  
  # Tissue space (V) colored by celltype
  SVD.v <- data.frame(
    cell_type = gsub(".+,", "", colnames(DSS)),
    cond      = gsub(",.+", "", colnames(DSS)),
    x         = SVD$v[, i]
  )
  
  p_tissue <- ggplot(SVD.v, aes(x = Arg(x) %% (2 * pi) * 12 / pi,
                                y = Mod(x))) +
    geom_point(aes(shape = cond, col = cell_type), size = 2) +
    geom_text_repel(aes(label = cell_type, col = cell_type), size = g.s,
                    max.overlaps = 10) +
    coord_polar() +
    scale_x_continuous(breaks = seq(0, 24, by = 4), expand = c(0, 0),
                       limits = c(0, 24)) +
    scale_color_manual(values = color_variations) +
    ylim(c(0, NA)) +
    theme_minimal() +
    theme(text = element_text(size = 12), legend.position=l.pos) +
    labs(title = paste0("Tissue space"),
         x = "ZT", y = "Scaling factor ")
  
  cond_colors <- c("Control_DC" = "#EC6166",
                   "Regenerating_DC" = "#b8070d",
                   "Regenerating_PC" = "#6C9FCA")
  
  p_tissue_cond <- ggplot(SVD.v, aes(x = Arg(x) %% (2 * pi) * 12 / pi,
                                     y = Mod(x))) +
    geom_point(aes(shape = cond, col = cond), size = 2) +
    geom_text_repel(aes(label = cell_type, col = cond), size = 3,
                    max.overlaps = 10) +
    coord_polar() +
    scale_x_continuous(breaks = seq(0, 24, by = 4), expand = c(0, 0),
                       limits = c(0, 24)) +
    scale_color_manual(values = cond_colors) +
    scale_shape_manual(values = c(16, 17, 15)) +
    ylim(c(0, NA)) +
    theme_minimal() +
    theme(text = element_text(size = 12), legend.position= l.pos) +
    labs(title = "Tissue space",
         x = "ZT", y = "Scaling factor")
  
  list(gene = p_gene, tissue = p_tissue, tissue_cond = p_tissue_cond,
       svd = SVD, varexp = varexp)
  }


#############################
# Resolve repo root: REPO_DIR env var > here::here() > getwd().
# Run this script from the repo root, or set REPO_DIR explicitly.
base_dir <- Sys.getenv("REPO_DIR", unset = NA)
if (is.na(base_dir)) {
  base_dir <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
}
out_dir  <- file.path(base_dir, "data")
plot_dir <- file.path(base_dir, "plot", "ctrl")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

data <- readRDS(file.path(out_dir, "data_integrated_celltype_cellbender_df.rds"))

data@meta.data$condition[data@meta.data$condition=='DSS'] = 'Regenerating'
data@meta.data$condition[data@meta.data$condition=='CTRL'] = 'Control'
data@meta.data$cond = paste(data@meta.data$condition, data@meta.data$region, sep="_")

data$ZT <- as.character(data$ZT)
data$ZT[data$ZT == "1"]  <- "01"
data$ZT[data$ZT == "7"]  <- "07"

# -- Colors (assuming col_vec and color_variations are loaded) -----------------

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


col_vec <- c("#D55E00", "grey", "#56B4E9", "#009E73", "#0072B2",
             "#E69F00", "#CC79A7", "#000000", "#9E0142")
names(col_vec) <- sort(unique(data$celltype.broad))

# Build a per-celltype color palette by generating light/base/dark variations
# of each broad-category color, one shade per cell type within the category.
make_color_variations <- function(celltype_to_broad, base_palette) {
  variations_for <- function(base_color, n) {
    if (n == 1) return(base_color)
    if (n == 2) return(c(base_color, lighten(base_color, amount = 0.6)))
    shades <- c(lighten(base_color, amount = 0.6),
                base_color,
                darken(base_color, amount = 0.6))
    if (n > 3) shades <- colorRampPalette(shades)(n)
    shades
  }

  broad_levels <- sort(unique(celltype_to_broad))
  cv <- unlist(setNames(
    lapply(seq_along(broad_levels), function(i) {
      variations_for(base_palette[i], sum(celltype_to_broad == broad_levels[i]))
    }),
    broad_levels
  ))
  cv <- cv[order(names(cv))]
  names(cv) <- names(celltype_to_broad)[order(celltype_to_broad)]
  cv
}

color_variations <- make_color_variations(celltype_to_broad, col_vec)

# -- Figure 3 scheme: A (scheme) + B (scheme) - leave empty, done in Illustrator -----------------

# -- Without labels (high-res PNG for background) -----------------------------
p_c_nolabel <- DimPlot(data, reduction = "umap", group.by = "celltype.broad",
                       label = FALSE, cols = col_vec,
                       raster = FALSE, alpha=0.2) +
  NoLegend() + NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

p_d_nolabel <- DimPlot(data, reduction = "umap", group.by = "celltype",
                       label = FALSE, cols = color_variations,
                       raster = FALSE, alpha=0.2) +
  NoLegend() + NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

cols <- c(
  "01"  = "#FFD700",  
  "07"  = "#FF8C00",  
  "13" = "#6A5ACD",  
  "19" = "#1B1B3A"   
)
p_e_zt_nolabel <- DimPlot(data, reduction = "umap", group.by = "ZT",
                          alpha = 0.2,
                          raster = FALSE, cols=cols) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

p_e_cond_nolabel <- DimPlot(data, reduction = "umap", group.by = "cond",
                            cols = c("#EC6166", "#b8070d", "#6C9FCA"),
                            alpha = 0.2,
                            raster = FALSE) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())



ggsave(file.path(plot_dir, "Fig1_C_nolabel.png"),
       p_c_nolabel, width = 12, height = 12,  dpi = 400)
ggsave(file.path(plot_dir, "Fig1_D_nolabel.png"),
       p_d_nolabel, width = 12, height = 12, dpi = 400)
ggsave(file.path(plot_dir, "Fig1_E_zt_nolabel.png"),
       p_e_zt_nolabel, width = 12, height = 12,  dpi = 400)
ggsave(file.path(plot_dir, "Fig1_E_cond_nolabel.png"),
       p_e_cond_nolabel, width = 12, height = 12,  dpi = 400)

# -- With labels (PDF vector for Illustrator) ----------------------------------
p_c_label <- DimPlot(data, reduction = "umap", group.by = "celltype.broad",
                     label = TRUE, repel = TRUE, label.size = 4,
                     cols = col_vec, pt.size = 0.1,
                     raster = TRUE, raster.dpi = c(512, 512)) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

p_d_label <- DimPlot(data, reduction = "umap", group.by = "celltype",
                     label = TRUE, repel = TRUE, label.size = 3,
                     cols = color_variations, pt.size = 0.1,
                     raster = TRUE, raster.dpi = c(512, 512)) +
   NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

ggsave(file.path(plot_dir, "Fig1_C_labels.pdf"),
       p_c_label, width = 12, height = 12)
ggsave(file.path(plot_dir, "Fig1_D_labels.pdf"),
       p_d_label, width = 12, height = 12)

# -- Panel F: Stacked barplot proportions --------------------------------------

pt <- as.data.frame(table(data$celltype, paste(data$ZT, data$sex), data$cond))
colnames(pt) <- c("celltype", "sample", "cond", "Freq")
pt$celltype <- factor(pt$celltype, levels = names(celltype_to_broad)[order(celltype_to_broad)])

pt_fraction <- pt %>%
  group_by(sample, cond) %>%
  mutate(Fraction = Freq / sum(Freq)) %>%
  ungroup()

# Split by region for DC and PC facets
p_f <- ggplot(pt_fraction, aes(x = sample, y = Fraction, fill = celltype)) +
  geom_col(width = 0.7) +
  facet_grid(. ~ cond, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_variations) +
  theme_classic() +
  theme(axis.text.x     = element_text(angle = 90, hjust = 1, size = 6),
        legend.position  = "none",
        panel.grid       = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 10)) +
  ylab("Fraction") + xlab("")

ggsave(file.path(plot_dir, "Fig1_F_fraction.pdf"),
       p_f, width = 10, height = 5)


# ==============================================================================
# S5.1 Pseudobulk aggregation and normalization
# ==============================================================================

DefaultAssay(data) <- "RNA"

# Filter lowly expressed genes
data.2 <- data[rowSums(data) > 1000, ]

# Create grouping variable
data.2$combined <- paste(data.2$sex, data.2$cond, data.2$celltype, data.2$ZT,
                         sep = ",")

# Aggregate to pseudobulk
pseudo_bulk <- AggregateExpression(
  data.2,
  group.by = "combined",
  assays   = "RNA",
  return.seurat = FALSE
)

fetched_data <- as.data.frame(pseudo_bulk$RNA)
colnames(fetched_data) <- gsub("-", "_", colnames(fetched_data))

# Parse metadata from column names
fetched_data_meta <- as.data.frame(
  do.call(rbind, strsplit(colnames(fetched_data), split = ","))
)
colnames(fetched_data_meta) <- c("sex", "cond", "celltype", "ZT")

# TMM normalization + log2 CPM
dge <- DGEList(counts = fetched_data)
dge <- calcNormFactors(dge, method = "TMM")
fetched_data_bn <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)


# ==============================================================================
# S5.2 Harmonic regression (24h cosinor) per celltype x condition
# ==============================================================================

group <- paste0(fetched_data_meta$cond, ",", fetched_data_meta$celltype)
group.l <- as.list(unique(group))
time <- as.numeric(fetched_data_meta$ZT)

n_cores <- min(12, max(1, parallel::detectCores() - 1))
DS_bn <- mclapply(group.l, function(x) {
  as.data.frame(t(apply(fetched_data_bn, 1, f24_R2_cycling,
                        t = time, period = 24,
                        group = group, single = x)))
}, mc.cores = n_cores)

# BH correction per group
for (i in seq_along(DS_bn)) {
  DS_bn[[i]]$qval <- p.adjust(DS_bn[[i]]$pval, method = "BH")
}

# Set genes with insufficient coverage to NA
for (i in seq_along(group.l)) {
  dat.sub <- fetched_data[, grep(group.l[i], colnames(fetched_data), fixed = TRUE)]
  rem <- names(which(rowSums(dat.sub) < 50 |
                       apply(dat.sub, 1, function(x) sum(x != 0)) < 5))
  DS_bn[[i]][rem, "qval"] <- NA
  DS_bn[[i]][rem, "pval"] <- NA
}

# ==============================================================================
# S5.3 Export supplementary table (Table S1)
# ==============================================================================

#--- Cell count filtering

# Level 1: Remove globally rare celltypes
nb_cell <- table(data$celltype)
k_rem <- names(nb_cell[nb_cell < 800])
cat("Removed celltypes (< 800 cells total):", paste(k_rem, collapse = ", "), "\n")

# Level 2: Check per-group coverage (celltype x condition x ZT)
cells_per_group <- as.data.frame(
  table(data.2$celltype, data.2$cond, data.2$ZT)
)
colnames(cells_per_group) <- c("celltype", "cond", "ZT", "n_cells")

# Flag celltype x condition combos where any ZT has < 20 cells
low_coverage <- cells_per_group %>%
  filter(!celltype %in% k_rem) %>%
  group_by(celltype, cond) %>%
  summarize(min_cells = min(n_cells),
            n_timepoints = sum(n_cells > 0),
            .groups = "drop") %>%
  filter(min_cells < 20 | n_timepoints < 4)

if (nrow(low_coverage) > 0) {
  cat("Warning: low coverage groups:\n")
  print(low_coverage)
}

k_rem = c(k_rem, as.character(low_coverage$celltype))
uniq <- unlist(group.l)

wb <- createWorkbook()

for (i in seq_along(uniq)) {
  ct <- unlist(strsplit(uniq[i], split = ","))[2]
  if (!ct %in% k_rem) {
    nam <- gsub("/", "_", uniq[i])
    nam <- substr(nam, 1, 31)
    df <- subset(DS_bn[[i]], pval < 0.05 & amp > 0.5)
    df <- df[order(df$pval), ]
    df <- df[is.finite(df$relamp), ]
    df$gene <- rownames(df)
    addWorksheet(wb, nam)
    writeData(wb, nam, df)
  }
}

all_rhythmic_genes <- c()
for (i in seq_along(uniq)) {
  ct <- unlist(strsplit(uniq[i], split = ","))[2]
  if (!ct %in% k_rem) {
    df <- subset(DS_bn[[i]], pval < 0.05 & amp > 0.5)
    df <- df[is.finite(df$relamp), ]
    all_rhythmic_genes <- c(all_rhythmic_genes, rownames(df))
  }
}
cat("Total unique rhythmic genes:", length(unique(all_rhythmic_genes)), "\n")


saveWorkbook(wb, file.path(out_dir, "Table_S1.xlsx"), overwrite = TRUE)
cat("Table S1 saved.\n")

# ==============================================================================
# Figure S5: Rhythmic gene counts vs thresholds
# ==============================================================================


# -- Condition labels with color coding ----------------------------------------

cond_title <- c(
  "Control_DC"      = "<span style='color:#EC6166'>Control Distal Colon</span>",
  "Regenerating_DC" = "<span style='color:#b8070d'>Regenerating Distal Colon</span>",
  "Regenerating_PC" = "<span style='color:#6C9FCA'>Regenerating Proximal Colon</span>"
)

# -- S5.1 Number of rhythmic genes vs amplitude threshold ----------------------

DS_all_amp <- NULL
for (g in seq(0, 3, 0.1)) {
  DS.r.tmp <- lapply(DS_bn, function(x) subset(x, pval < 0.05 & amp > g))
  DS_all_amp <- rbind(DS_all_amp, unlist(lapply(DS.r.tmp, nrow)))
}
colnames(DS_all_amp) <- group.l
rownames(DS_all_amp) <- seq(0, 3, 0.1)
DS_all_amp <- melt(DS_all_amp)
DS_all_amp$condition <- gsub(",.+", "", DS_all_amp$Var2)
DS_all_amp$cell_type <- gsub(".+,", "", DS_all_amp$Var2)
DS_all_amp <- DS_all_amp[!DS_all_amp$cell_type %in% k_rem, ]

# -- S5.2 Number of rhythmic genes vs p-value threshold ------------------------

pval_seq <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)

DS_all_pv <- NULL
for (p in pval_seq) {
  DS.r.tmp <- lapply(DS_bn, function(x) subset(x, pval < p & amp > 0.5))
  DS_all_pv <- rbind(DS_all_pv, unlist(lapply(DS.r.tmp, nrow)))
}
colnames(DS_all_pv) <- group.l
rownames(DS_all_pv) <- pval_seq
DS_all_pv <- melt(DS_all_pv)
DS_all_pv$condition <- gsub(",.+", "", DS_all_pv$Var2)
DS_all_pv$cell_type <- gsub(".+,", "", DS_all_pv$Var2)
DS_all_pv <- DS_all_pv[!DS_all_pv$cell_type %in% k_rem, ]

# -- Shared theme and y-axis ---------------------------------------------------

common_nb_theme <- theme_bw() +
  theme(aspect.ratio     = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "none")

y_max <- max(max(DS_all_amp$value, na.rm = TRUE),
             max(DS_all_pv$value, na.rm = TRUE))

# -- Plot functions ------------------------------------------------------------

make_amp_plot <- function(cond_name) {
  ggplot(subset(DS_all_amp, condition == cond_name),
         aes(x = Var1, y = value)) +
    geom_line(aes(col = cell_type)) +
    xlab(expression(Log[2] ~ amplitude)) + ylab("Nb. rhythmic genes") +
    ggtitle(cond_title[cond_name]) +
    scale_color_manual(values = color_variations) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey",
               linewidth = 0.5) +
    common_nb_theme + ylim(0, y_max_amp) +
    theme(plot.title = element_markdown(size = 11, hjust = 0.5))
}

make_pv_plot <- function(cond_name) {
  ggplot(subset(DS_all_pv, condition == cond_name),
         aes(x = -log10(Var1), y = log2(1+value))) +
    geom_line(aes(col = cell_type)) +
    xlab(expression(-log[10] ~ p - value)) + ylab("Nb. rhythmic genes (log2)") +
    ggtitle(cond_title[cond_name]) +
    scale_color_manual(values = color_variations) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey",
               linewidth = 0.5) +
    common_nb_theme + ylim(0, log2(1+y_max_pv)) +
    theme(plot.title = element_markdown(size = 11, hjust = 0.5))
}

# -- Build panels --------------------------------------------------------------

S5A <- make_amp_plot("Control_DC")
S5B <- make_amp_plot("Regenerating_DC")
S5C <- make_amp_plot("Regenerating_PC")
S5D <- make_pv_plot("Control_DC")
S5E <- make_pv_plot("Regenerating_DC")
S5F <- make_pv_plot("Regenerating_PC")

# -- Assemble ------------------------------------------------------------------

p_S5 <- (S5A | S5B | S5C) /
  (S5D | S5E | S5F) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(plot_dir, "FigS5_rhythmic_genes.pdf"),
       p_S5, width = 297, height = 200, units = "mm")

# ==============================================================================
# Figure 4: Rhythmic gene landscape
# ==============================================================================

# -- Rhythmic gene counts -----------------------------------------------------

DS.r <- lapply(DS_bn, function(x) subset(x, pval < 0.05 & amp > 0.5))

nb.rhythmic <- data.frame(
  nb   = unlist(lapply(DS.r, nrow)),
  name = unlist(group.l)
)
nb.rhythmic$cond    <- gsub(",.+", "", nb.rhythmic$name)
nb.rhythmic$cluster <- gsub(".+,", "", nb.rhythmic$name)
nb.rhythmic <- nb.rhythmic[!nb.rhythmic$cluster %in% k_rem, ]

ctrl_nb <- nb.rhythmic[nb.rhythmic$cond == "Control_DC", ]
ctrl_nb$broad <- celltype_to_broad[as.character(ctrl_nb$cluster)]

# Define broad category order
broad_order <- c("Epithelial", "Fibroblast", "SMC", "Immune",
                 "Endothelial", "Lymphatic", "ICC", "Enteric_Nervous", "Mesothelial")

ctrl_nb$broad <- factor(ctrl_nb$broad, levels = rev(broad_order))
ctrl_nb <- ctrl_nb[order(ctrl_nb$broad, ctrl_nb$nb), ]
ct_order <- as.character(ctrl_nb$cluster)

# Apply fixed order to all conditions
nb.rhythmic$cluster <- factor(nb.rhythmic$cluster, levels = ct_order)
nb.rhythmic$broad <- celltype_to_broad[as.character(nb.rhythmic$cluster)]
nb.rhythmic$broad <- factor(nb.rhythmic$broad, levels = rev(broad_order))

# -- Phase data ----------------------------------------------------------------

DS.r.phase <- lapply(DS.r, function(x) x[, "phase"])
names(DS.r.phase) <- group.l
DS.r.phase <- stack(DS.r.phase)
DS.r.phase$cell_type <- gsub(".+,", "", DS.r.phase$ind)
DS.r.phase$cond      <- gsub(",.+", "", DS.r.phase$ind)
DS.r.phase$phase     <- 2 * pi * DS.r.phase$values / 24
DS.r.phase$broad     <- celltype_to_broad[DS.r.phase$cell_type]
DS.r.phase <- DS.r.phase[!DS.r.phase$cell_type %in% k_rem, ]

nbins    <- 48
binwidth <- (2 * pi) / nbins

# -- Helper functions ----------------------------------------------------------

y_max_bar <- max(nb.rhythmic$nb, na.rm = TRUE)

make_barplot <- function(cond_name) {
  df <- nb.rhythmic[nb.rhythmic$cond == cond_name, ]
  
  ggplot(df, aes(x = cluster, y = nb, fill = cluster)) +
    geom_col() +
    scale_fill_manual(values = color_variations) +
    coord_flip(ylim = c(0, y_max_bar)) +
    theme_classic() +
    theme(axis.text.x      = element_text(size = 8),
          axis.text.y      = element_text(size = 8),
          legend.position  = "none",
          plot.title       = element_text(hjust = 0.5, size = 11)) +
    xlab("") + ylab("Nb. rhythmic genes") +
    ggtitle(gsub("_", " ", cond_name))
}

make_phase_plot <- function(cond_name) {
  ggplot(subset(DS.r.phase, cond == cond_name),
         aes(x = phase,
             y = after_stat(count / tapply(..count.., ..PANEL.., sum)[..PANEL..]),
             fill = cell_type, colour = cell_type)) +
    geom_histogram(binwidth = binwidth, position = "identity",
                   alpha = 0.6, linewidth = 0.1) +
    facet_wrap(~ broad, ncol = 4, scales = "free_y") +
    scale_x_continuous(limits = c(0, 2 * pi),
                       breaks = seq(0, 2 * pi, by = pi / 2)[-5],
                       labels = c("0", "6", "12", "18")) +
    scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    scale_fill_manual(values = color_variations) +
    scale_color_manual(values = color_variations) +
    theme_bw() +
    theme(aspect.ratio     = 1,
          axis.text        = element_text(size = 10),
          axis.text.y      = element_blank(),
          axis.ticks       = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey60"),
          panel.border     = element_blank(),
          strip.background = element_blank(),
          strip.text       = element_text(size = 10),
          legend.position  = "none") +
    xlab("ZT") + ylab("")
}



# ==============================================================================
# Figure 4D/E/F: KEGG dotplot by broad category, one per condition
# ==============================================================================

# -- Enrichment on rhythmic genes per broad category x condition ---------------

DS.r <- lapply(DS_bn, function(x) subset(x, pval < 0.05 & amp > 0.5))
names(DS.r) <- unlist(group.l)

ds_info <- data.frame(
  group    = names(DS.r),
  cond     = gsub(",.+", "", names(DS.r)),
  celltype = gsub(".+,", "", names(DS.r))
)
ds_info <- ds_info[!ds_info$celltype %in% k_rem, ]
ds_info$broad <- celltype_to_broad[ds_info$celltype]

broad_cond_groups <- split(ds_info$group, paste(ds_info$cond, ds_info$broad, sep = "|"))

DS.r_broad <- lapply(broad_cond_groups, function(groups) {
  unique(unlist(lapply(DS.r[groups], rownames)))
})

dbb <- c("KEGG_2019_Mouse")

GOL_broad <- list()
for (i in seq_along(DS.r_broad)) {
  if (length(DS.r_broad[[i]]) < 10) next
  GOL_broad[[names(DS.r_broad)[i]]] <- enrichr(DS.r_broad[[i]], databases = dbb)
}

# -- Collect all results -------------------------------------------------------

ALL_broad <- NULL
for (condi in names(GOL_broad)) {
  for (f in names(GOL_broad[[condi]])) {
    g.kegg <- GOL_broad[[condi]][[f]]
    if (nrow(g.kegg) == 0) next
    g.kegg <- g.kegg[, c("Term", "Adjusted.P.value", "Genes", "Overlap", "Odds.Ratio")]
    g.kegg$tot        <- as.numeric(sapply(strsplit(g.kegg$Overlap, "/"), "[[", 2))
    g.kegg$nb_gene    <- as.numeric(sapply(strsplit(g.kegg$Overlap, "/"), "[[", 1))
    g.kegg$odds_ratio <- g.kegg$Odds.Ratio
    g.kegg$type  <- f
    g.kegg$group <- condi
    colnames(g.kegg)[2] <- "adj.Pval"
    ALL_broad <- rbind(ALL_broad, g.kegg)
  }
}

ALL_broad$condition  <- gsub("\\|.+", "", ALL_broad$group)
ALL_broad$broad      <- gsub(".+\\|", "", ALL_broad$group)
ALL_broad$score      <- pmin(-log10(ALL_broad$adj.Pval), 5)
ALL_broad$odds_ratio <- pmin(ALL_broad$odds_ratio, 20)

# -- Filter --------------------------------------------------------------------

disease_pattern <- paste(
  "disease", "cancer", "Hepatitis", "Influenza", "leukemia", "infection",
  "carcinogenesis", "carcinoma", "diabetic", "Legionellosis", "syndrome",
  "Melanoma", "Glioma", "Allograft rejection", "Amoebiasis", "Toxoplasmosis",
  "Amphetamine addiction", "Oocyte meiosis", "Fanconi anemia",
  sep = "|"
)

ALL_broad <- ALL_broad[!grepl(disease_pattern, ALL_broad$Term, ignore.case = TRUE), ]
ALL_broad <- subset(ALL_broad, tot < 750 & nb_gene > 3)

# -- Broad category ordering and colors ----------------------------------------

broad_order <- c("Epithelial", "Immune", "Fibroblast", "SMC",
                 "Endothelial", "Lymphatic", "ICC", "Enteric_Nervous", "Mesothelial")

broad_colors <- c(
  "Epithelial"      = "#56B4E9",
  "Immune"          = "#E69F00",
  "Fibroblast"      = "#009E73",
  "SMC"             = "#D55E00",
  "Endothelial"     = "#CC79A7",
  "Lymphatic"       = "#0072B2",
  "ICC"             = "#000000",
  "Enteric_Nervous" = "#9E0142",
  "Mesothelial"     = "grey"
)

# -- Condition titles ----------------------------------------------------------

cond_titles <- c(
  "Control_DC"      = "<span style='color:#EC6166'>Control Distal Colon</span>",
  "Regenerating_DC" = "<span style='color:#b8070d'>Regenerating Distal Colon</span>",
  "Regenerating_PC" = "<span style='color:#6C9FCA'>Regenerating Proximal Colon</span>"
)

# -- Helper: dotplot for one condition -----------------------------------------

make_kegg_dotplot <- function(cond_chosen, odds_range, odds_breaks) {
  
  df_cond <- subset(ALL_broad, condition == cond_chosen)
  
  # Terms significant in THIS condition only
  sig_terms <- unique(df_cond$Term[df_cond$adj.Pval < 0.1])
  if (length(sig_terms) == 0) return(NULL)
  
  # Keep all broad categories for those terms (show non-significant as faint)
  df_dot <- df_cond[df_cond$Term %in% sig_terms, ]
  
  # Hierarchical clustering by gene overlap
  all_genes <- unique(unlist(strsplit(df_dot$Genes, ";")))
  term_gene_mat <- matrix(0, nrow = length(sig_terms), ncol = length(all_genes),
                          dimnames = list(sig_terms, all_genes))
  
  for (i in seq_len(nrow(df_dot))) {
    genes <- unlist(strsplit(df_dot$Genes[i], ";"))
    term_gene_mat[df_dot$Term[i], intersect(genes, all_genes)] <- 1
  }
  
  term_gene_mat <- term_gene_mat[!duplicated(term_gene_mat), , drop = FALSE]
  sig_terms_clean <- rownames(term_gene_mat)
  
  if (length(sig_terms_clean) < 3) {
    term_order <- sig_terms_clean
    term_representatives <- sig_terms_clean
  } else {
    term_dist <- dist(term_gene_mat, method = "binary")
    term_clust <- hclust(term_dist, method = "ward.D2")
    
    # Cut and keep representative per cluster
    n_clusters <- max(2, round(length(sig_terms_clean) / 3))
    n_clusters <- min(n_clusters, length(sig_terms_clean))
    term_groups <- cutree(term_clust, k = n_clusters)
    
    term_representatives <- df_dot %>%
      filter(Term %in% sig_terms_clean) %>%
      mutate(term_group = term_groups[as.character(Term)]) %>%
      group_by(term_group) %>%
      slice_min(adj.Pval, n = 1, with_ties = FALSE) %>%
      pull(Term) %>%
      unique()
    
    term_order <- sig_terms_clean[term_clust$order]
    term_order <- intersect(term_order, term_representatives)
  }
  
  df_dot <- df_dot[df_dot$Term %in% term_representatives, ]
  df_dot$Term <- factor(df_dot$Term, levels = term_order)
  df_dot$Term <- droplevels(df_dot$Term)
  
  # Broad ordering
  df_dot$broad <- factor(df_dot$broad,
                         levels = broad_order[broad_order %in% df_dot$broad])
  
  broad_labels <- setNames(
    paste0("<span style='color:", broad_colors[levels(df_dot$broad)], "'>**",
           levels(df_dot$broad), "**</span>"),
    levels(df_dot$broad)
  )
  
  cat(cond_chosen, ":", length(term_representatives), "representative terms\n")
  
  ggplot(df_dot, aes(x = broad, y = Term)) +
    geom_point(aes(size = odds_ratio, fill = score),
               shape = 21, stroke = 0.3, color = "black") +
    scale_fill_gradient(low = "white", high = "darkred",
                        limits = c(0, 5),
                        name = expression(-log[10] ~ p[adj])) +
    scale_size_continuous(range = c(1, 6), name = "Odds ratio",
                          limits = odds_range,
                          breaks = odds_breaks) +
    scale_x_discrete(labels = broad_labels) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x      = element_markdown(angle = 45, hjust = 1, size = 9),
          axis.text.y      = element_text(size = 8),
          axis.title       = element_blank(),
          plot.title       = element_markdown(size = 12, hjust = 0.5),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          legend.position  = "right") +
    ggtitle(cond_titles[cond_chosen])
}

# -- Build plots ---------------------------------------------------------------
odds_range <- range(ALL_broad$odds_ratio[ALL_broad$adj.Pval < 0.1], na.rm = TRUE)
odds_breaks <- c(2, 5, 10, 20)
odds_breaks <- odds_breaks[odds_breaks <= max(odds_range)]

p_kegg_ctrl     <- make_kegg_dotplot("Control_DC", odds_range, odds_breaks)
p_kegg_regen_dc <- make_kegg_dotplot("Regenerating_DC", odds_range, odds_breaks)
p_kegg_regen_pc <- make_kegg_dotplot("Regenerating_PC", odds_range, odds_breaks)


# ==============================================================================
# Assemble full Figure 4
# ==============================================================================

# Row 1: Control DC
row_1 <- wrap_elements(make_barplot("Control_DC")) |
  wrap_elements(make_phase_plot("Control_DC")) |
  wrap_elements(p_kegg_ctrl)

# Row 2: Regenerating DC
row_2 <- wrap_elements(make_barplot("Regenerating_DC")) |
  wrap_elements(make_phase_plot("Regenerating_DC")) |
  wrap_elements(p_kegg_regen_dc)

# Row 3: Regenerating PC
row_3 <- wrap_elements(make_barplot("Regenerating_PC")) |
  wrap_elements(make_phase_plot("Regenerating_PC")) |
  wrap_elements(p_kegg_regen_pc)

p_fig4 <- (row_1 / row_2 / row_3) +
  plot_layout(widths = c(1, 2, 1.5)) +
  plot_annotation(tag_levels = list(c("A", "B", "C",
                                      "D", "E", "F",
                                      "G", "H", "I")))

ggsave(file.path(plot_dir, "Fig4_rhythmic_landscape.pdf"),
       p_fig4, width = 500, height = 350, units = "mm", dpi = 300)

cat("Figure 4 saved.\n")
# ==============================================================================
# Assemble full Figure 6
# ==============================================================================

# -- Fig 6A: Core clock gene examples -----------------------------------------

main_celltypes <- c("Goblet", "Colonocytes", "MM", "MP", "Telocytes",
                    "Trophocytes", "B", "T", "Macrophages", "ICC", "ISC",
                    "TA/ISC_Proliferative", "Trans_Colonocytes", "Trans_Goblet")

main_celltypes <- c("MM", "MP", "Telocytes", "Trophocytes",  "ISC",
                    "TA/ISC_Proliferative", "Trans_Colonocytes", "Trans_Goblet","Goblet", "Colonocytes")


p_per2 <- plot_data("Arntl", main_celltypes, 9)
p_dbp  <- plot_data("Dbp",  main_celltypes, 9)

# -- Fig 6B/C: cSVD on core clock genes ---------------------------------------

clocks <- c("Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3",
            "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Ciart",
            "Rora", "Rorc", "Dbp", "Nfil3", "Tef", "Hlf")
names(DS_bn) = unlist(group.l)
res_clock <- run_csvd(clocks, DS_bn, k_rem, color_variations, l.pos =c(0.8,1), g.s=1)

# -- Assemble ------------------------------------------------------------------

p_fig6 <- (p_per2 | res_clock$gene) /
  (p_dbp  | res_clock$tissue) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = list(c("A", "D", "B", "C")))

ggsave(file.path(plot_dir, "Fig6_clock_SVD.pdf"),
       p_fig6, width = 280, height = 297, units = "mm")


# ==============================================================================
# Figure 7: Pan-rhythmic genes
# ==============================================================================

# -- Fig 7A/B/C: Example pan-rhythmic genes ------------------------------------

epi_celltypes <- c("ISC", "TA/ISC_Proliferative", "Trans_Colonocytes",
                   "Colonocytes", "Trans_Goblet", "Goblet")

p_7a <- plot_data("Hspa8",  epi_celltypes, 11)
p_7b <- plot_data("Cirbp",  epi_celltypes, 11)
p_7c <- plot_data("Ndufa8", epi_celltypes, 11)

# -- Fig 7D/E: cSVD on pan-rhythmic genes --------------------------------------

pan_genes <- names(sort(table(unlist(lapply(DS.r, rownames))), decreasing = TRUE))
pan_genes <- pan_genes[!pan_genes %in% clocks]
pan_genes <- head(pan_genes, 200)
pan_genes <- pan_genes[!pan_genes %in% c('Muc2','Tppp3')]
res_pan <- run_csvd(pan_genes, DS_bn, k_rem, color_variations)

# -- Assemble Figure 7 --------------------------------------------------------

left_col  <- p_7a / p_7b / p_7c
right_col <- res_pan$gene / res_pan$tissue_cond

p_fig7 <- left_col | right_col +
  plot_layout(widths = c(2, 1))

# Need parentheses for correct operator precedence
p_fig7 <- (left_col | right_col) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = list(c("A", "B")))

ggsave(file.path(plot_dir, "Fig7_pan_rhythmic_SVD.pdf"),
       p_fig7, width = 350, height = 250, units = "mm")


# ==============================================================================
# Figure 5 & S6: GO enrichment per cell type group (dotplot with clustering)
# ==============================================================================

# -- Enrichment on rhythmic genes merged across conditions per celltype --------

DS.r <- lapply(DS_bn, function(x) subset(x, pval < 0.05 & amp > 0.5))
names(DS.r) <- unlist(group.l)

cell_types <- gsub(".+,", "", names(DS.r))
cell_type_groups <- split(seq_along(DS.r), cell_types)
DS.R_cell_types <- lapply(cell_type_groups, function(idxs) {
  unique(unlist(lapply(DS.r[idxs], rownames)))
})
DS.R_cell_types <- DS.R_cell_types[!names(DS.R_cell_types) %in% k_rem]

dbb <- c("KEGG_2019_Mouse", "WikiPathways_2024_Mouse", "Reactome_Pathways_2024")

GOL <- list()
for (i in seq_along(DS.R_cell_types)) {
  GOL[[i]] <- enrichr(DS.R_cell_types[[i]], databases = dbb)
}
names(GOL) <- names(DS.R_cell_types)

# -- Shared odds ratio range across all cell types -----------------------------

# Collect all results first to get global range
ALL_ct <- NULL
for (condi in names(GOL)) {
  for (f in names(GOL[[condi]])) {
    if (nrow(GOL[[condi]][[f]]) == 0) next
    g.kegg <- GOL[[condi]][[f]][, c("Term", "Adjusted.P.value", "Genes",
                                    "Overlap", "Odds.Ratio")]
    g.kegg$tot        <- as.numeric(sapply(strsplit(g.kegg$Overlap, "/"), "[[", 2))
    g.kegg$nb_gene    <- as.numeric(sapply(strsplit(g.kegg$Overlap, "/"), "[[", 1))
    g.kegg$odds_ratio <- pmin(g.kegg$Odds.Ratio, 20)
    g.kegg$type <- f
    g.kegg$cond <- condi
    colnames(g.kegg)[2] <- "adj.Pval"
    g.kegg <- subset(g.kegg, tot < 750 & nb_gene > 3)
    ALL_ct <- rbind(ALL_ct, g.kegg)
  }
}

odds_range_ct <- range(ALL_ct$odds_ratio[ALL_ct$adj.Pval < 0.1], na.rm = TRUE)
odds_breaks_ct <- c(2, 5, 10, 20)
odds_breaks_ct <- odds_breaks_ct[odds_breaks_ct <= max(odds_range_ct)]

# -- Lollipop/dotplot function with clustering and odds ratio ------------------

plot_go_multi <- function(condi_list, tt, divisor = 3) {
  All_cond <- NULL
  
  for (condi in condi_list) {
    if (!condi %in% names(GOL)) next
    g <- GOL[[condi]]
    
    for (f in names(g)) {
      if (nrow(g[[f]]) == 0) next
      
      g.kegg <- g[[f]][, c("Term", "Adjusted.P.value", "Genes",
                           "Overlap", "Odds.Ratio")]
      g.kegg$tot        <- as.numeric(sapply(strsplit(g.kegg$Overlap, "/"), "[[", 2))
      g.kegg$nb_gene    <- as.numeric(sapply(strsplit(g.kegg$Overlap, "/"), "[[", 1))
      g.kegg$odds_ratio <- pmin(g.kegg$Odds.Ratio, 20)
      g.kegg$type <- f
      g.kegg$cond <- condi
      colnames(g.kegg)[2] <- "adj.Pval"
      
      g.kegg <- subset(g.kegg, tot < 750 & nb_gene > 3)
      g.kegg <- g.kegg[!duplicated(g.kegg$Term), ]
      
      if (nrow(g.kegg) > 0) All_cond <- rbind(All_cond, g.kegg)
    }
  }
  
  if (is.null(All_cond) || nrow(All_cond) == 0) return(NULL)
  
  All_cond$Term <- gsub("\\s*(\\(GO:\\d+\\)|WP\\d+)", "", All_cond$Term)
  All_cond$Term <- trimws(All_cond$Term)
  All_cond$score <- pmin(-log10(All_cond$adj.Pval), 5)
  
  disease_kw <- c("disease", "cancer", "tumor", "leukemia", "carcinoma",
                  "infection", "virus", "viral", "pathogen", "neoplasm",
                  "spike", "anemia","Toxoplasmosis")
  All_cond <- All_cond[!grepl(paste(disease_kw, collapse = "|"),
                              All_cond$Term, ignore.case = TRUE), ]
  
  # Select significant terms (top 5 per db per celltype)
  top_terms <- All_cond %>%
    filter(as.numeric(adj.Pval) < 0.1) %>%
    group_by(type, cond) %>%
    arrange(adj.Pval) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    distinct(Term, type)
  
  All_cond_complete <- All_cond[All_cond$Term %in% top_terms$Term, ]
  if (nrow(All_cond_complete) == 0) return(NULL)
  
  # -- Hierarchical clustering and pruning ------------------------------------
  
  sig_terms <- unique(All_cond_complete$Term)
  
  all_genes <- unique(unlist(strsplit(All_cond_complete$Genes, ";")))
  term_gene_mat <- matrix(0, nrow = length(sig_terms), ncol = length(all_genes),
                          dimnames = list(sig_terms, all_genes))
  
  for (i in seq_len(nrow(All_cond_complete))) {
    genes <- unlist(strsplit(All_cond_complete$Genes[i], ";"))
    term_gene_mat[All_cond_complete$Term[i], intersect(genes, all_genes)] <- 1
  }
  
  term_gene_mat <- term_gene_mat[!duplicated(term_gene_mat), , drop = FALSE]
  sig_terms_clean <- rownames(term_gene_mat)
  
  if (length(sig_terms_clean) < 3) {
    term_order <- sig_terms_clean
    term_representatives <- sig_terms_clean
  } else {
    term_dist <- dist(term_gene_mat, method = "binary")
    term_clust <- hclust(term_dist, method = "ward.D2")
    
    n_clusters <- max(2, round(length(sig_terms_clean) / divisor))
    n_clusters <- min(n_clusters, length(sig_terms_clean))
    term_groups <- cutree(term_clust, k = n_clusters)
    
    term_representatives <- All_cond_complete %>%
      filter(Term %in% sig_terms_clean) %>%
      mutate(term_group = term_groups[as.character(Term)]) %>%
      group_by(term_group) %>%
      slice_min(adj.Pval, n = 1, with_ties = FALSE) %>%
      pull(Term) %>%
      unique()
    
    term_order <- sig_terms_clean[term_clust$order]
    term_order <- intersect(term_order, term_representatives)
  }
  
  All_cond_complete <- All_cond_complete[All_cond_complete$Term %in% term_representatives, ]
  
  # Wrap long names
  wrap_term <- function(term) {
    if (nchar(term) > 40) paste(strwrap(term, width = 40), collapse = "\n")
    else term
  }
  All_cond_complete$Term <- sapply(All_cond_complete$Term, wrap_term)
  term_order <- sapply(term_order, wrap_term)
  All_cond_complete$Term <- factor(All_cond_complete$Term, levels = term_order)
  
  # -- Plot -------------------------------------------------------------------
  ct_in_plot <- unique(All_cond_complete$cond)
  ct_labels <- setNames(
    paste0("<span style='color:", color_variations[ct_in_plot], "'>**",
           ct_in_plot, "**</span>"),
    ct_in_plot
  )
  
  ggplot(All_cond_complete,
         aes(y = Term, x = cond)) +
    geom_point(aes(size = odds_ratio, fill = score),
               shape = 21, stroke = 0.3, color = "black") +
    facet_grid(type ~ ., space = "free", scales = "free_y") +
    scale_fill_gradient(low = "white", high = "darkred",
                        limits = c(0, 5),
                        name = expression(-log[10] ~ p[adj])) +
    scale_size_continuous(range = c(1, 6), name = "Odds ratio",
                          limits = odds_range_ct,
                          breaks = odds_breaks_ct) +
    scale_x_discrete(labels = ct_labels) +
    labs(x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y          = element_text(size = tt),
      axis.text.x          = element_markdown(size = 10, angle = 45, hjust = 1),
      legend.position      = 'none',
      strip.text           = element_blank(),
      panel.background     = element_rect(fill = "grey95", color = NA),
      strip.background     = element_blank(),
      panel.grid.major     = element_line(color = "white", linewidth = 0.5),
      panel.grid.minor     = element_blank()
    )
}

# -- Save helper ---------------------------------------------------------------

save_go_panel <- function(plot, filename, width = 120, per_term_height = 8, min_height = 80) {
  n_terms <- length(levels(plot$data$Term))
  n_facets <- length(unique(plot$data$type))
  h <- max(min_height, n_terms * per_term_height + n_facets * 15 + 20)
  ggsave(file.path(plot_dir, filename), plot,
         width = width, height = h, units = "mm")
  cat(filename, ":", n_terms, "terms,", h, "mm height\n")
}

# ==============================================================================
# Figure 5: Epithelial and stromal rhythmic pathways
# ==============================================================================

p_5a <- plot_go_multi(c("ISC", "TA/ISC_Proliferative"), 11)
p_5b <- plot_data(c("Cdc6", "Atp5g1"), c("ISC", "TA/ISC_Proliferative"), 10)
p_5c <- plot_go_multi(c("Goblet", "Trans_Goblet", "Colonocytes", "Trans_Colonocytes"), 11)
p_5d <- plot_go_multi(c("Telocytes", "Trophocytes", "Interstitial"), 10)
p_5e <- plot_data(c("Pten", "Aifm2"), "Telocytes", 10)
p_5f <- plot_data(c("Grb2", "Git2"), "Trophocytes", 10)
p_5g <- plot_go_multi(c("Macrophages", "T", "T_Proliferative",
                        "B", "B_Proliferative", "Dendritic"), 11)

# GO panels
save_go_panel(p_5a, "Fig5_A_GO_ISC.pdf")
save_go_panel(p_5c, "Fig5_C_GO_epithelial.pdf")
save_go_panel(p_5d, "Fig5_D_GO_stromal.pdf")
save_go_panel(p_5g, "Fig5_G_GO_immune.pdf")

# Gene panels
ggsave(file.path(plot_dir, "Fig5_B_genes_ISC.pdf"),
       p_5b, width = 200, height = 100, units = "mm")
ggsave(file.path(plot_dir, "Fig5_E_genes_telocytes.pdf"),
       p_5e, width = 200, height = 100, units = "mm")
ggsave(file.path(plot_dir, "Fig5_F_genes_trophocytes.pdf"),
       p_5f, width = 200, height = 100, units = "mm")

# ==============================================================================
# Figure S6: GO enrichment for remaining cell types
# ==============================================================================

p_s6b <- plot_go_multi(c("MM", "MP", "LPM"), 11)
p_s6c <- plot_go_multi("Endothelial", 11)
p_s6d <- plot_go_multi("Lymphatic", 11)
p_s6e <- plot_go_multi("ICC", 11)
p_s6f <- plot_data(c("Tcf7l1", "Lrp5"), "ICC", 11)

save_go_panel(p_s6b, "FigS6_B_GO_smc.pdf")
save_go_panel(p_s6c, "FigS6_C_GO_endothelial.pdf")
save_go_panel(p_s6d, "FigS6_D_GO_lymphatic.pdf")
save_go_panel(p_s6e, "FigS6_E_GO_icc.pdf")

ggsave(file.path(plot_dir, "FigS6_F_genes_icc.pdf"),
       p_s6f, width = 200, height = 100, units = "mm")

cat("Figure 5 and S6 saved.\n")



# ==============================================================================
# Figure S7: Clock gene polar plots per gene and per cell type
# ==============================================================================

# -- Build polar data ----------------------------------------------------------

clocks_polar <- c("Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3",
                  "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Ciart", "Rora", "Rorc", "Dbp")

DSS.pv    <- t(do.call(rbind, lapply(DS_bn, function(x) x[clocks_polar, "pval"])))
DSS.amp   <- t(do.call(rbind, lapply(DS_bn, function(x) x[clocks_polar, "amp"])))
DSS.phase <- t(do.call(rbind, lapply(DS_bn, function(x) x[clocks_polar, "phase"])))

colnames(DSS.pv) <- colnames(DSS.amp) <- colnames(DSS.phase) <- unlist(group.l)
rownames(DSS.pv) <- rownames(DSS.amp) <- rownames(DSS.phase) <- clocks_polar

DSS.pv    <- melt(DSS.pv)
DSS.amp   <- melt(DSS.amp)
DSS.phase <- melt(DSS.phase)

DF.polar <- data.frame(DSS.pv, amp = DSS.amp$value, phase = DSS.phase$value)
DF.polar$cell_type <- gsub(".+,", "", DF.polar$Var2)
DF.polar$cond      <- gsub(",.+", "", DF.polar$Var2)
DF.polar$Var1      <- as.character(DF.polar$Var1)

DF.polar <- DF.polar[!DF.polar$cell_type %in% k_rem, ]
DF.polar <- subset(DF.polar, !is.na(amp) & !is.na(value))

DF.polar <- DF.polar %>%
  mutate(
    pval_category = case_when(
      -log10(value) < 1 ~ "< 1",
      -log10(value) < 2 ~ "1-2",
      -log10(value) < 3 ~ "2-3",
      TRUE ~ ">= 3"
    ),
    significant = value <= 0.1
  )

# -- Aesthetics ----------------------------------------------------------------

cond_colors <- c(
  Control_DC      = "#EC6166",
  Regenerating_DC = "#b8070d",
  Regenerating_PC = "#6C9FCA"
)

cond_shapes <- c(
  Control_DC      = 16,
  Regenerating_DC = 17,
  Regenerating_PC = 15
)

pval_sizes <- c("< 1" = 1.2, "1-2" = 2.2, "2-3" = 3.8, ">= 3" = 5.5)

polar_x <- scale_x_continuous(
  breaks = c(0, 6, 12, 18),
  labels = c("ZT00", "ZT06", "ZT12", "ZT18"),
  limits = c(0, 24),
  expand = c(0, 0)
)

polar_coord <- coord_polar(theta = "x", start = 0, direction = -1, clip = "off")

polar_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
    axis.title.x       = element_blank(),
    axis.text.x        = element_text(size = 10, color = "grey25"),
    axis.text.y        = element_text(size = 10, color = "grey35")
  )


# -- Plot functions ------------------------------------------------------------

plot_clock_gene <- function(gene_name) {
  
  df <- filter(DF.polar, Var1 == gene_name)

  ggplot(df, aes(phase, amp)) +
    geom_point(
      data = filter(df, significant),
      aes(col = cell_type, size = pval_category),
      alpha = 0.95
    ) +
    geom_point(
      data = filter(df, !significant),
      aes(col = cell_type, size = pval_category),
      shape = 1,
      stroke = 0.7,
      alpha = 0.45    ) +
    scale_color_manual(values = color_variations) +
    scale_size_manual(values = pval_sizes) +
    polar_x + polar_coord +
    facet_wrap(~ cond, ncol = 3, scales = 'free_y') +
    polar_theme +
    labs(title = gene_name, x = NULL, y = expression(Log[2] ~ fold-change))
}

plot_clock_celltype <- function(ct_name) {
  
  df <- filter(DF.polar, cell_type == ct_name)
  ggplot(df, aes(phase, amp)) +
    geom_point(
      data = filter(df, significant),
      aes(col = cond, size = pval_category),
      alpha = 0.95
    ) +
    geom_point(
      data = filter(df, !significant),
      aes(col = cond, size = pval_category),
      shape = 1,
      stroke = 0.7,
      alpha = 0.45
    ) +
    geom_text_repel(data = subset(df, Var1 %in% c("Arntl", "Nr1d1", "Dbp")),
      aes(label = Var1, col = cond, alpha = significant),
      size = 3,
      max.overlaps = 15,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15,
      segment.alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cond_colors) +
    scale_size_manual(values = pval_sizes) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
    polar_x + polar_coord +
    polar_theme +
    labs(title = ct_name, x = NULL, y = expression(Log[2] ~ fold-change))
}

# -- Build figure --------------------------------------------------------------

p_figS7 <- (
  plot_clock_gene("Arntl") /
    plot_clock_gene("Nr1d1") /
    plot_clock_gene("Dbp")
) / (
  plot_clock_celltype("MP") |
    plot_clock_celltype("Trans_Colonocytes") |
    plot_clock_celltype("TA/ISC_Proliferative")
) +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(plot_dir, "FigS7_clock_polar.pdf"),
  p_figS7,
  width = 420,
  height = 500,
  units = "mm"
)




