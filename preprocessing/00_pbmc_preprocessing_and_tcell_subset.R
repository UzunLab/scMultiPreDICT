# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Signac)
library(readr)
library(Matrix)
library(patchwork)
library(rtracklayer)
library(biovizBase)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

# ============================================================================
# CONFIGURATION - UPDATE THESE PATHS FOR YOUR SYSTEM
# ============================================================================
# Download the PBMC 10k Multiome dataset from 10x Genomics:
# https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0

# Path to the filtered feature barcode matrix folder
INPUT_MATRIX_DIR <- "/path/to/filtered_feature_bc_matrix/"

# Define data file paths
mtx <- file.path(INPUT_MATRIX_DIR, "matrix.mtx.gz")
features <- file.path(INPUT_MATRIX_DIR, "features.tsv.gz")
barcodes <- file.path(INPUT_MATRIX_DIR, "barcodes.tsv.gz")

# Path to ATAC fragments file (must be sorted and indexed with .tbi file)
fragments_path <- "/path/to/atac_fragments.tsv.gz"

# Output directory for T cell subset
OUTPUT_DIR <- "~/scMultiPreDICT_output/preprocessing/"

# Load RNA
counts_all <- ReadMtx(
  mtx = mtx,
  features = features,
  cells = barcodes
)

# load features 
feat <- read_tsv(features, col_names = FALSE, show_col_types = FALSE)

# Split into RNA and ATAC by feature type
# The features.tsv is expected to have at least 3 columns: feature_id, feature_name, feature_type
rna_rows <- which(feat$X3 %in% c("Gene Expression", "Gene expression", "Gene"))
atac_rows <- which(feat$X3 %in% c("Peaks", "Peak"))

if (length(rna_rows) == 0) stop("No RNA features found in features file (check column 3 values)")

rna_counts <- counts_all[rna_rows, , drop = FALSE]
atac_counts <- counts_all[atac_rows, , drop = FALSE]

# Create Seurat Object for RNA
seurat_obj_full <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = "Human_PBMC_10k"
)

# check dim
dim(seurat_obj_full)

# Extract genome annotations from EnsDB
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- "UCSC"
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
genome(annotations) <- "hg38"


frags <- CreateFragmentObject(path = fragments_path, cells = colnames(seurat_obj_full), validate.fragments = TRUE)

# Create chromatin assay and add to the seurat_obj
seurat_obj_full[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  fragments = frags,
  genome = "hg38",
  sep = c(":", "-"),
  annotation = annotations
)

# check dim
dim(seurat_obj_full)


#======================Perform QC on ATAC=========================================
# Check default assay
DefaultAssay(seurat_obj_full)

# Change default assay to ATAC
DefaultAssay(seurat_obj_full) <- "ATAC"

seurat_obj_full <- NucleosomeSignal(object = seurat_obj_full)

# look at the fragment length periodicity for all cells
seurat_obj_full$nucleosome_group <- ifelse(seurat_obj_full$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = seurat_obj_full, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

# Calculate TSS Enrichment Scores
seurat_obj_full <-  TSSEnrichment(seurat_obj_full)



# Set assay to RNA and add percentage of reads that map to the mitochondrial genome (will later filter reads with high mitochondria counts out)
DefaultAssay(seurat_obj_full) <-  "RNA"
seurat_obj_full[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_full, pattern = "^MT-") 

dim(seurat_obj_full)


# Manually PLot to view

# Load the ggplot library
library(ggplot2) 

# Extrat metadata features to plot
qc_df <- seurat_obj_full@meta.data 

# nCount_RNA 
ggplot(qc_df, aes(x = "", y = nCount_RNA)) + 
  geom_violin() + 
  labs(x = NULL, y = "nCount_RNA") + 
  theme_bw() 

# percent_mt 
ggplot(qc_df, aes(x = "", y = percent.mt)) + 
  geom_violin() + 
  labs(x = NULL, y = "percent.mt") + 
  theme_bw()

# nFeature_RNA 
ggplot(qc_df, aes(x = "", y = nFeature_RNA)) + 
  geom_violin() + 
  labs(x = NULL, y = "nFeature_RNA") + 
  theme_bw()


# PLot side by side
library(tidyr)
library(dplyr)
library(ggplot2) 


df <- FetchData(seurat_obj_full, vars = c("nCount_RNA", "percent.mt", "nFeature_RNA")) 

df$group <- Idents(seurat_obj_full) 


df_long <- df %>% 
  mutate(.cell = rownames(.)) %>% 
  pivot_longer(cols = c("nCount_RNA", "percent.mt", "nFeature_RNA"), 
               names_to = "metric", values_to = "value") %>% 
  dplyr::filter(is.finite(value)) 

ggplot(df_long, aes(x = group, y = value)) + 
  geom_violin(trim = TRUE, scale = "width") + 
  geom_jitter(width = 0.1, size = 0.2, alpha = 0.4) + 
  facet_wrap(~ metric, scales = "free_y", ncol = 4) 



# Filter low quality cells
seurat_obj_full <- subset(seurat_obj_full, subset = nFeature_RNA > 500 & nFeature_RNA < 3200 & percent.mt < 15)

#check dimensions
dim(seurat_obj_full)

seurat_obj_full <- subset(
  x = seurat_obj_full,
  subset =
    nCount_RNA < 9000 &
    nCount_RNA > 500 &
    nCount_ATAC < 90000 &
    nCount_ATAC > 4000 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 4
)


dim(seurat_obj_full)

# Normalization and Dimensionality Reduction

seurat_obj_full <- NormalizeData(seurat_obj_full)

# select highly variable genes
seurat_obj_full <- FindVariableFeatures(seurat_obj_full,  selection.method = "vst", nfeatures = 3000)

# scale data for PCA
all.genes <- rownames(seurat_obj_full)
seurat_obj_full <- ScaleData(seurat_obj_full, features = all.genes)


# Carry out PCA
seurat_obj_full <- RunPCA(seurat_obj_full, features = VariableFeatures(object = seurat_obj_full))
ElbowPlot(seurat_obj_full)

# Clustering
seurat_obj_full <- FindNeighbors(seurat_obj_full, dims = 1:10)
seurat_obj_full <- FindClusters(seurat_obj_full, resolution = 0.5)

# Generate UMPA for Visualization purpsoses
seurat_obj_full <- RunUMAP(seurat_obj_full, dims = 1:10, verbose = F)

# Take a look at cluster sizes
table(seurat_obj_full@meta.data$seurat_clusters)

# Visualize clsuters
DimPlot(seurat_obj_full, label.size =  4, repel = T, label = T)

# Differential Expression and  Marker selection
all.markers <- FindAllMarkers(seurat_obj_full, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

# Take a quick glance at the markers
dim(all.markers)
table(all.markers$cluster)
top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

# Cell Type annotation using SingleR
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()
monaco.ref <- celldex::MonacoImmuneData()

# Convert seurat obj to single cell experiment for convenience use of SingleR for cell annotation
sce <- as.SingleCellExperiment(seurat_obj_full)
sce

# Carry out annotations using different databases
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)

# View  Table
table(hpca.main$pruned.labels)
table(dice.main$pruned.labels)
table(monaco.main$pruned.labels)

# Add annotations to seurat obj
seurat_obj_full@meta.data$hpca.main   <- hpca.main$pruned.labels
seurat_obj_full@meta.data$dice.main   <- dice.main$pruned.labels
seurat_obj_full@meta.data$monaco.main <- monaco.main$pruned.labels
seurat_obj_full@meta.data$hpca.fine   <- hpca.fine$pruned.labels
seurat_obj_full@meta.data$dice.fine   <- dice.fine$pruned.labels
seurat_obj_full@meta.data$monaco.fine <- monaco.fine$pruned.labels

# Save seurat obj with cell annotations
seurat_output_dir <- file.path(OUTPUT_DIR, "seurat_obj")
dir.create(seurat_output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(seurat_obj_full, file.path(seurat_output_dir, "Human_PBMC_With_Cell_Annotations_Multiome_processed.rds"))

# Read in the saved seurat obj
seurat_obj <- readRDS(file.path(seurat_output_dir, "Human_PBMC_With_Cell_Annotations_Multiome_processed.rds"))

# Visualise cell annotations
seurat_obj_sub <- subset(seurat_obj, subset = !is.na(monaco.fine))
print(DimPlot(seurat_obj_sub, reduction = "umap", group.by = 'monaco.fine', repel = T, label = TRUE))

# Inspect Label and Counts
table(seurat_obj$hpca.main)
unique(seurat_obj$hpca.main)

# replace "T cell" below with the exact label(s) from the table above
t_labels <- c("T_cells") 
tcells_barcodes <- colnames(seurat_obj)[seurat_obj$hpca.main %in% t_labels]
length(tcells_barcodes)  # how many cells

mono_labels <- c("Monocyte")
mcells_barcodes <- colnames(seurat_obj)[seurat_obj$hpca.main %in% mono_labels]
length(mcells_barcodes)

# Save barcodes for downstream analysis
barcodes_output_dir <- file.path(OUTPUT_DIR, "celltype_specific_barcodes")
dir.create(barcodes_output_dir, recursive = TRUE, showWarnings = FALSE)

write.table(tcells_barcodes, file=file.path(barcodes_output_dir, "Tcells_barcodes.txt"),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(mcells_barcodes, file=file.path(barcodes_output_dir, "monocytes_barcodes.txt"),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# 1. Force join layers for the RNA assay 

# This ensures "counts" is one unified matrix containing ALL cells 
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]]) 

# 2. Update your cells_to_keep list 
# (Re-running this to be absolutely sure) 
cells_to_keep <- colnames(seurat_obj)[seurat_obj$hpca.main == "T_cells"] 

# Debug check: Does the object actually have these cells? 
print(paste("Cells found:", length(cells_to_keep))) 

# 3. Get the unified counts matrix 
# We use slot = "counts" to get raw data (compatible across Seurat versions)
raw_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts") 

# 4. CRITICAL STEP: Find the intersection 
# This prevents "subscript out of bounds" by only asking for cells that exist in BOTH lists 
valid_cells <- intersect(colnames(raw_counts), cells_to_keep) 

# Check how many we lost (should be 0 if JoinLayers worked) 
print(paste("Original request:", length(cells_to_keep))) 
print(paste("Valid cells in matrix:", length(valid_cells))) 

# 5. Now subset the matrix 
subset_counts <- raw_counts[, valid_cells] 

# 6. Create the fresh object 
seurat_subset <- CreateSeuratObject( 
  counts = subset_counts, 
  meta.data = seurat_obj@meta.data[valid_cells, ] 
) 

# Get ATAC Data
raw_atac <- GetAssayData(seurat_obj, assay = "ATAC", layer = "counts")
subset_atac <- raw_atac[, valid_cells]

# Create a fragments object restricted to the T-cell barcodes (faster and safer)
frags_subset <- CreateFragmentObject(path = fragments_path, cells = valid_cells, validate.fragments = TRUE)

chrom_assay <-  CreateChromatinAssay(
  counts = subset_atac,
  fragments = frags_subset,
  sep = c(":","-"),
  genome = "hg38"
)

# 7. Attach the fragments (Required for calling peaks later) 
seurat_subset[["ATAC"]] <- chrom_assay
class(seurat_subset)
class(seurat_subset[["ATAC"]])

print("New clean object created successfully.") 

print(seurat_subset) 


# Call peaks
peaks_output_dir <- file.path(OUTPUT_DIR, "peaks_t_cells_output")
dir.create(peaks_output_dir, recursive = TRUE, showWarnings = FALSE)

peaks_subset <- CallPeaks(
  object = seurat_subset,
  assay = "ATAC",
  name = "t_cells",
  macs2.path = "/ri/shared/modules8/MACS2/2.2.9.1/bin/macs2",  # Update this path for your system
  outdir = peaks_output_dir
)

peaks_subset <- keepStandardChromosomes(peaks_subset, pruning.mode = "coarse")
if (exists("blacklist_hg38_unified")) {
  peaks_subset <- subsetByOverlaps(x = peaks_subset, ranges = blacklist_hg38_unified, invert = TRUE)
} else {
  message("blacklist_hg38_unified not found — skipping blacklist filtering")
}

frag <- Fragments(seurat_subset[["ATAC"]])[[1]]

peak_counts <- FeatureMatrix(
  fragments = frag,
  features = peaks_subset,
  cells = colnames(seurat_subset)
)

seurat_subset[["ATAC"]] <- CreateChromatinAssay(
  counts = peak_counts,
  fragments = frag
)

DefaultAssay(seurat_subset) <- "ATAC"

class(seurat_subset)
class(seurat_subset[["ATAC"]])
Fragments(seurat_subset[["ATAC"]])
table(seurat_subset$hpca.main)
ncol(seurat_subset)

seurat_subset@assays$ATAC@counts
head(rownames(seurat_subset[["ATAC"]]))

# Filter genes after cell filtering in at least 10% of cells
cat("\n=== Filtering low-expression genes ===\n\n")
DefaultAssay(seurat_subset) <- "RNA"
num_cells <- ncol(seurat_subset)
gene_counts <- GetAssayData(seurat_subset, slot = "counts")
gene_frac <- Matrix::rowSums(gene_counts > 0) / num_cells
genes_keep <- names(gene_frac[gene_frac >= 0.10])
n_genes_before <- nrow(seurat_subset[["RNA"]])
seurat_subset[["RNA"]] <- subset(seurat_subset[["RNA"]], features = genes_keep)
n_genes_after <- nrow(seurat_subset[["RNA"]])
cat(sprintf("Gene filter: Keep genes expressed in >= 10%% of filtered cells\n"))
cat(sprintf("Genes before: %d\n", n_genes_before))
cat(sprintf("Genes after:  %d\n", n_genes_after))
cat(sprintf("Genes removed: %d (%.1f%%)\n", n_genes_before - n_genes_after, 100 * (n_genes_before - n_genes_after) / n_genes_before))

dim(seurat_subset)
Assays(seurat_subset)

dim(peak_counts) 
all(colnames(peak_counts)%in% valid_cells)
head(colnames(Fragments(seurat_subset[["ATAC"]])[[1]]))


saveRDS(seurat_subset, file.path(seurat_output_dir, "T_cells_seurat_multiome_processed.rds"))

message("T cell subset saved to: ", file.path(seurat_output_dir, "T_cells_seurat_multiome_processed.rds"))
