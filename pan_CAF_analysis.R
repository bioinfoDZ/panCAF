#CAF Integration
mel_human_caf <- readRDS(file.choose())
Melanoma_cafs_ids <- c("Mel CAFs 1", "Mel CAFs 2")
names(Melanoma_cafs_ids) <- levels(mel_human_caf)
mel_human_caf <-RenameIdents(mel_human_caf, Melanoma_cafs_ids)
DimPlot(mel_human_caf)

hnsc_human_caf <- readRDS(file.choose())
HNSCC_cafs_ids <- c("HNSCC CAFs 1", "HNSCC CAFs 2", "HNSCC CAFs 3", "HNSCC CAFs 4", "HNSCC CAFs 5")
names(HNSCC_cafs_ids) <- levels(hnsc_human_caf)
hnsc_human_caf <- RenameIdents(hnsc_human_caf, HNSCC_cafs_ids)
DimPlot(hnsc_human_caf)

lung_human_caf <- readRDS(file.choose())
Lung_cafs_ids <- c("LC CAFs 1", "LC CAFs 2", "LC CAFs 3", "LC CAFs 4")
names(Lung_cafs_ids) <- levels(lung_human_caf)
lung_human_cancer_resolution0.1 <-RenameIdents(lung_human_caf, Lung_cafs_ids)
DimPlot(lung_human_caf)

caf.anchors <- FindIntegrationAnchors(object.list = list(mel_human_caf, hnsc_human_caf, lung_human_caf), dims = 1:20, k.filter = 150)
caf.combined <- IntegrateData(anchorset = caf.anchors, dims = 1:20)
DefaultAssay(caf.combined) <- "integrated"
caf.combined <- ScaleData(caf.combined, verbose = FALSE)
caf.combined <- RunPCA(caf.combined, features = VariableFeatures(object = caf.combined))
DimPlot(caf.combined, reduction = "pca")
DimHeatmap(caf.combined, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 11:20, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 21:30, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 31:40, cells = 500, balanced = TRUE)
DimHeatmap(caf.combined, dims = 41:50, cells = 500, balanced = TRUE)
caf.combined <- JackStraw(caf.combined, num.replicate = 100)
caf.combined <- ScoreJackStraw(caf.combined, dims = 1:20)
JackStrawPlot(caf.combined, dims = 1:20)
ElbowPlot(caf.combined, ndims = 50)
caf.combined <- RunUMAP(caf.combined, reduction = "pca", dims = 1:30)
caf.combined <- FindNeighbors(caf.combined, reduction = "pca", dims = 1:30)
#Resolution 0.2 PCs 1 through 30
caf.combined.resolution.02 <- FindClusters(caf.combined, resolution = 0.2)
DimPlot(caf.combined.resolution.02, reduction = "umap", label = TRUE)
caf.combined.resolution.02.markers <- FindAllMarkers(caf.combined.resolution.02, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
caf.combined.resolution.02.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
caf.combined.resolution.02.markers.top20 <- caf.combined.resolution.02.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(caf.combined.resolution.02, features = caf.combined.resolution.05.markers.top20$gene, label = F) + NoLegend()
integrated_cafs_ids <- c("Pan-myCAF", "Pan-dCAF", "Pan-iCAF", "Pan-iCAF-2", "Pan-nCAF", "LQ-CAF", "Pan-pCAF")
names(integrated_cafs_ids) <- levels(caf.combined.resolution.02)
caf.combined.resolution.02 <-RenameIdents(caf.combined.resolution.02, integrated_cafs_ids)
features <- c('ACTA2', 'MYH11', 'TAGLN', 'MCAM', 'MYLK', 'COL1A1', 'COL3A1', 'COL10A1', 'MMP1', 'MMP11', 'STC1', 'CFD', 'C3', 'CXCL14', 'CXCL12', 'BDKRB1', 'CLU', 'CXCL2', 'ICAM1', 'TNFAIP3', 'APOC1', 'TPD52L1', 'TPD52', 'CXCR4', 'CCL19', 'C10orf10', 'BIRC5', 'TOP2A', 'CDK1', 'CDC25C', 'CDC45')
DoHeatmap(caf.combined.resolution.02, features = 'STC1, IGF1', label = F) + NoLegend()
#cell cycle analysis of 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
l1 <- CellCycleScoring(caf.combined.resolution.02, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(caf.combined.resolution.02[[]])
DimPlot(l1, reduction = "umap", group.by = 'Phase')
#Qian et al functional gene sets
pan_caf_file <- read.table(file.choose(), header = T, row.names = 1, sep = '\t')
proteins <- c('COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL6A1', 'COL7A1', 'COL8A1', 'COL10A1', 'COL11A1', 'COL12A1', 'COL13A1', 'COL14A1', 'COL15A1', 'COL16A1', 'COL18A1', 'BGN', 'DCN', 'LUM', 'TAGLN','ELN', 'FN1', 'MMP1', 'MMP2', 'MMP3', 'MMP9', 'MMP10', 'MMP11', 'MMP14', 'MMP19', 'SERPINE1', 'CTHRC1', 'THBS2', 'SULF1', 'TGFBI', 'COMP', 'INHBA', 'EGFL6', 'ANGPT2', 'PDGFA', 'PDGFC', 'VEGFA', 'ACTA2', 'MYL6', 'MYH9', 'MYH11', 'PLN', 'TPM1', 'TMP2', 'SORBS2', 'RRAS', 'RASL12', 'RASGRP2', 'CFD', 'CFI', 'C3', 'C7', 'CCL21', 'CXCL14', 'CXCL12', 'IL33', 'CXCL3', 'CXCL2', 'CXCL1', 'CCL2', 'CCXL26',  'IL6', 'IL7')
pan_caf_file_proteins_for_heatmap <- pan_caf_file[match(proteins, row.names(pan_caf_file)),]
pan_caf_file_proteins_for_heatmap <- na.omit(pan_caf_file_proteins_for_heatmap)
pheatmap(pan_caf_file_proteins_for_heatmap[,-6], scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = T, show_colnames = F,  color=colorRampPalette(c("blue","white","red"))(50))
#Cell surface
cell_surface_protein_atlas <- read.table(file.choose(), header = T, row.names = 1)
pan_caf_marker_genes <- read.table(file.choose(), header = T, row.names = 1)
pan_caf_marker_genes_cell_surface <- pan_caf_marker_genes[match(row.names(cell_surface_protein_atlas), pan_caf_marker_genes$gene),]
pan_caf_marker_genes_cell_surface <- na.omit(pan_caf_marker_genes_cell_surface)
write.table(pan_caf_marker_genes_cell_surface, file = '/Users/phillipgalbo/Desktop/pan_caf_marker_genes_cell_surface.txt', sep = '\t', col.names = T, row.names = T, quote = F)
cell_surface_markers <- c('PARM1', 'APOC3', 'CSPG4', 'CD36', 'SUSD2', 'SDC1', 'TREM1', 'THY1', 'BAMBI', 'P2RY6', 'CD34', 'PI16', 'GPC3', 'RAMP2', 'LRRN3', 'CLDN1', 'ICAM1', 'MUSK', 'CXCR4', 'NCAM1', 'EFNB1', 'EBP')
pan_caf_file <- read.table(file.choose(), header = T, row.names = 1, sep = '\t')
colnames_to_declare <- c('Pan-myCAF', 'Pan-dCAF', 'Pan-iCAF', 'Pan-iCAF-2', 'Pan-nCAF', 'Pan-nCAF-2', 'Pan-pCAF')
colnames(pan_caf_file) <- colnames_to_declare
pan_caf_file_receptors <- pan_caf_file[match(cell_surface_markers, row.names(pan_caf_file)),]
pan_caf_file_receptors <- na.omit(pan_caf_file_receptors)
pan_caf_file_receptors <- pan_caf_file_receptors[, -6]
pheatmap(pan_caf_epic_file_receptors, scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = T, angle_col = c('45'), color=colorRampPalette(c("blue","white","red"))(20))
#TF analysis
integrated <- readRDS(file.choose()) #load integared pan-caf seurat analysis
integrated_pan_caf_ids <- c("Pan-mCAF", "Pan-dCAF", "Pan-iCAF", "Pan-iCAF-2", "Pan-nCAF", "LQ-CAF", "Pan-pCAF")
names(integrated_pan_caf_ids) <- levels(integrated)
integrated <- RenameIdents(integrated, integrated_pan_caf_ids)
DimPlot(integrated, reduction = "umap", label = TRUE)
integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
integrated.markers.top20 <- integrated.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(integrated, features = integrated.markers.top20$gene) + NoLegend()

tfs <- readRDS(file.choose()) #load Geneie3 RDS file from SCENIC analysis
marker_genes <- read.table(file.choose(), header = T, row.names = 1, sep = '\t')
marker_genes_mCAF <- subset(marker_genes, marker_genes$cluster == 0)
meft2C_tf_target <- subset(tfs, tfs$TF == 'MEF2C')
marker_genes_mef2c_target_genes <- marker_genes_mCAF[match(meft2C_tf_target$Target, marker_genes_mCAF$gene),]
marker_genes_mef2c_target_genes <- subset(marker_genes_mef2c_target_genes, marker_genes_mef2c_target_genes$avg_logFC >= 1 & marker_genes_mef2c_target_genes$p_val_adj < 0.05)
marker_genes_mef2c_target_genes <- na.omit(marker_genes_mef2c_target_genes)
mef2c_target_genes <- marker_genes_mef2c_target_genes$gene
pan_cafs_epic <- read.table(file.choose(), header = T, row.names = 1, sep = '\t')
pan_cafs_epic_mef2c_target_genes <- pan_cafs_epic[match(mef2c_target_genes, row.names(pan_cafs_epic)),]
pheatmap(pan_cafs_epic_mef2c_target_genes[,-6], scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = F,  color=colorRampPalette(c("blue","white","red"))(20))

marker_genes_dCAF <- subset(marker_genes, marker_genes$cluster == 1)
twist1_tf_target <- subset(tfs, tfs$TF == 'TWIST1')
marker_genes_twist1_target_genes <- marker_genes_dCAF[match(twist1_tf_target$Target, marker_genes_dCAF$gene),]
marker_genes_twist1_target_genes <- subset(marker_genes_twist1_target_genes, marker_genes_twist1_target_genes$avg_logFC >= 1 & marker_genes_twist1_target_genes$p_val_adj < 0.05)
marker_genes_twist1_target_genes <- na.omit(marker_genes_twist1_target_genes)
twist1_target_genes <- marker_genes_twist1_target_genes$gene
pan_cafs_epic_twist1_target_genes <- pan_cafs_epic[match(twist1_target_genes, row.names(pan_cafs_epic)),]
pheatmap(pan_cafs_epic_twist1_target_genes[,-6], scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = F,  color=colorRampPalette(c("blue","white","red"))(20))
pheatmap(pan_cafs_epic_twist1_target_genes, scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = T, fontsize = .5,  color=colorRampPalette(c("blue","white","red"))(20))

marker_genes_iCAF <- subset(marker_genes, marker_genes$cluster == 2)
nr1h3_tf_target <- subset(tfs, tfs$TF == 'NR1H3')
marker_genes_nr1h3_target_genes <- marker_genes_iCAF[match(nr1h3_tf_target$Target, marker_genes_dCAF$gene),]
marker_genes_nr1h3_target_genes <- subset(marker_genes_nr1h3_target_genes, marker_genes_nr1h3_target_genes$avg_logFC >= 1 & marker_genes_nr1h3_target_genes$p_val_adj < 0.05)
marker_genes_nr1h3_target_genes <- na.omit(marker_genes_nr1h3_target_genes)
nr1h3_target_genes <- marker_genes_nr1h3_target_genes$gene
pan_cafs_epic_nr1h3_target_genes <- pan_cafs_epic[match(nr1h3_target_genes, row.names(pan_cafs_epic)),]
pheatmap(pan_cafs_epic_nr1h3_target_genes[,-6], scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = F,  color=colorRampPalette(c("blue","white","red"))(20))
pheatmap(pan_cafs_epic_nr1h3_target_genes, scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = T, fontsize = .5,  color=colorRampPalette(c("blue","white","red"))(20))

marker_genes_iCAF2 <- subset(marker_genes, marker_genes$cluster == 3)
relb_tf_target <- subset(tfs, tfs$TF == 'RELB')
marker_genes_relb_target_genes <- marker_genes_iCAF2[match(relb_tf_target$Target, marker_genes_iCAF2$gene),]
marker_genes_relb_target_genes <- subset(marker_genes_relb_target_genes, marker_genes_relb_target_genes$avg_logFC >= 1 & marker_genes_relb_target_genes$p_val_adj < 0.05)
marker_genes_relb_target_genes <- na.omit(marker_genes_relb_target_genes)
relb_target_genes <- marker_genes_relb_target_genes$gene
pan_cafs_epic_relb_target_genes <- pan_cafs_epic[match(relb_target_genes, row.names(pan_cafs_epic)),]
pheatmap(pan_cafs_epic_relb_target_genes[,-6], scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = F,  color=colorRampPalette(c("blue","white","red"))(20))

marker_genes_pCAF <- subset(marker_genes, marker_genes$cluster == 6)
foxm1_tf_target <- subset(tfs, tfs$TF == 'FOXM1')
marker_genes_foxm1_target_genes <- marker_genes_pCAF[match(foxm1_tf_target$Target, marker_genes_pCAF$gene),]
marker_genes_foxm1_target_genes <- subset(marker_genes_foxm1_target_genes, marker_genes_foxm1_target_genes$avg_logFC >= 1 & marker_genes_foxm1_target_genes$p_val_adj < 0.05)
marker_genes_foxm1_target_genes <- na.omit(marker_genes_foxm1_target_genes)
foxm1_target_genes <- marker_genes_foxm1_target_genes$gene
pan_cafs_epic_foxm1_target_genes <- pan_cafs_epic[match(foxm1_target_genes, row.names(pan_cafs_epic)),]
pheatmap(pan_cafs_epic_foxm1_target_genes[,-6], scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = F,  color=colorRampPalette(c("blue","white","red"))(20))
pheatmap(pan_cafs_epic_foxm1_target_genes, scale = "row", cluster_rows = F, cluster_cols = F,  show_rownames = T, fontsize = .5,  color=colorRampPalette(c("blue","white","red"))(20))

tfs_names <- unique(tfs$TF)
DoHeatmap(integrated, features = tfs_names) + NoLegend()
tfs_data_frame <- data.frame(tfs)
write.table(tfs_data_frame,file = '/Users/phillipgalbo/Desktop/tfs_data_frame.txt', col.names = T, row.names = T, quote = F, sep = '\t')

mef2c <- subset(tfs_data_frame, tfs_data_frame$TF == 'MEF2C')
DoHeatmap(integrated, features = mef2c$Target) + NoLegend()

regulon_auc <- readRDS(file.choose())
VlnPlot(integrated, features = 'MEF2C', pt.size = F)
VlnPlot(integrated, features = 'TWIST1', pt.size = F)
VlnPlot(integrated, features = 'NR1H3', pt.size = F)
VlnPlot(integrated, features = 'RELB', pt.size = F)
VlnPlot(integrated, features = 'FOXM1', pt.size = F)

