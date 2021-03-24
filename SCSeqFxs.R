# Functions for plotting and other single cell sequencing analysis
# SCSeqFxs.R

# plot_huang2021
# plot_filioDCMarkers

# get the name of an object
get_obj_name <- function(x) {
  print(as.character(substitute(x)))
}

# General use function to make plots with whatever set of markers you want
# Works best with 4 or 6 markers, but can change ncol to fit others better. 
# dataset is the seurat object
# objname is obj_name(dataset)
# marker_panel is the list of markers
# example function call:
# plot_markers(combined.DCs, get_obj_name(combined.DCs), c("Irf4", "Il12a", "Il15", "Cxcl12", "Cxcl16", "Ccl19", "Mki67", "Xcr1", "Eadem1"))
plot_markers <- function(dataset, objname, marker_panel) {
  # create subdirectory to put the outputs into
  foldername <- paste("custom_markers", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  # change default assay to RNA, but create a temporary variable to change it back at the end. 
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # Make plots
  # plot UMAPs for context
  pdf(paste(objname, "UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  
  # Violin plots
  pdf(paste(filenames, "_Violins.pdf", sep=""), width = 10, height = 10)
  print(VlnPlot(dataset, features = marker_panel, pt.size = 0, ncol = 2))
  dev.off()
  # Feature plots
  pdf(paste(filenames, "_FeaturePlot.pdf", sep=""), width = 10, height = 10)
  print(FeaturePlot(dataset, features = marker_panel, ncol = 2))
  dev.off()
  # Violin plots split by condition
  pdf(paste(filenames, "_SplitViolins.pdf", sep=""), width = 13, height = 10)
  plots <- VlnPlot(dataset, features = marker_panel, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}

# For all clusters, find the DE genes between conditions
# Example, what's different between tumor vs. egressed in cluster 3? (for each cluster)
# dataset is the seurat object
# objname is obj_name(dataset)
# example function call:
# FindMarkersByConditionEachCluster(combined.DCs, get_obj_name(combined.DCs))
FindMarkersByConditionEachCluster <- function(dataset, objname) {
  # create subdirectory to put the outputs into
  oldwd <- getwd()
  foldername <- paste("MarkersByCondition", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "integrated"
  
  # make a new metadata column "cluster_condition" and fill it with the cluster (Idents()) and condition info
  # Idents() is the seurat cluster, orig.idents is conditions (ie tumor or egressed)
  # examples: "0_egressed" or "0_tumor" or "1_egressed" or "1_tumor" etc...
  #dataset$cluster_condition <- paste(Idents(dataset), dataset$Condition, sep = "_")
  # make a new column called celltype to store the old cluster identities
  #dataset$celltype <- Idents(dataset)
  # Fill the identities column with the values in the new "cluster_condition" column
  #Idents(dataset) <- "cluster_condition"
  # now we can find the markers distinguishing the conditions within a cluster
  
  # cycle through all clusters
  # assumes cluster number starts at 0 and goes up as integers
  # assumes clusters have not been renamed... could get around this by just using seurat clusters... Implemented below
  ClusterList <- sort(c(unique(dataset$seurat_clusters)))-1
  # assumes only 2 conditions
  ConditionList <- c(unique(dataset$orig.ident))
  # ClusterList <- sort(c(unique(Idents(dataset))))-1
  # save a text file describing the conditions compared
  write.table(ConditionList[1:2], sep="\t", file="Condition1vsCondition2.txt", quote = FALSE, row.names = TRUE)
  for (ClusterID in ClusterList) {
    #FindMarkersByCondition(subset(dataset, Idents(dataset)==ClusterID), obj_name(dataset), ClusterID)
    # ^error, subset function doesn't really work programmatically at this point. Here's a workaround: 
    cells.use <- colnames(dataset)[which(dataset[[]]['seurat_clusters'] == ClusterID)]
    cell_subset <- subset(dataset, cells = cells.use)
    # make a new column called celltype to store the old cluster identities
    cell_subset$celltype <- Idents(cell_subset)
    # Fill the identities column with the values in the "Condition" column ("egressed" or "tumor")
    Idents(cell_subset) <- "Condition"
    # If the number of cells in a certain group is 3 or less, skip it and print a note.
    # cells.cond1.use <- colnames(cell_subset)[which(cell_subset[[]]['Condition'] == "egressed")]
    # cells.cond2.use <- colnames(cell_subset)[which(cell_subset[[]]['Condition'] == "tumor")]
    # try using orig.ident instead:
    cells.cond1.use <- colnames(cell_subset)[which(cell_subset[[]]['orig.ident'] == ConditionList[1])]
    cells.cond2.use <- colnames(cell_subset)[which(cell_subset[[]]['orig.ident'] == ConditionList[2])]
    if (length(cells.cond1.use) < 4 || length(cells.cond2.use) < 4) {
      print('Number of cells in a condition group of cluster: ')
      print(ClusterID)
      print('is 3 or less. Skipping this cluster.')
      next
    }
    # print('head(orig.ident)')
    # print(head(dataset$orig.ident, 30))
    FindMarkersByCondition(cell_subset, obj_name(dataset), ClusterID)
    
    DefaultAssay(dataset) <- tempdftassay
  }
  print('resetting to old wd: ')
  print(oldwd)
  setwd(oldwd)
}

# For one clusters, find the DE genes between conditions
# Example, what's different between tumor vs. egressed in cluster 3?
# data should already be subset to only include the cluster of interest. 
FindMarkersByCondition <- function(dataset, objname, ClusterID) {
  markers <- FindMarkers(dataset, ident.1 = "egressed", ident.2 = "tumor", min.pct = 0.25, verbose = FALSE)
  # save the results to a text file
  clustername <- paste("Cluster", ClusterID, sep="")
  write.table(markers, sep="\t", file=paste(clustername, "MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(markers['avg_log2FC'], sep="\t", file=paste(clustername, "Markers&LFCByCondition.rnk", sep="_"), quote = FALSE, row.names = TRUE, col.names = FALSE)
  # top 50 most up- and down- regulated genes
  top50 <- markers %>% top_n(n = 50, wt = avg_log2FC)
  bot50 <- markers %>% top_n(n = -50, wt = avg_log2FC)
  write.table(top50, sep="\t", file=paste(clustername, "top50UP_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(bot50, sep="\t", file=paste(clustername, "top50DOWN_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(row.names(top50), sep="\t", file=paste(clustername, "Names_top50UP_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(row.names(bot50), sep="\t", file=paste(clustername, "Names_top50DOWN_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
}

# Markers from Huang et al 2021 - Innervated Lymph Nodes ------------------
plot_huang2021 <- function(dataset, objname) {
  # create subdirectory to put the outputs into
  foldername <- paste("Huang2021_markers", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # plot UMAPs for context
  pdf(paste(objname, "UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  
  # Neutrophil panel
  neutrophil_panel1 <- c("Clec4d", "Csf3r", "Cxcr2", "Ngp", "Camp", "S100a9")
  pdf(paste(objname, "Neutrophil_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = neutrophil_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Neutrophil_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
  print(FeaturePlot(dataset, features = neutrophil_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "Neutrophil_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = neutrophil_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # CD4+ Cell panel
  CD4_panel1 <- c("Foxp3", "Ctla4", "Cd5", "Trbc2", "Cd4", "Lef1")
  pdf(paste(objname, "CD4_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = CD4_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "CD4_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
  print(FeaturePlot(dataset, features = CD4_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "CD4_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = CD4_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # NK Cell panel
  NK_panel1 <- c("Tox", "Cd8a", "Nkg7", "Birc5", "Kif11", "Ncr1", "Klrb1c")
  pdf(paste(objname, "NK_Panel_Violins1.pdf", sep="_"), width = 13, height = 10)
  print(VlnPlot(dataset, features = NK_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "NK_Panel_FeaturePlot1.pdf", sep="_"), width = 13, height = 10)
  print(FeaturePlot(dataset, features = NK_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "NK_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = NK_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Mitotic Cell panel
  Mitotic_panel1 <- c("Cenpe", "Tubb5", "Thy1", "Top2a", "Mki67", "", "")
  pdf(paste(objname, "Mitotic_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = Mitotic_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Mitotic_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
  print(FeaturePlot(dataset, features = Mitotic_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "Mitotic_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Mitotic_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # "Tissue T" Cell panel
  TissueTcells_panel1 <- c("Gata3", "Itgae", "Cxcr6", "Cd3d", "", "", "")
  pdf(paste(objname, "TissueTcells_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = TissueTcells_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "TissueTcells_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
  print(FeaturePlot(dataset, features = TissueTcells_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "TissueTcells_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = TissueTcells_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Mast Cell panel
  Mast_panel1 <- c("Gata2", "Kit", "Tpsb2", "Ccl9", "Fcgr3", "Csf1r", "C1qb", "Irf8")
  pdf(paste(objname, "Mast_Panel_Violins1.pdf", sep="_"), width = 13, height = 10)
  print(VlnPlot(dataset, features = Mast_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Mast_Panel_FeaturePlot1.pdf", sep="_"), width = 13, height = 10)
  print(FeaturePlot(dataset, features = Mast_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "Mast_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Mast_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # DC Cell panel
  DC_panel4 <- c("Clec9a", "Tlr11", "Xcr1", "Ccr7", "Arc", "Il15ra", "Ccl22", "")
  pdf(paste(objname, "DC_Panel_Violins1.pdf", sep="_"), width = 13, height = 10)
  print(VlnPlot(dataset, features = DC_panel4, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "DC_Panel_FeaturePlot1.pdf", sep="_"), width = 13, height = 10)
  print(FeaturePlot(dataset, features = DC_panel4, ncol = 3))
  dev.off()
  pdf(paste(objname, "DC_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = DC_panel4, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Random Cell panel
  Random_panel1 <- c("Ccr9", "Siglech", "Jchain", "Igkc", "Cd274", "", "", "")
  pdf(paste(objname, "Random_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = Random_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Random_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
  print(FeaturePlot(dataset, features = Random_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "Random_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Random_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Langerhans Cell panel
  Langerhans_panel1 <- c("Cd1d1", "Cd207", "Flt3l", "Batf3", "Irf8", "", "", "")
  pdf(paste(objname, "Langerhans_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = Langerhans_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Langerhans_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
  print(FeaturePlot(dataset, features = Langerhans_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "Langerhans_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Langerhans_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # DC5 Cell panel
  DC5_panel1 <- c("Axl", "Ptprc", "Csf1r", "", "", "", "", "")
  pdf(paste(obj_name(dataset), "DC5_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = DC5_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(obj_name(dataset), "DC5_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
  print(FeaturePlot(dataset, features = DC5_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "DC5_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = DC5_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}


# ------------------------------------------------
# plot Filio's DC Markers
# dataset is the dataset, objname is the name of the seurat object
# DC & Monocyte Marker Subsets from Filio...
# pDCs: Lin- (Ter119-(Ly76), Cd3-(Cd3e), Cd19-, Nk1.1-(Klrb1c))
#       Cd11c+(Itgax), Siglech+, CD45R/B220+(Ptprc),
plot_FilioDCMarkers <- function(dataset, objname) {
  # create subdirectory to put the outputs into
  foldername <- paste("FilioDC_markers", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # plot UMAPs for context
  pdf(paste(objname, "DC_UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "DC_UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  
  # Lineage Marker panel. Ter-119 is Ly76
  Lin_panel1 <- c("Ter-119", "Cd3e", "Cd19", "Klrb1c", "", "", "", "")
  pdf(paste(objname, "Lin_Panel_Violins1.pdf", sep="_"), width = 13, height = 10)
  print(VlnPlot(dataset, features = Lin_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Lin_Panel_FeaturePlot1.pdf", sep="_"), width = 13, height = 10)
  print(FeaturePlot(dataset, features = Lin_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "Lin_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Lin_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # pDC panel
  pDC_panel1 <- c("Itgax", "Siglech", "Ptprc", "Bst2", "Batf3", "Irf8", "Flt3l", "Irf7")
  pdf(paste(objname, "pDC_Panel_Violins1.pdf", sep="_"), width = 14, height = 12)
  print(VlnPlot(dataset, features = pDC_panel1, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "pDC_Panel_FeaturePlot1.pdf", sep="_"), width = 14, height = 12)
  print(FeaturePlot(dataset, features = pDC_panel1, ncol = 3))
  dev.off()
  pdf(paste(objname, "pDC_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = pDC_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # classical DC panel (cDC1 and cDC2)
  # cDC1:
  # CD11c hi (Itgax), MHCII hi (H2-Ab1), Cd24+, Cd8a+, Cd370+ (Clec9a), Cd11b low (Itgam), Cd172a low (Sirpa)
  # cDC2: 
  # CD11c hi (Itgax), MHCII hi (H2-Ab1), Cd24+, Cd8a low, Cd370 low (Clec9a), Cd11b+ (Itgam), Cd172a+ (Sirpa)
  # Tissue Resident DC: 2types?
  # one is: CD103+(Itgae), CD45+(Ptprc), CD11c high (Itgax), MHCII+ (H2-Ab)
  # other is: CD11b+(Itgam), CD45+(Ptprc), CD11c high (Itgax), MHCII+ (H2-Ab)
  cDC_panel1 <- c("Itgax", "H2-Ab1", "Cd24a", "Cd8a", "Clec9a", "Itgam", "Sirpa", "Itgae", "Ccr7")
  pdf(paste(objname, "cDC_Panel_Violins1.pdf", sep="_"), width = 14, height = 12)
  print(VlnPlot(dataset, features = cDC_panel1, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "cDC_Panel_FeaturePlot1.pdf", sep="_"), width = 14, height = 12)
  print(FeaturePlot(dataset, features = cDC_panel1, ncol = 3))
  dev.off()
  pdf(paste(objname, "cDC_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = cDC_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Neutrophil panel
  # Ly6G+, CD11b+(Itgam)
  Neutrophil_panel2 <- c("Ly6g", "Itgam", "Ly6c1", "Fcgr3", "", "", "", "")
  pdf(paste(objname, "Neutrophil_Panel2_Violins1.pdf", sep="_"), width = 13, height = 10)
  print(VlnPlot(dataset, features = Neutrophil_panel2, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Neutrophil_Panel2_FeaturePlot1.pdf", sep="_"), width = 13, height = 10)
  print(FeaturePlot(dataset, features = Neutrophil_panel2, ncol = 2))
  dev.off()
  pdf(paste(objname, "Neutrophil_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Neutrophil_panel2, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Monocyte panel
  # 2 types: Ly6C hi and low
  # Gr-1 low (Ly6g1), Cd115+(Csf1r), Cd11b+(Itgam), F4/80 int(Adgre1), Cxcr1+, Ccr2+, Cd43+(Spn)
  Mono_panel1 <- c("Ly6c1", "Ly6g", "Csf1r", "Itgam", "Adgre1", "Cxcr1", "Ccr2", "Spn")
  pdf(paste(objname, "Mono_Panel_Violins1.pdf", sep="_"), width = 14, height = 12)
  print(VlnPlot(dataset, features = Mono_panel1, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "Mono_Panel_FeaturePlot1.pdf", sep="_"), width = 14, height = 12)
  print(FeaturePlot(dataset, features = Mono_panel1, ncol = 3))
  dev.off()
  pdf(paste(objname, "Mono_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Mono_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots, nrow = 3))
  dev.off()
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}

# ------------------------------------------------
# plot Taylor's CD8 T-cell Markers
# Exploring and defining the CD8+ clusters ------------------------------------
# Taylor CD8+ T-cell subtype panels -------------------------------------------
# Vps37b, Icos, Tcf7, S1pr1, Cxcr6, Itgae, Jun, Fos, Nr4a3, Klf2, Ifng, Cd69
# Gzmb and Sell
# Eomes, Tbx21, Runx3, Id2, Id3, Tox, Il7r, Klrg1, Ccr7, Pdcd1, Havcr2, and Lag3
# dataset is the dataset, objname is the name of the seurat object
plot_TaylorCD8Markers <- function(dataset, objname) {
  # create subdirectory to put the outputs into
  foldername <- paste("TaylorCD8_markers", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # plot UMAPs for context
  pdf(paste(objname, "UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75))
  dev.off()
  
  # Taylor panel1
  Taylor_panel1 <- c("Vps37b", "Icos", "Tcf7", "S1pr1", "Cxcr6", "Cxcr4")
  pdf(paste(objname, "Taylor_Panel_Violins1.pdf", sep="_"), width = 13, height = 10)
  print(VlnPlot(dataset, features = Taylor_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_FeaturePlot1", sep="_"), width = 13, height = 10)
  print(FeaturePlot(dataset, features = Taylor_panel1, ncol = 2))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_SplitViolins1.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Taylor_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Taylor panel2
  Taylor_panel2 <- c("Jun", "Fos", "Nr4a3", "Klf2", "Ifng", "Cd69")
  pdf(paste(objname, "Taylor_Panel_Violins2.pdf", sep="_"), width = 14, height = 12)
  print(VlnPlot(dataset, features = Taylor_panel2, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_FeaturePlot2.pdf", sep="_"), width = 14, height = 12)
  print(FeaturePlot(dataset, features = Taylor_panel2, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_SplitViolins2.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Taylor_panel2, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Taylor panel3
  Taylor_panel3 <- c("Gzmb", "Sell", "Eomes", "Tbx21", "Runx3", "Tox")
  pdf(paste(objname, "Taylor_Panel_Violins3.pdf", sep="_"), width = 14, height = 12)
  print(VlnPlot(dataset, features = Taylor_panel3, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_FeaturePlot3.pdf", sep="_"), width = 14, height = 12)
  print(FeaturePlot(dataset, features = Taylor_panel3, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_SplitViolins3.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Taylor_panel3, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Taylor panel4
  Taylor_panel4 <- c("Id2", "Id3", "Il7r", "Klrg1", "Hsf1", "")
  pdf(paste(objname, "Taylor_Panel_Violins4.pdf", sep="_"), width = 13, height = 10)
  print(VlnPlot(dataset, features = Taylor_panel4, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_FeaturePlot4", sep="_"), width = 13, height = 10)
  print(FeaturePlot(dataset, features = Taylor_panel4, ncol = 2))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_SplitViolins4.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Taylor_panel4, split.by = "orig.ident", 
                   group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 2))
  dev.off()
  
  # Taylor panel5
  Taylor_panel5 <- c("Ccr7", "Pdcd1", "Havcr2", "Lag3", "", "")
  pdf(paste(objname, "Taylor_Panel_Violins5.pdf", sep="_"), width = 14, height = 12)
  print(VlnPlot(dataset, features = Taylor_panel5, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_FeaturePlot5.pdf", sep="_"), width = 14, height = 12)
  print(FeaturePlot(dataset, features = Taylor_panel5, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_SplitViolins5.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Taylor_panel5, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 2))
  dev.off()
  
  # Taylor panel6
  Taylor_panel6 <- c("Cd44", "Mki67", "Fosb", "Nr4a1", "Ccr7", "")
  pdf(paste(objname, "Taylor_Panel_Violins6.pdf", sep="_"), width = 14, height = 12)
  print(VlnPlot(dataset, features = Taylor_panel6, pt.size = 0, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_FeaturePlot6.pdf", sep="_"), width = 14, height = 12)
  print(FeaturePlot(dataset, features = Taylor_panel6, ncol = 3))
  dev.off()
  pdf(paste(objname, "Taylor_Panel_SplitViolins6.pdf", sep="_"), width = 18, height = 14)
  plots <- VlnPlot(dataset, features = Taylor_panel6, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 2))
  dev.off()
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}

