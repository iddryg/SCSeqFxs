# Functions for plotting and other single cell sequencing analysis
# SCSeqFxs.R


# get_obj_name
# plot_markers
# FindMarkersByConditionEachCluster
# FindMarkersByCondition
# compare_clusters
# plot_huang2021
# plot_filioDCMarkers
# Plot_TaylorMarkers
# doPseudotime
# doGSEA
# doAnalysis


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
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  
  # Violin plots
  pdf(paste(filenames, "_Violins.pdf", sep=""), width = 12, height = 10)
  print(VlnPlot(dataset, features = marker_panel, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(filenames, "_FeaturePlot.pdf", sep=""), width = 12, height = 10)
  print(FeaturePlot(dataset, features = marker_panel, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(filenames, "_SplitViolins.pdf", sep=""), width = 16, height = 10)
  plots <- VlnPlot(dataset, features = marker_panel, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}


# For all clusters, find the DE genes between conditions
# Note: Now includes GSEA! Does this for each cluster. 
# Example, what's different between tumor vs. egressed in cluster 3? (for each cluster)
# dataset is the seurat object
# objname is obj_name(dataset)
# example function call:
# FindMarkersByConditionEachCluster(combined.DCs, get_obj_name(combined.DCs))
FindMarkersByConditionEachCluster <- function(dataset, objname, runGsea = TRUE) {
  # create subdirectory to put the outputs into
  oldwd <- getwd()
  string1 <- "MarkersByCondition"
  if (runGsea == TRUE) {
    string1 <- "GSEA_&_DE_ByCondition"
  }
  foldername <- paste(string1, objname, sep="_")
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
  
  # I'm gonna assume this is already installed since it's in the readme to do so. 
  ## install organism library for GSEA
  #if (runGsea == TRUE) {
  #  # SET THE DESIRED ORGANISM HERE
  #  # organism = "org.Mm.eg.db"
  #  BiocManager::install(organism, character.only = TRUE)
  #  library(organism, character.only = TRUE)
  #}
  
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
    
    # If we're running the doGSEA function, FindMarkersByCondition will be run anyway. 
    if (runGsea == FALSE) {
      FindMarkersByCondition(cell_subset, objname, ClusterID)
    }
    
    # use a trycatch here since GSEA gives errors so often.
    #doGSEA(cell_subset, objname, ClusterID)
    tryGSEA(cell_subset, objname, ClusterID)
  }
  # Reset to original working directory and default array
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}


# For one cluster, find the DE genes between conditions
# Example, what's different between tumor vs. egressed in cluster 3?
# data should already be subset to only include the cluster of interest. 
# Can define the number of genes to include in the "top" and "bot" comparisons
# If return_table is TRUE, it returns the results table, which can be used for GSEA and other things. 
FindMarkersByCondition <- function(dataset, objname, ClusterID, num_genes = 100, return_table = FALSE) {
  # create subdirectory to put the outputs into
  oldwd <- getwd()
  clustername <- paste("Cluster", ClusterID, sep="") # Cluster0
  gene_amount <- paste("Top", num_genes, sep="") # Top50
  clusterandgene <- paste(clustername, gene_amount, sep="_") # Cluster0_Top50
  foldername <- paste(clustername, objname, sep="_") # Cluster0__Top50_combined.DCs
  filenames <- paste(foldername, Sys.time(), sep="_") # Cluster0__Top50_combined.DCs_datetime...
  dir.create(file.path(filenames))
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "integrated"
  
  
  # make a list of the different conditions, or "orig.ident"s
  cond_list <- unique(dataset$orig.ident)
  # make a new column called celltype to store the old cluster identities
  dataset$celltype <- Idents(dataset)
  # Fill the identities column with the values in the "Condition" column ("egressed" or "tumor")
  Idents(dataset) <- "orig.ident"
  # find markers between the two conditions
  markers <- FindMarkers(dataset, ident.1 = cond_list[1], ident.2 = cond_list[2], min.pct = 0.25, verbose = FALSE)
  # save a text file describing the conditions compared
  write.table(cond_list[1:2], sep="\t", file="Condition1vsCondition2.txt", quote = FALSE, row.names = TRUE)
  # save the results to a text file
  clustername <- paste("Cluster", ClusterID, sep="")
  write.table(markers, sep="\t", file=paste(clustername, "MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(markers['avg_log2FC'], sep="\t", file=paste(clustername, "Markers&LFCByCondition.rnk", sep="_"), quote = FALSE, row.names = TRUE, col.names = FALSE)
  # top 50 most up- and down- regulated genes
  top_genes <- markers %>% top_n(n = num_genes, wt = avg_log2FC)
  bot_genes <- markers %>% top_n(n = (-1)*num_genes, wt = avg_log2FC)
  write.table(top_genes, sep="\t", file=paste(clusterandgene, "UP_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(bot_genes, sep="\t", file=paste(clusterandgene, "DOWN_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(row.names(top_genes), sep="\t", file=paste(clusterandgene, "UP_Names_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(row.names(bot_genes), sep="\t", file=paste(clusterandgene, "DOWN_Names_MarkersByCondition.txt", sep="_"), quote = FALSE, row.names = TRUE)
  
  # Heatmaps of the above markers
  pdf(paste(clusterandgene, "UP_MarkersByCondition.pdf", sep="_"), height = 12)
  print(DoHeatmap(dataset, features = row.names(top_genes), size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 5)))
  dev.off()
  pdf(paste(clusterandgene, "DOWN_MarkersByCondition.pdf", sep="_"), height = 12)
  print(DoHeatmap(dataset, features = row.names(bot_genes), size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 5)))
  dev.off()
  
  # Reset to original working directory and default array
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
  
  # Return the results table if return_table is TRUE
  if (return_table == TRUE) {
    return(markers)
  }
}


# Finds the DE genes between two clusters, outputs tables and heatmaps. 
# Example, what's different between clusters 9 and 1?
# ClusterID1 and ClusterID2 should be strings containing the cluster IDs
# Example call:
# compare_clusters(combined.DCs, get_obj_name(combined.DCs), "0", "1")
compare_clusters <- function(dataset, objname, ClusterID1, ClusterID2) {
  # create subdirectory to put the outputs into
  comparison <- paste(ClusterID1, ClusterID2, sep="_vs_")
  oldwd <- getwd()
  name1 <- paste(comparison, objname, sep="_")
  foldername <- paste(name1, "DEMarkers", sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "integrated"
  
  markers <- FindMarkers(dataset, ident.1 = ClusterID1, ident.2 = ClusterID2, min.pct = 0.25, verbose = FALSE)
  # save the results to a text file
  write.table(markers, sep="\t", file=paste(comparison, "Markers.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(markers['avg_log2FC'], sep="\t", file=paste(comparison, "Markers&LFC.rnk", sep="_"), quote = FALSE, row.names = TRUE, col.names = FALSE)
  # top 50 most up- and down- regulated genes
  top50 <- markers %>% top_n(n = 50, wt = avg_log2FC)
  bot50 <- markers %>% top_n(n = -50, wt = avg_log2FC)
  write.table(top50, sep="\t", file=paste(comparison, "top50UP_Markers.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(bot50, sep="\t", file=paste(comparison, "top50DOWN_Markers.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(row.names(top50), sep="\t", file=paste(comparison, "Names_top50UP_Markers.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(row.names(bot50), sep="\t", file=paste(comparison, "Names_top50DOWN_Markers.txt", sep="_"), quote = FALSE, row.names = TRUE)
  
  cells_compared <- subset(dataset, seurat_clusters == ClusterID1 | seurat_clusters == ClusterID2)
  # Heatmaps of the above markers
  # show all clusters but only the DE genes between the two: 
  pdf(paste(comparison, "top50UP_Markers_AllClusters.pdf", sep="_"), height = 12)
  print(DoHeatmap(dataset, features = row.names(top50), size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 5)))
  dev.off()
  pdf(paste(comparison, "top50DOWN_Markers_AllClusters.pdf", sep="_"), height = 12)
  print(DoHeatmap(dataset, features = row.names(bot50), size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 5)))
  dev.off()
  # show only the two clusters compared
  pdf(paste(comparison, "top50UP_Markers.pdf", sep="_"), height = 12)
  print(DoHeatmap(cells_compared, features = row.names(top50), size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 5)))
  dev.off()
  pdf(paste(comparison, "top50DOWN_Markers.pdf", sep="_"), height = 12)
  print(DoHeatmap(cells_compared, features = row.names(bot50), size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 5)))
  dev.off()
  
  # Reset to original working directory and default array
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
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
  pdf(paste(objname, "DC5_Panel_Violins1.pdf", sep="_"), width = 10, height = 10)
  print(VlnPlot(dataset, features = DC5_panel1, pt.size = 0, ncol = 2))
  dev.off()
  pdf(paste(objname, "DC5_Panel_FeaturePlot1.pdf", sep="_"), width = 10, height = 10)
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


# ------------------------------------------------
# plot the Rgs (1-21) and Grk (1-7) genes
# These play a role in activating/deactivating CXCR4
# CXCR4 is a g-protein coupled receptor
# These can modify it from inside the cell and affect signaling,
# thereby altering cell egress via CXCR4 / CXCL12 activity. 
plot_Rgs_Grk_Markers <- function(dataset, objname) {
  # create subdirectory to put the outputs into
  foldername <- paste("Rgs_Grk_Markers", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  # change default assay to RNA, but create a temporary variable to change it back at the end. 
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # plot UMAPs for context
  pdf(paste(objname, "UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  
  # Define panels
  Rgs_panel1 <- c("Rgs1", "Rgs2", "Rgs3", "Rgs4", "Rgs5", "Rgs6", "Rgs7", "Rgs8", "Rgs9")
  Rgs_panel2 <- c("Rgs10", "Rgs11", "Rgs12", "Rgs13", "Rgs14", "Rgs15")
  Rgs_panel3 <- c("Rgs16", "Rgs17", "Rgs18", "Rgs19", "Rgs20", "Rgs21")
  Grk_panel <- c("Grk1", "Grk2", "Grk3", "Grk4", "Grk5", "Grk6", "Grk7")
  
  # Rgs_panel1
  # Violin plots
  pdf(paste(objname, "Rgs_Violins1.pdf", sep=""), width = 12, height = 10)
  print(VlnPlot(dataset, features = Rgs_panel1, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Rgs_FeaturePlot1.pdf", sep=""), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Rgs_panel1, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Rgs_SplitViolins1.pdf", sep=""), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Rgs_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Rgs_panel2
  # Violin plots
  pdf(paste(objname, "Rgs_Violins2.pdf", sep=""), width = 12, height = 10)
  print(VlnPlot(dataset, features = Rgs_panel2, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Rgs_FeaturePlot2.pdf", sep=""), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Rgs_panel2, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Rgs_SplitViolins2.pdf", sep=""), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Rgs_panel2, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Rgs_panel3
  # Violin plots
  pdf(paste(objname, "Rgs_Violins3.pdf", sep=""), width = 12, height = 10)
  print(VlnPlot(dataset, features = Rgs_panel3, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Rgs_FeaturePlot3.pdf", sep=""), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Rgs_panel3, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Rgs_SplitViolins3.pdf", sep=""), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Rgs_panel3, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Grk_panel
  # Violin plots
  pdf(paste(objname, "Grk_Violins1.pdf", sep=""), width = 12, height = 10)
  print(VlnPlot(dataset, features = Grk_panel, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Grk_FeaturePlot1.pdf", sep=""), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Grk_panel, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Grk_SplitViolins1.pdf", sep=""), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Grk_panel, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}


# ------------------------------------------------
# explore Egress Markers for Maria's paper
# Exploring and defining the CD8+ clusters ------------------------------------
# CXCR4, its ligand CXCL12, the decoy receptor Ackr3
# Intracellular regulators of CXCR4, the Rgs genes
# Cell cycle/proliferation marker Mki67
# dataset is the dataset, objname is the name of the seurat object
explore_EgressMarkers <- function(dataset, objname) {
  # create subdirectory to put the outputs into
  foldername <- paste("ExploreEgressMarkers", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # plot UMAPs for context
  pdf(paste(objname, "UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  
  # Define panels
  Egress_panel1 <- c("Cxcr4", "Cxcl12", "Ackr3", "Rgs1", "Rgs2", "Rgs3", "Rgs10", "Rgs16")
  Egress_panel2 <- c("Ccr7", "S1pr1", "Cd69", "", "", "", "", "")
  Exhaustion_panel1 <- c("Pdcd1", "Havcr2", "Lag3", "Tox", "Mki67", "H2-Ab1")
  # Thymic_Egress <- c("Sphk1", "Sphk2", "Spns2", "Sp1", "Lpp3", "LTBR", "MST1", "MST2", "Coro1a")
  Thymic_Egress <- c("Sphk1", "Sphk2", "Spns2", "Sp1", "LPP3", "Ltbr", "Mst1", "Stk3", "Coro1a")
  # Stk3 is Mst2
  # Egress_modulators are from James, Jenkinson, and Anderson 2018, it's thymic egress but still checking
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6174998/
  # Taylor CD8+ T-cell subtype panels -------------------------------------------
  # Vps37b, Icos, Tcf7, S1pr1, Cxcr6, Itgae, Jun, Fos, Nr4a3, Klf2, Ifng, Cd69
  # Gzmb and Sell
  # Eomes, Tbx21, Runx3, Id2, Id3, Tox, Il7r, Klrg1, Ccr7, Pdcd1, Havcr2, and Lag3
  
  # Egress Panel1
  # Violin plots
  pdf(paste(objname, "Egress_Violins1.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Egress_panel1, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Egress_FeaturePlot1.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Egress_panel1, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Egress_SplitViolins1.pdf", sep="_"), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Egress_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Egress Panel1
  # Violin plots
  pdf(paste(objname, "Egress_Violins2.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Egress_panel2, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Egress_FeaturePlot2.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Egress_panel2, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Egress_SplitViolins2.pdf", sep="_"), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Egress_panel2, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Exhaustion Panel1
  # Violin plots
  pdf(paste(objname, "Exhaustion_Violins1.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Exhaustion_panel1, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Exhaustion_FeaturePlot1.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Exhaustion_panel1, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Exhaution_SplitViolins1.pdf", sep="_"), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Exhaustion_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Thymic Egress Panel1
  # Violin plots
  pdf(paste(objname, "ThymicEgress_Violins1.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Thymic_Egress, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "ThymicEgress_FeaturePlot1.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Thymic_Egress, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "ThymicEgress_SplitViolins1.pdf", sep="_"), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Thymic_Egress, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)

}

# ------------------------------------------------
# plot Egress Markers for Maria's paper
# Making more legit plots with select markers
# Exploring and defining the CD8+ clusters ------------------------------------
# CXCR4, its ligand CXCL12, the decoy receptor Ackr3
# Intracellular regulators of CXCR4, the Rgs genes
# Cell cycle/proliferation marker Mki67
# dataset is the dataset, objname is the name of the seurat object
plot_EgressMarkers <- function(dataset, objname) {
  # create subdirectory to put the outputs into
  foldername <- paste("PlotEgressMarkers", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # plot UMAPs for context
  pdf(paste(objname, "UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  
  # Define panels
  Main_panel1 <- c("S1pr1", "Rgs3", "Pdcd1", "Rgs1", "Rgs2", "Havcr2", "Rgs10", "Rgs16", "Lag3")
  Egress_panel1 <- c("Cxcr4", "Cxcl12", "Ackr3", "Ccr7", "S1pr1", "Sp1", "Rgs1", "Rgs2", "Rgs3", "Rgs10", "Rgs16")
  Exhaustion_panel1 <- c("Pdcd1", "Havcr2", "Lag3", "Tox", "Mki67", "H2-Ab1")
  Taylor_panel1 <- c("Tcf7", "Cxcr6", "Ifng", "Gzmb", "Sell", "Eomes", "Il7r", "Tbx21", "Id2")
  # Taylor CD8+ T-cell subtype panels -------------------------------------------
  # Vps37b, Icos, Tcf7, S1pr1, Cxcr6, Itgae, Jun, Fos, Nr4a3, Klf2, Ifng, Cd69
  # Gzmb and Sell
  # Eomes, Tbx21, Runx3, Id2, Id3, Tox, Il7r, Klrg1, Ccr7, Pdcd1, Havcr2, and Lag3
  
  # Main Panel1
  # Violin plots
  pdf(paste(objname, "Main_Violins1.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Main_panel1, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Main_FeaturePlot1.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Main_panel1, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Main_SplitViolins1.pdf", sep="_"), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Main_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Egress Panel1
  # Violin plots
  pdf(paste(objname, "Egress_Violins1.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Egress_panel1, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Egress_FeaturePlot1.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Egress_panel1, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Egress_SplitViolins1.pdf", sep="_"), width = 21, height = 13)
  plots <- VlnPlot(dataset, features = Egress_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 4))
  dev.off()
  
  # Exhaustion Panel1
  # Violin plots
  pdf(paste(objname, "Exhaustion_Violins1.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Exhaustion_panel1, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Exhaustion_FeaturePlot1.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Exhaustion_panel1, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Exhaution_SplitViolins1.pdf", sep="_"), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Exhaustion_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  # Taylor Panel1
  # Violin plots
  pdf(paste(objname, "Taylor_Violins1.pdf", sep="_"), width = 12, height = 10)
  print(VlnPlot(dataset, features = Taylor_panel1, pt.size = 0, ncol = 3))
  dev.off()
  # Feature plots
  pdf(paste(objname, "Taylor_FeaturePlot1.pdf", sep="_"), width = 12, height = 10)
  print(FeaturePlot(dataset, features = Taylor_panel1, ncol = 3))
  dev.off()
  # Violin plots split by condition
  pdf(paste(objname, "Taylor_SplitViolins1.pdf", sep="_"), width = 20, height = 10)
  plots <- VlnPlot(dataset, features = Taylor_panel1, split.by = "orig.ident", 
                   pt.size = 0, combine = FALSE)
  print(wrap_plots(plots = plots, nrow = 3))
  dev.off()
  
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
  
}

# doPseudotime
doPseudotime <- function(dataset, objname, partition_option = FALSE) {
  # create subdirectory to put the outputs into
  foldername <- paste("doPseudotime", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  # temporarily change the default assay during the use of this function
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # Building trajectories with Monocle3
  cds <- as.cell_data_set(dataset)
  cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
  # if the user defines partition_option = TRUE, it will separate into the two partition (as long as they're there)
  cds <- learn_graph(cds, use_partition = partition_option)

  ## Calculate size factors using built-in function in monocle3
  cds <- estimate_size_factors(cds)
  
  ## Add gene names into CDS
  cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(dataset[["RNA"]])
  
  # This step opens a pop-up window for user input. 
  # The user has to pick a starting point for pseudotime. 
  # order cells
  cds <- order_cells(cds, reduction_method = "UMAP")
  
  # plot trajectories colored by pseudotime
  pseudotime_plt <- plot_cells(
    cds = cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    cell_size = 1.25,
    group_label_size = 4,
    graph_label_size = 8,
    trajectory_graph_segment_size = 1.75
  )
  
  # Plot UMAPs for context
  umap_plt <- DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE)
  umap_by_cond <- DimPlot(dataset, reduction = "umap", group.by = "orig.ident")
  
  pdf(paste(objname, "UMAPs.pdf", sep="_"), width=14)
  print(umap_plt + labs(title = "Clusters") + theme(plot.title = element_text(hjust = 0.5)) + 
    umap_by_cond + labs(title = "By Condition"))
  dev.off()
  
  pdf(paste(objname, "UMAP.pdf", sep="_"), width=10, height = 10)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE))
  dev.off()
  
  # Plot pseudotime by itself
  pdf(paste(objname, "pseudotime.pdf", sep="_"), width=12, height=12)
  print(pseudotime_plt + labs(title = "Pseudotime") + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  # try a complex plotting layout with pseudotime, general UMAP, and UMAP split by condition: 
  lay1 <- rbind(c(1,1,2),
               c(1,1,3))
  
  pdf(paste(objname, "pseudotime+UMAPs.pdf", sep="_"), width=21, height=14)
  grid.arrange(pseudotime_plt + labs(title = "Pseudotime") + theme(plot.title = element_text(hjust = 0.5)), 
               umap_plt + labs(title = "Clusters") + theme(plot.title = element_text(hjust = 0.5)),  
               umap_by_cond + labs(title = "By Condition"), 
               layout_matrix = lay1)
  dev.off()
  
  # Need to add the rest of the Monocle3 functionality next once I learn it. 
  
  # return to defaults
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd)
}
  
# doGSEA
# For mouse, organism is "org.Mm.eg.db"
# numsets is the number of GSEA plots to make (for the top up and down most enriched sets)
# num_genesets is the number of gene sets to include in the table output. 
doGSEA <- function(dataset, objname, ClusterID = "All", organism = "org.Mm.eg.db", numsets = 20, num_genesets = 200, installOrganism = FALSE) {
  # GSEA with ClusterProfiler
  # https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
  # create subdirectory to put the outputs into
  clustername <- paste("Cluster", ClusterID, sep="") # Cluster0
  folderstring <- paste(clustername, "doGSEA", sep = "_") # Cluster0_doGSEA
  foldername <- paste(folderstring, objname, sep="_") # Cluster0_doGSEA_combined.CD8s
  filenames <- paste(foldername, Sys.time(), sep="_") # Cluster0_doGSEA_combined.CD8s_datetime
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  # temporarily change the default assay during the use of this function
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # make a list of the different conditions, or "orig.ident"s
  cond_list <- unique(dataset$orig.ident)
  # make a string storing the comparison made, to include in file names later
  clustername <- paste("Cluster", ClusterID, sep="") # Cluster0
  comparison_str <- paste(cond_list[1], cond_list[2], sep = "_vs_") # Egress_vs_Tumor
  comparison_str2 <- paste(comparison_str, objname, sep = "_") # Egress_vs_Tumor_combined.CD8s
  filename_str <- paste(clustername, comparison_str2, sep = "_") # Cluster0_Egress_vs_Tumor_combined.CD8s
  
  #if (installOrganism == TRUE) {
  #  # SET THE DESIRED ORGANISM HERE
  #  # organism = "org.Mm.eg.db"
  #  BiocManager::install(organism, character.only = TRUE)
  #  library(organism, character.only = TRUE)
  #}
  
  # get the DE gene table using my custom function (essentially uses seurat's FindMarkers)
  DE.genes.by.condition <- FindMarkersByCondition(dataset, objname, ClusterID, num_genes = 100, return_table = TRUE)
  # get the log2 fold change values
  ranked.list.prep <- DE.genes.by.condition$avg_log2FC
  # reconnect with the gene names
  names(ranked.list.prep) <- row.names(DE.genes.by.condition)
  # omit any NA values
  ranked.list <- na.omit(ranked.list.prep)
  # sort the list in decreasing order (required for clusterProfiler)
  ranked.list = sort(ranked.list, decreasing = TRUE)
  
  # check which options are available...
  # keytypes(org.Mm.eg.db)
  # For Lund Lab 2021, our genes are in "SYMBOL"
  
  # Run GSEA... 
  # ont - one of “BP”, “MF”, “CC” or “ALL”
  # nPerm - the higher the number of permutations you set, the more accurate your result will, but the longer the analysis will take.
  # note: I got errors when using nPerm. They say it's not recommended to use it anymore. So I took it out. 
  # used to be nPerm = 10000
  # minGSSize - minimum number of genes in set (gene sets with lower than this many genes in your dataset will be ignored).
  # maxGSSize - maximum number of genes in set (gene sets with greater than this many genes in your dataset will be ignored).
  # pvalueCutoff - pvalue Cutoff.
  # pAdjustMethod - one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
  
  gse <- gseGO(geneList=ranked.list, 
               ont ="ALL", 
               keyType = "SYMBOL", 
               minGSSize = 5, 
               maxGSSize = 1000, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "none")
  
  # outputs
  # dotplot
  require(DOSE)
  pdf(paste(filename_str, "GSEA_dotplot.pdf", sep="_"), height=16, width=14)
  print(dotplot(gse, showCategory=30, split=".sign") + facet_grid(.~.sign))
  dev.off()
  
  # category Netplot
  # categorySize can be either 'pvalue' or 'geneNum'
  pdf(paste(filename_str, "GSEA_netplot.pdf", sep="_"), height=16, width=14)
  print(cnetplot(gse, categorySize="pvalue", foldChange=ranked.list, showCategory = 3))
  dev.off()
  
  # ridgeplot
  pdf(paste(filename_str, "GSEA_ridgeplot.pdf", sep="_"), height=16, width=14)
  print(ridgeplot(gse) + labs(x = "enrichment distribution"))
  dev.off()
  
  # Filter the data to save tables and make GSEA plots for the highest NES & p-values
  # create a temp variable to experiment with
  gse_sorted <- gse
  # order it on NES & p-value
  hi_ndx <- order(-gse_sorted$NES, gse_sorted$pvalue)
  lo_ndx <- order(gse_sorted$NES, gse_sorted$pvalue)
  gse_sorted_enrichment_hi <- gse_sorted[hi_ndx,]
  gse_sorted_enrichment_lo <- gse_sorted[lo_ndx,]
  # Filter out genesets with 5 or fewer genes
  # did this in the gseGO() call so commenting these out. 
  # gse_sorted_enrichment_hi <- subset(gse_sorted_enrichment_hi, gse_sorted_enrichment_hi$setSize > 5)
  # gse_sorted_enrichment_lo <- subset(gse_sorted_enrichment_lo, gse_sorted_enrichment_lo$setSize > 5)
  # these should still be sorted, so good to go for next steps. I think. 
  
  # Save txt files of the sorted data to look at
  # top 50 most up- and down- regulated genes
  # Only include ones with p-values less than 0.1
  up_genesets <- subset(gse_sorted_enrichment_hi, pvalue <= 0.1)
  down_genesets <- subset(gse_sorted_enrichment_lo, pvalue <= 0.1)
  write.table(up_genesets, sep="\t", file=paste(filename_str, "UP_GSEA.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(down_genesets, sep="\t", file=paste(filename_str, "DOWN_GSEA.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(up_genesets$description, sep="\t", file=paste(filename_str, "UP_Names_GSEA.txt", sep="_"), quote = FALSE, row.names = TRUE)
  write.table(down_genesets$description, sep="\t", file=paste(filename_str, "DOWN_Names_GSEA.txt", sep="_"), quote = FALSE, row.names = TRUE)
  
  # GSEA Plots
  # UP-regulated genesets
  # create subdirectory to put the GSEA plots into
  # clustername # Cluster0
  
  gsea_folder <- paste(filenames,"upGSEAPlots", sep = "_") # Cluster0_doGSEA_combined.CD8s_datetime
  dir.create(file.path(gsea_folder))
  maindir <- getwd()
  setwd(file.path(gsea_folder))
  
  # cycle through the top N=20 and plot
  ct <- 1
  # cycles through every geneset ID in the sorted dataframe
  for (geneset_id in gse_sorted_enrichment_hi$ID[1:numsets]) {
    # reformat file naming for each cycle
    gsea_filename_str <- paste(filename_str, ct, sep="_")
    # plot and save gsea enrichment plot for up-regulated genes
    # Create an annotation for the NES and p-value with 5 sig figs
    grob1 <- grobTree(textGrob(paste("NES: ", trimws(format(round(gse_sorted_enrichment_hi$NES[ct],5), nsmall=5))), x=0.7,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=13)))
    grob2 <- grobTree(textGrob(paste("p-value: ", trimws(format(round(gse_sorted_enrichment_hi$pvalue[ct],5), nsmall=5))), x=0.7,  y=0.90, hjust=0,
                              gp=gpar(col="red", fontsize=13)))
    gplt <- enrichplot::gseaplot(gse, geneSetID = geneset_id, by = "runningScore", title = gse_sorted_enrichment_hi$Description[ct])
    pdf(paste(gsea_filename_str, "GSEA_UP.pdf", sep="_"))
    #print(enrichplot::gseaplot(gse, geneSetID = geneset_id, by = "runningScore", title = gse_sorted_enrichment_hi$Description[ct]))
    print(gplt + annotation_custom(grob1) + annotation_custom(grob2))
    dev.off()
    # plot and save gsea enrichment plot for down-regulated genes
    #pdf(paste(gsea_filename_str, "GSEA_DOWN.pdf", sep="_"))
    # print(gseaplot(gse_sorted_enrichment_lo, by = "all", title = gse_sorted_enrichment_lo$Description[ct], geneSetID = ct))
    # gseaplot(gse_sorted_enrichment_lo, by = "all", title = gse_sorted_enrichment_lo$Description[ct], geneSetID = ct)
    #enrichplot::gseaplot(gse_sorted_enrichment_lo, geneSetID = ct, by = "runningScore", title = gse_sorted_enrichment_lo$Description[ct])
    #dev.off()
    ct <- ct + 1
  }
  
  # DOWN-retulated genesets...
  # return to main directory
  setwd(maindir) # come out of gsea plots folder
  # DOWN-regulated genes
  # create subdirectory to put the GSEA plots into
  gsea_folder <- paste(filenames,"downGSEAPlots", sep = "_")
  dir.create(file.path(gsea_folder))
  setwd(file.path(gsea_folder))
  
  # cycle through the top N=20 and plot
  ct <- 1
  # cycles through every geneset ID in the sorted dataframe
  for (geneset_id in gse_sorted_enrichment_lo$ID[1:numsets]) {
    # reformat file naming for each cycle
    gsea_filename_str <- paste(filename_str, ct, sep="_")
    # plot and save gsea enrichment plot for up-regulated genes
    # Create an annotation for the NES and p-value
    # format(round(x, k), nsmall = k) # cuts x to only show k decimal places
    grob1 <- grobTree(textGrob(paste("NES: ", trimws(format(round(gse_sorted_enrichment_lo$NES[ct],3), nsmall=3))), x=0.7,  y=0.95, hjust=0,
                               gp=gpar(col="red", fontsize=13)))
    grob2 <- grobTree(textGrob(paste("p-value: ", trimws(format(round(gse_sorted_enrichment_lo$pvalue[ct],3), nsmall=3))), x=0.7,  y=0.90, hjust=0,
                               gp=gpar(col="red", fontsize=13)))
    gplt <- enrichplot::gseaplot(gse, geneSetID = geneset_id, by = "runningScore", title = gse_sorted_enrichment_lo$Description[ct])
    pdf(paste(gsea_filename_str, "GSEA_DOWN.pdf", sep="_"))
    # print(enrichplot::gseaplot(gse, geneSetID = geneset_id, by = "runningScore", title = gse_sorted_enrichment_lo$Description[ct]))
    print(gplt + annotation_custom(grob1) + annotation_custom(grob2))
    dev.off()
    # plot and save gsea enrichment plot for down-regulated genes
    #pdf(paste(gsea_filename_str, "GSEA_DOWN.pdf", sep="_"))
    # print(gseaplot(gse_sorted_enrichment_lo, by = "all", title = gse_sorted_enrichment_lo$Description[ct], geneSetID = ct))
    # gseaplot(gse_sorted_enrichment_lo, by = "all", title = gse_sorted_enrichment_lo$Description[ct], geneSetID = ct)
    #enrichplot::gseaplot(gse_sorted_enrichment_lo, geneSetID = ct, by = "runningScore", title = gse_sorted_enrichment_lo$Description[ct])
    #dev.off()
    ct <- ct + 1
  }
  
  
  # return to defaults
  setwd(maindir) # come out of gsea plots folder
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd) # come out of the folder for this process
}

# doAnalysis
# This will run most of the other functions here.
# 
doAnalysis <- function(dataset, objname, annotation_plots = TRUE) {
  
  # create subdirectory to put the outputs into
  foldername <- paste("SC-Analysis", objname, sep="_")
  filenames <- paste(foldername, Sys.time(), sep="_")
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  # temporarily change the default assay during the use of this function
  tempdftassay <- DefaultAssay(dataset)
  DefaultAssay(dataset) <- "RNA"
  
  # Run Functions --------------------------------------------------------------
  # For the entire dataset, runs pseudotime analysis. 
  # requires user input for selecting the starting point. 
  doPseudotime(dataset, objname)# For each cluster, finds DE genes between conditions. 
  
  # Now includes GSEA comparing conditions for each cluster
  FindMarkersByConditionEachCluster(dataset, objname)
  
  # note: decided to skip this for now... not that worth it. 
  # For the entire dataset, find DE genes between conditions and find enriched gene sets
  # requires user input for installing the species library
  # doGSEA(dataset, objname)
  
  # Call all the annotation plotting functions. 
  if (annotation_plots == TRUE) {
    plot_huang2021(dataset, objname)
    
    plot_FilioDCMarkers(dataset, objname)
    
    plot_TaylorCD8Markers(dataset, objname)
    
    plot_Rgs_Grk_Markers(dataset, objname)
    
    plot_EgressMarkers(dataset, objname)
  }
  # ----------------------------------------------------------------------------
  
  # return to defaults
  DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd) # come out of the folder for this process
}


# Subset out the desired cluster first, then call this function to do clustering. 
# Example:
# # Subset CD8+ T-cells...
# combined.cd8s <- subset(immune.combined, idents = c("2: CD8+ T-cells", "5: Mitotic T-cells"))
# combined.cd8s <- doSeuratClustering(combined.cd8s)
doSeuratClustering <- function(dataset, objname, num_pca_dims=50, evaluatejackstraw=TRUE) {
  # Make subdirectory for this clustering
  foldername <- paste(objname, "clustering", sep="_") # combined.CD8s_clustering
  filenames <- paste(foldername, Sys.time(), sep="_") # combined.CD8s_clustering_datetime
  dir.create(file.path(filenames))
  oldwd <- getwd()
  setwd(file.path(filenames))
  
  # Begin Clustering
  print('~~FindVariableFeatures~~')
  dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000)
  print('~~ScaleData~~')
  dataset <- ScaleData(dataset, verbose = FALSE)
  print('~~RunPCA~~')
  dataset <- RunPCA(dataset, npcs = num_pca_dims, verbose = FALSE)
  
  if (evaluatejackstraw == TRUE) {
    print('~~RunJackStraw~~')
    # Which PCAs to select? Save these plots for later just to see. 
    dataset <- JackStraw(dataset, num.replicate = 100, dims = num_pca_dims)
    dataset <- ScoreJackStraw(dataset, dims = 1:num_pca_dims)
    pdf(paste(objname, "Check_Jackstar_dim50.pdf"), width = 14)
    print(JackStrawPlot(dataset, dims = 1:num_pca_dims))
    dev.off()
    pdf(paste(objname, "Elbow_PCA_dim50.pdf"))
    print(ElbowPlot(dataset))
    dev.off()
  }
  
  # Finish clustering
  print('~~RunUMAP~~')
  dataset <- RunUMAP(dataset, reduction = "pca", dims = 1:num_pca_dims)
  print('~~FindNeighbors~~')
  dataset <- FindNeighbors(dataset, reduction = "pca", dims = 1:num_pca_dims)
  print('~~FindClusters~~')
  dataset <- FindClusters(dataset, resolution = 0.5)
  
  # temporarily change the default assay during the use of this function
  # tempdftassay <- DefaultAssay(dataset)
  # DefaultAssay(dataset) <- "RNA"
  
  print('~~Saving UMAPs~~')
  # Visualize and save UMAPs
  pdf(paste(objname, "UMAP1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  # split by orig.ident
  pdf(paste(objname, "UMAP_byCondition1.pdf", sep="_"), width=12, height = 12)
  print(DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6))
  dev.off()
  # Both plots together
  p1 <- DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6)
  p2 <- DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.75, label.size = 6)
  pdf(paste(objname, "UMAPs.pdf", sep="_"), width=20, height = 12)
  print(p1 + p2)
  dev.off()
  
  # return to defaults
  # DefaultAssay(dataset) <- tempdftassay
  setwd(oldwd) # come out of the folder for this process
  
  endmsg <- paste("~~ Finished with ", objname, sep = "")
  print(paste(endmsg, "~~", sep = ""))
  return(dataset)
}

# a GSEA try catch function to continue with the code if GSEA fails, which it
# often does depending on the results of the DE gene analysis... If there are very few
# DE expressed genes, for example. 
tryGSEA <- function(dataset, objname, ClusterID) {
  tryCatch(
    
    ########################################################
    # Try part: define the expression(s) you want to "try" #
    ########################################################
    
    {
      # Try...
      message(paste("Trying GSEA for cluster ", ClusterID))
      oldwd <- getwd() # save the current wd. If it errors out, reset to that. 
      # Write some stats to a log
      cond_list <- unique(dataset$orig.ident)
      
      tryGSEA_log <- paste(oldwd, "tryGSEA_log.txt", sep = "/")
      cat("-----------------------------------------------------------", file = tryGSEA_log, append = TRUE, sep = "\n")
      cat("-----------------------------------------------------------", file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(paste("Cluster ", ClusterID), file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(paste("dataset size: ", as.character(dim(dataset))), file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(paste("number of cond1: ", as.character(dim(dataset[ which(dataset$orig.ident==cond_list[1])])[1])), file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(paste("number of cond2: ", as.character(dim(dataset[ which(dataset$orig.ident==cond_list[2])])[1])), file = tryGSEA_log, append = TRUE, sep = "\n")
      
      #tryGSEA_log <- file(paste(oldwd, "tryGSEA_log.txt", sep = "/"), open = "wt")
      #writeLines(c("\n",
      #             "-----------------------------------------------------------",
      #             "-----------------------------------------------------------",
      #             paste("Cluster ", ClusterID), 
      #             paste("dataset size: ", as.character(dim(dataset))),
      #             paste("number of cond1: ", as.character(dim(dataset[ which(dataset$orig.ident==cond_list[1])]))),
      #             paste("number of cond2: ", as.character(dim(dataset[ which(dataset$orig.ident==cond_list[2])])))),
      #           con = tryGSEA_log)
      #close(tryGSEA_log)
      
      doGSEA(dataset, objname, ClusterID)
    },
    
    ########################################################################
    # Condition handler part: define how you want conditions to be handled #
    ########################################################################
    
    # Handler when a warning occurs:
    warning = function(cond) {
      tryGSEA_log <- paste(oldwd, "tryGSEA_log.txt", sep = "/")
      cat("- - - - - - - - - - - - - - - - - - - -", file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(paste(" WARNING IN CLUSTER ", ClusterID), file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(as.character(cond), file = tryGSEA_log, append = TRUE, sep = "\n")
      
      #tryGSEA_log <- file(paste(oldwd, "tryGSEA_log.txt", sep = "/"), open = "wt")
      #writeLines(c("\n",
      #             "- - - - - - - - - - - - - - - - - - - -",
      #             paste(" WARNING IN CLUSTER ", ClusterID), 
      #             "Warning message: ",
      #             as.character(cond)),
      #           con = tryGSEA_log)
      #close(tryGSEA_log)
      
      message(paste("A warning was caused for cluster ", ClusterID))
      message("Here's the original warning message: ")
      message(cond)
    },
    
    # Handler when an error occurs:
    error = function(cond) {
      tryGSEA_log <- paste(oldwd, "tryGSEA_log.txt", sep = "/")
      cat("- - - - - - - - - - - - - - - - - - - -", file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(paste(" ERROR IN CLUSTER ", ClusterID), file = tryGSEA_log, append = TRUE, sep = "\n")
      cat(as.character(cond), file = tryGSEA_log, append = TRUE, sep = "\n")
      
      #tryGSEA_log <- file(paste(oldwd, "tryGSEA_log.txt", sep = "/"), open = "wt")
      #writeLines(c("\n",
      #             "- - - - - - - - - - - - - - - - - - - -",
      #             paste(" ERROR IN CLUSTER ", ClusterID), 
      #             "Error message: ",
      #             as.character(cond)),
      #           con = tryGSEA_log)
      #close(tryGSEA_log)
      
      message(paste("An error was caused for cluster ", ClusterID))
      message("Here's the original error message: ")
      message(cond)
    },
    
    ###############################################
    # Final part: define what should happen AFTER #
    # everything has been tried and/or handled    #
    ###############################################
    
    finally = {
      setwd(oldwd) # come out of the folder for this process
      message(paste("Processed GSEA for cluster ", ClusterID))
    }
  )    
}
