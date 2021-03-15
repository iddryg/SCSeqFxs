## SCSeqFxs
Some functions for streamlined single cell sequencing analysis. Uses seurat objects and functions. 
Written by Ian Dryg for use in Amanda Lund's lab (www.thelundlab.com) in March 2021. 

## Documentation
# To use these functions:
1. Download SCSeqFxs.R and move to whatever folder you want. 
2. In your own R script, you have to source SCSeqFxs.R to be able to use the functions included. To do this, run:
souce('filepath...../SCSeqFxs.R')
example: 
source('/Users/drygi01/Documents/R/SCSeqFxs.R')
3. Then you can use any of the functions listed below. 

# Dependencies: 
Requires the following libraries:
Seurat

# List of functions:
1. get_obj_name(x)
2. plot_markers(dataset, objname, marker_panel)
3. FindMarkersByConditionEachCluster(dataset, objname)
4. FindMarkersByCondition(dataset, objname, ClusterID)
5. plot_huang2021(dataset, objname)
6. plot_FilioDCMarkers(datset, objname)
7. plot_TaylorCD8Markers(dataset, objname)

## Function Descriptions
# get_obj_name(x)
Gets the name of an object. You can use this to provide the "objname" parameter for most of the other functions. 
Example:
plot_markers(dataset, get_obj_name(dataset), marker_panel)

# plot_markers(dataset, objname, marker_panel)
Makes some basic plots from a seurat object. You need to provide the object name and a list of strings with the gene names to plot. 
Plots will be saved at pdf files. You can edit the R file to change the size of the pdf files if you wish. 
Plots created: 
- Violin plots showing gene expression in each cluster (seurat VlnPlot() function)
- Feature Plot showing gene expression overlaid on the UMAP clustering (seurat FeaturePlot() function)
- Violin plots showing gene expression in each cluster, split by condition
Example call:
plot_markers(combined.DCs, get_obj_name(combined.DCs), c("Irf4", "Il12a", "Il15", "Cxcl12", "Cxcl16", "Ccl19", "Mki67", "Xcr1", "Eadem1"))

# FindMarkersByConditionEachCluster(dataset, objname)
Finds the top differentially expressed genes between conditions within a given cluster, for every cluster in the dataset. 
For example, "What's differentially expressed between tumor vs. egressed populations in cluster 3? (for each cluster)"
For each cluster, it will save text files with the top 50 up- and down- regulated genes between conditions. Additionally, ALL DE genes will be saved as a text file. Also, a .rnk file will be saved with only the top 50 up- and down- regulated genes and their L2FC values to use for GSEA. 
This function uses FindMarkersByCondition() within it. 
Example call: 
FindMarkersByConditionEachCluster(combined.DCs, get_obj_name(combined.DCs))

# FindMarkersByCondition(dataset, objname, ClusterID)
Finds the top differentially expressed genes between conditions within a given cluster. 
For example, "What's differentially expressed between tumor vs. egressed populations in cluster 3?"
It will save text files with the top 50 up- and down- regulated genes between conditions. Additionally, ALL DE genes will be saved as a text file. Also, a .rnk file will be saved with only the top 50 up- and down- regulated genes and their L2FC values to use for GSEA. 
ClusterID is assumed to be the value in seurat_clusters, but is only used for file names. 
This function is used by FindMarkersByCondtitionEachCluster(). 

# plot_huang2021(dataset, objname)
Makes a bunch of plots with all of the markers and subpopulations defined in Huang et al 2021, Figure 6 (innervated lymph node paper). 
Saves the plots in a subdirectory of your current working directory. 
https://www.sciencedirect.com/science/article/pii/S0092867420315646?via%3Dihub
Makes violin plots, feature plots, and violin plots split by condition. (same as plot_markers() above)

# plot_FilioDCMarkers(datset, objname)
Makes a bunch of plots with various DC markers that Filio sent me in 2021. 
Saves the plots in a subdirectory of your current working directory. 
Makes violin plots, feature plots, and violin plots split by condition. (same as plot_markers() above)

# plot_TaylorCD8Markers(datset, objname)
Makes a bunch of plots with various markers of interest for CD8 cells that Taylor sent me in 2021. 
Saves the plots in a subdirectory of your current working directory. 
Makes violin plots, feature plots, and violin plots split by condition. (same as plot_markers() above)


