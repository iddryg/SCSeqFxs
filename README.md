# SCSeqFxs
Some functions for streamlined single cell sequencing analysis. Uses seurat objects and functions. 
Written by Ian Dryg for use in Amanda Lund's lab (www.thelundlab.com) in March 2021. 

# Documentation
## To use these functions:
1. Download SCSeqFxs.R and move to whatever folder you want. 
2. In your own R script, you have to source SCSeqFxs.R to be able to use the functions included. To do this, run:
source('filepath...../SCSeqFxs.R')
example: 
source('/Users/drygi01/Documents/R/SCSeqFxs.R')
3. Then you can use any of the functions listed below. 

## Dependencies: 
R Studio recommended (https://www.rstudio.com/)
Requires the following libraries:
- Seurat: https://satijalab.org/seurat/articles/install.html
- cowplot: https://cran.r-project.org/web/packages/cowplot/readme/README.html
- tidyverse (includes ggplot2, dplyr, others): https://ggplot2.tidyverse.org/
- patchwork: https://cran.r-project.org/web/packages/patchwork/readme/README.html
- dittoSeq: http://www.bioconductor.org/packages/release/bioc/html/dittoSeq.html
- clusterProfiler: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
- enrichplot: http://bioconductor.org/packages/release/bioc/html/enrichplot.html
- DOSE: https://bioconductor.org/packages/release/bioc/html/DOSE.html
- org.Mm.eg.db: https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html

## List of functions:
1. doAnalysis(dataset, objname)
2. get_obj_name(x)
3. plot_markers(dataset, objname, marker_panel)
4. FindMarkersByConditionEachCluster(dataset, objname)
5. FindMarkersByCondition(dataset, objname, ClusterID)
6. compare_clusters(datset, objname, ClusterID1, ClusterID2)
7. plot_huang2021(dataset, objname)
8. plot_FilioDCMarkers(datset, objname)
9. plot_TaylorCD8Markers(dataset, objname)
10. doPseudotime(dataset, objname)
11. doGSEA(dataset, objname)

# Function Descriptions
## doAnalysis(dataset, objname)
Runs several of the functions below to streamline cluster annotation and a general single cell analysis. 
- doPseudotime 
- FindmarkersByConditionEachCluster (DE gene analysis and GSEA by condition for each cluster)
- doGSEA (DE gene analysis and GSEA for the entire dataset together)
- plot_huang2021
- plot_FilioDCMarkers
- plot_TaylorCD8Markers
- plot_Rgs_Grk_Markers
- plot_EgressMarkers

## get_obj_name(x)
Gets the name of an object. You can use this to provide the "objname" parameter for most of the other functions. 
Example:
plot_markers(dataset, get_obj_name(dataset), marker_panel)

## plot_markers(dataset, objname, marker_panel)
Makes some basic plots from a seurat object. You need to provide the object name and a list of strings with the gene names to plot. 
Plots will be saved at pdf files. You can edit the R file to change the size of the pdf files if you wish. 
Plots created: 
- Violin plots showing gene expression in each cluster (seurat VlnPlot() function)
- Feature Plot showing gene expression overlaid on the UMAP clustering (seurat FeaturePlot() function)
- Violin plots showing gene expression in each cluster, split by condition
Example call:
plot_markers(combined.DCs, get_obj_name(combined.DCs), c("Irf4", "Il12a", "Il15", "Cxcl12", "Cxcl16", "Ccl19", "Mki67", "Xcr1", "Eadem1"))

## FindMarkersByConditionEachCluster(dataset, objname)
Finds the top differentially expressed genes between conditions within a given cluster, for every cluster in the dataset. Also, performs GSEA using these differentially expressed genes using the doGSEA() function below. 
For example, "What's differentially expressed between tumor vs. egressed populations in cluster 3? (for each cluster)"
For each cluster, it will save text files with the top 50 up- and down- regulated genes between conditions. Additionally, ALL DE genes will be saved as a text file. Also, a .rnk file will be saved with only the top 50 up- and down- regulated genes and their L2FC values to use for GSEA. 
This function uses FindMarkersByCondition() within it. 
Example call: 
FindMarkersByConditionEachCluster(combined.DCs, get_obj_name(combined.DCs))

## FindMarkersByCondition(dataset, objname, ClusterID)
Finds the top differentially expressed genes between conditions within a given cluster. 
For example, "What's differentially expressed between tumor vs. egressed populations in cluster 3?"
It will save text files with the top 50 up- and down- regulated genes between conditions. Additionally, ALL DE genes will be saved as a text file. Also, a .rnk file will be saved with only the top 50 up- and down- regulated genes and their L2FC values to use for GSEA. 
ClusterID is assumed to be the value in seurat_clusters, but is only used for file names. 
Also plots heatmaps with the top 50 up and down genes. 
This function is used by FindMarkersByCondtitionEachCluster(). 

## compare_clusters(dataset, objname, ClusterID1, ClusterID2)
Finds the top differentially expressed genes between two clusters. 
For example, "What's differentially expressed between clusters 9 and 1?"
It will save text files with the top 50 up- and down- regulated genes between those clusters. Additionally, All DE genes will be saved as a text file and a .rnk file for GSEA. Also plots heatmaps with the top 50 up and down genes. 
Example call: 
compare_clusters(combined.DCs, get_obj_name(combined.DCs), "9", "1")

## plot_huang2021(dataset, objname)
Makes a bunch of plots with all of the markers and subpopulations defined in Huang et al 2021, Figure 6 (innervated lymph node paper). 
Saves the plots in a subdirectory of your current working directory. 
https://www.sciencedirect.com/science/article/pii/S0092867420315646?via%3Dihub
Makes violin plots, feature plots, and violin plots split by condition. (same as plot_markers() above)

## plot_FilioDCMarkers(datset, objname)
Makes a bunch of plots with various DC markers that Filio sent me in 2021. 
Saves the plots in a subdirectory of your current working directory. 
Makes violin plots, feature plots, and violin plots split by condition. (same as plot_markers() above)

## plot_TaylorCD8Markers(datset, objname)
Makes a bunch of plots with various markers of interest for CD8 cells that Taylor sent me in 2021. 
Saves the plots in a subdirectory of your current working directory. 
Makes violin plots, feature plots, and violin plots split by condition. (same as plot_markers() above)

## doPseudotime(datset, objname)
Performs pseudotime analysis using Monocle3. Requires user input to choose the starting point for pseudotime. Outputs plots showing trajectories on the UMAP. 

## doGSEA(datset, objname)
Finds differentially expressed genes between conditions and performs GSEA on a seurat object. This function is used by FindMarkersByConditionEachCluster() 

