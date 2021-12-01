print('==========================================================')
print("Loading libraries: optparse, Matrix, dplyr, pagoda2, conos")
library("optparse")
library("Matrix")
library("dplyr")
library("pagoda2")
library("conos")
print("Loaded libraries: optparse, Matrix, dplyr, pagoda2, conos")
print('==========================================================')
options(ggrepel.max.overlaps = Inf)

# # ====================================

# # PARAMETERS for Leiden communities and walktrap
# runleiden <- TRUE
# leiden_resolution = 1
# runwalktrap <- TRUE
# walktrap_steps = 10
# # ====================================
# # PARAMETERS for UMAP
# umap_distance = 0.05
# umap_spread = 5
# # ====================================

# Parse input arguments
parser = OptionParser()
# ====================================
parser <- add_option(parser, c("--conos_object"), help = "RDS file created by Conos.Preprocess.")
# ====================================
# PARAMETERS for Leiden
parser <- add_option(parser, c("--runleiden"),type='character',default='True', help = "Whether or not to perform Leithen community clustering")
parser <- add_option(parser, c("--leiden_resolution"),type='double', default=1.0, help = "Resolution for Leiden Community Clustering (default:1). Generally best performance is achieved  with a resolution between 1 and 2")
# ====================================
# PARAMETERS for walktrap
parser <- add_option(parser, c("--runwalktrap"),type='character',default='True', help = "Whether or not to perform runwalktrap community clustering")
parser <- add_option(parser, c("--walktrap_steps"),type='integer', default=10, help = "Stepts to be taken for the walktrap algorithm (default:10).")
# ====================================
# PARAMETERS for UMPAP
parser <- add_option(parser, c("--umap_distance"),type='double',default=0.05, help = "UMAP Minimum Distance (default: 0.05)")
parser <- add_option(parser, c("--umap_spread"),type='double', default=5.0, help = "UMAP Spread (default: 5.0)")
# ====================================
print('==========================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('==========================================================')
# ====================================

# Restore the object
conos_object <- readRDS(file = args$conos_object)  # reads in a varaible called 'conos_object'
con <- conos_object$con
con_space <- conos_object$con_space
data_source <- conos_object$data_source

if (args$runleiden == "True") {
 runleiden = TRUE
 resol <- args$leiden_resolution
} else {
 runleiden = FALSE
 resol <- FALSE
}

if (args$runwalktrap == "True") {
 runwalktrap = TRUE
 stepnum <- args$walktrap_steps
} else {
 runwalktrap = FALSE
 stepnum <- FALSE
}

udist <- args$umap_distance
uspread <- args$umap_spread

# ===

print(paste("Begin UMAP Embedding:", Sys.time()))
con$embedGraph(method = "UMAP", min.dist = udist, spread = uspread)
print(paste("Finished UMAP Embedding:", Sys.time()))


if (runleiden) {
 print(paste("Finding Leiden Communities:", Sys.time()))
 con$findCommunities(method = leiden.community, resolution = resol, min.group.size=2)

 if (data_source == "matrix") {
  persample_emb = NULL
  emb_space = "tSNE"
 }

 if (data_source == "seurat") {
  # Check for which dimensionality reduction embeddings exist in the objects
  tsne_embeddings <- rep(NA, length(con$samples))
  for (i in seq_len(length(con$samples))) {
   tsne_embeddings[i] <- !is.null(con$samples[[i]]@reductions$tsne)
  }
  tsne_embeddings <- all(tsne_embeddings == TRUE)

  umap_embeddings <- rep(NA, length(con$samples))
  for (i in seq_len(length(con$samples))) {
   umap_embeddings[i] <- !is.null(con$samples[[i]]@reductions$umap)
  }
  umap_embeddings <- all(umap_embeddings == TRUE)

  if (tsne_embeddings == TRUE) {
   persample_emb = "tsne"
   emb_space = "tSNE"
  }
  if (umap_embeddings == TRUE) {
   persample_emb = "umap"
   emb_space = "UMAP"
  }
 }

#  ## Capture (C)PCA space embedded global leiden communities
#  png(paste0("DefaultVIS_", con_space, "_Leiden", resol, "_Clusters.png"), width = 16, 
#   height = 9, units = "in", res = 300)
#  print(cowplot::plot_grid(con$plotGraph(alpha = 0.1, clustering = "leiden"), con$plotGraph(alpha = 0.1, 
#   color.by = "sample", mark.groups = F, show.legend = T, legend.position = "bottom", 
#   legend.title = "")))
#  dev.off()

 ## Capture Leiden community composition
 png(paste0("Leiden", resol, "_Cluster_Composition.png"), width = 16, height = 9, 
  units = "in", res = 300)
 print(plotClusterBarplots(con, clustering = "leiden", legend.height = 0.2))
 dev.off()

 if (data_source == "matrix") {
  leiden.de <- con$getDifferentialGenes(clustering = "leiden", append.auc = TRUE, 
   groups = con$clusters$leiden$groups)
  capture.output(leiden.de, file = paste0("Leiden", resol, "_Cluster_Differential_Genes.txt"))
  png(paste0("Leiden", resol, "_Cluster_Top5DE_Heatmap.png"), width = 16, height = 9, 
   units = "in", res = 300)
  print(plotDEheatmap(con, as.factor(con$clusters$leiden$result$membership), 
   leiden.de, n.genes.per.cluster = 5, column.metadata = list(samples = con$getDatasetPerCell()), 
   row.label.font.size = 7))
  dev.off()
 } else {
  print(paste0("Harmonized cluster marker gene detection is not currently supported for Seurat objects. See: https://github.com/kharchenkolab/conos/issues/16 for details."))
 }
}
# ===

if (runwalktrap) {
 ## Generate Walktrap clusters
 print(paste("Begin Walktrap:", Sys.time()))
 con$findCommunities(method = igraph::walktrap.community, steps = stepnum, min.group.size=2)
 print(paste("Finished Walktrap:", Sys.time()))

 print("About to save figures.")
 ## Capture Walktrap community composition
 png(paste0("Walktrap", stepnum, "_Cluster_Composition.png"), width = 16, height = 9, 
  units = "in", res = 300)
 print(plotClusterBarplots(con, clustering = "walktrap", legend.height = 0.2))
 dev.off()
 print("Done saving figures.")

 if (data_source == "matrix") {
  walktrap.de <- con$getDifferentialGenes(clustering = "walktrap", append.auc = TRUE, 
   groups = con$clusters$walktrap$groups)
  capture.output(walktrap.de, file = paste0("Walktrap", stepnum, "_Cluster_Differential_Genes.txt"))
  png(paste0("Walktrap", stepnum, "_Cluster_Top5DE_Heatmap.png"), width = 16, 
   height = 9, units = "in", res = 300)
  print(plotDEheatmap(con, as.factor(con$clusters$walktrap$groups), leiden.de, 
   n.genes.per.cluster = 5, column.metadata = list(samples = con$getDatasetPerCell()), 
   row.label.font.size = 7))
  dev.off()
 } else {
  print(paste0("Harmonized cluster marker gene detection is not currently supported for Seurat objects. See: https://github.com/kharchenkolab/conos/issues/16 for details."))
 }
}

## =============================================================== Return
## runleiden, runwalktrap, and data_source

# Save an object to a file
saveRDS(list(con = con, runleiden = runleiden, runwalktrap = runwalktrap, data_source = data_source), 
 "conos_cluster_output.rds")
print("saved conos_cluster_output.rds")

# Save an object to a file saveRDS(con, 'conos_object.rds') print('saved
# conos_object.rds')

## ===============================================================
"Plot gene expression on the joint graph"

print("About to save figures.")
if (runleiden == TRUE) {
 ## Capture per-sample global leiden communities in UMAP space
 png(paste0("Per-sample_Global_Leiden", resol, "_Clusters_Individual_", emb_space, 
  ".png"), width = 16, height = 9, units = "in", res = 300)
 print(con$plotPanel(font.size = 4, clustering = "leiden", embedding = persample_emb))
 dev.off()

 png(paste0("Per-sample_Global_Leiden", resol, "_Clusters_Common_UMAP.png"), width = 16, 
  height = 9, units = "in", res = 300)
 print(con$plotPanel(font.size = 4, clustering = "leiden", use.common.embedding = TRUE))
 dev.off()

 ## Capture CPCA space embedded global leiden communities UMAP visualization
 png(paste0("Leiden", resol, "_Clusters_UMAP.png"), width = 16, height = 9, units = "in", 
  res = 300)
 print(cowplot::plot_grid(con$plotGraph(alpha = 0.1, clustering = "leiden", embedding = "UMAP"), 
  con$plotGraph(alpha = 0.1, color.by = "sample", embedding = "UMAP", mark.groups = F, 
   show.legend = T, legend.position = "bottom", legend.title = "")))
 dev.off()
}

if (runwalktrap == TRUE) {
 ## Capture per-sample global walktrap communities in UMAP space
 png(paste0("Per-sample_Global_Walktrap", stepnum, "_Clusters_Individual_", emb_space, 
  ".png"), width = 16, height = 9, units = "in", res = 300)
 print(con$plotPanel(font.size = 4, clustering = "walktrap", embedding = persample_emb))
 dev.off()

 png(paste0("Per-sample_Global_Walktrap", stepnum, "_Clusters_Common_UMAP.png"), 
  width = 16, height = 9, units = "in", res = 300)
 print(con$plotPanel(font.size = 4, clustering = "walktrap", use.common.embedding = TRUE))
 dev.off()

 ## Capture CPCA space embedded global Walktrap communities UMAP visualization
 png(paste0("Walktrap", stepnum, "_Clusters_UMAP.png"), width = 16, height = 9, 
  units = "in", res = 300)
 print(cowplot::plot_grid(con$plotGraph(alpha = 0.1, clustering = "walktrap", 
  embedding = "UMAP"), con$plotGraph(alpha = 0.1, color.by = "sample", embedding = "UMAP", 
  mark.groups = F, show.legend = T, legend.position = "bottom", legend.title = "")))
 dev.off()
}

if (runleiden == TRUE & runwalktrap == TRUE) {
 ## Capture comparison of Walktrap and Leiden clusters in UMAP Space
 png(paste0("Leiden", resol, "_vs_Walktrap", stepnum, "_Clusters_UMAP.png"), width = 16, 
  height = 9, units = "in", res = 300)
 print(cowplot::plot_grid(con$plotGraph(alpha = 0.1, clustering = "leiden", embedding = "UMAP"), 
  con$plotGraph(alpha = 0.1, clustering = "walktrap", embedding = "UMAP"), 
  con$plotGraph(alpha = 0.1, color.by = "sample", embedding = "UMAP", mark.groups = F, 
   show.legend = T, legend.position = "bottom", legend.title = "")))
 dev.off()
}

print("Done saving figures.")

print("Done!")
