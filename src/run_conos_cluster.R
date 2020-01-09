print('==========================================================')
print("Loading libraries: optparse, Matrix, dplyr, pagoda2, conos")
library("optparse")
library("Matrix")
library("dplyr")
library("pagoda2")
library("conos")
print("Loaded libraries: optparse, Matrix, dplyr, pagoda2, conos")
print('==========================================================')

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
conos_object <- readRDS(file = args$conos_object) # reads in a varaible called 'conos_object'
con <- conos_object$con
con_space <- conos_object$con_space

if(args$runleiden=='True'){
  runleiden = TRUE
  resol <- args$leiden_resolution
}else{
  runleiden = FALSE
  resol <- FALSE
}

if(args$runwalktrap=='True'){
  runwalktrap = TRUE
  stepnum <- args$walktrap_steps
}else{
  runwalktrap = FALSE
  stepnum <- FALSE
}

umap_distance <- args$umap_distance
umap_spread <- args$umap_spread


#===

if(runleiden){
print(paste("Finding Leiden Communities:", Sys.time()))
con$findCommunities(method=leiden.community, resolution=resol)

## Capture per-sample global leiden communities
png(paste0("Per-sample_Global_Leiden",resol,"_Clusters.png"), width=16, height=9, units = 'in', res=300)
print(con$plotPanel(font.size=4, clustering='leiden'))
dev.off()

## Capture CPCA space embedded global leiden communities
png(paste0("DefaultVIS_",con_space,"_Leiden",resol,"_Clusters.png"), width=16, height=9, units = 'in', res=300)
print(cowplot::plot_grid(
con$plotGraph(alpha=0.1, clustering='leiden'),
con$plotGraph(alpha=0.1, color.by='sample', mark.groups=F, show.legend=T, legend.position='bottom', legend.title = "")))
dev.off()

## Capture Leiden community composition
png(paste0("Leiden",resol,"_Cluster_Composition.png"), width=16, height=9, units = 'in', res=300)
print(plotClusterBarplots(con, clustering='leiden', legend.height = 0.2))
dev.off()

leiden.de <- con$getDifferentialGenes(clustering='leiden')
capture.output(leiden.de, file = paste0("harmonized_cluster_markers_leiden",resol,".txt"))
}

#===

if(runwalktrap){
## Generate Walktrap clusters
print(paste("Begin Walktrap:", Sys.time()))
con$findCommunities(method = igraph::walktrap.community, steps=stepnum)
print(paste("Finished Walktrap:", Sys.time()))

print("About to save figures.")
## Capture Walktrap community composition
png("Walktrap_Cluster_Composition.png", width=16, height=9, units = 'in', res=300)
print(plotClusterBarplots(con, clustering='walktrap', legend.height = 0.2))
dev.off()
print("Done saving figures.")

walktrap.de <- con$getDifferentialGenes(clustering='walktrap')
capture.output(walktrap.de, file = "harmonized_cluster_markers_walktrap.txt")
}

##===============================================================
##===============================================================
##===============================================================
"Plot gene expression on the joint graph"

udist <- umap_distance
uspread <- umap_spread

print(paste("Begin UMAP Embedding:", Sys.time()))
con$embedGraph(method="UMAP", min.dist=udist, spread=uspread, n.cores=3) # n.cores=1 led to some incredibly slow times on sample data locally.
print(paste("Finished UMAP Embedding:", Sys.time()))

print("About to save figures.")
if(runleiden == TRUE){
## Capture per-sample global leiden communities in UMAP space
png("UMAP_per-sample_Global_Leiden_Clusters.png", width=16, height=9, units = 'in', res=300)
print(con$plotPanel(font.size=4, clustering='leiden'))
dev.off()

## Capture CPCA space embedded global leiden communities UMAP visualization
png("UMAP_space_Leiden_Clusters.png", width=16, height=9, units = 'in', res=300)
print(cowplot::plot_grid(
con$plotGraph(alpha=0.1, clustering='leiden'),
con$plotGraph(alpha=0.1, color.by='sample', mark.groups=F, show.legend=T, legend.position='bottom', legend.title = "")))
dev.off()
}

if(runwalktrap == TRUE){
## Capture per-sample global walktrap communities in UMAP space
png("UMAP_per-sample_Global_Leiden_Clusters.png", width=16, height=9, units = 'in', res=300)
print(con$plotPanel(font.size=4, clustering='walktrap'))
dev.off()

## Capture CPCA space embedded global Walktrap communities UMAP visualization
png("UMAP_space_Walktrap_Clusters.png", width=16, height=9, units = 'in', res=300)
print(cowplot::plot_grid(
con$plotGraph(alpha=0.1, clustering='walktrap'),
con$plotGraph(alpha=0.1, color.by='sample', mark.groups=F, show.legend=T, legend.position='bottom', legend.title = "")))
dev.off()
}

if(runleiden == TRUE & runwalktrap == TRUE){
## Capture comparison of Walktrap and Leiden clusters in UMAP Space
png("Leiden_vs_Walktrap_Clusters.png", width=16, height=9, units = 'in', res=300)
print(cowplot::plot_grid(
con$plotGraph(alpha=0.1, clustering='leiden'),
con$plotGraph(alpha=0.1, clustering='walktrap'),
con$plotGraph(alpha=0.1, color.by='sample', mark.groups=F, show.legend=T, legend.position='bottom', legend.title = "")))
dev.off()
}

print("Done saving figures.")

# Return runleiden, runwalktrap

# Save an object to a file
saveRDS(list(con=con,runleiden=runleiden,runwalktrap=runwalktrap), "conos_cluster_output.rds")
print('saved conos_cluster_output.rds')

# Save an object to a file
# saveRDS(con, "conos_object.rds")
# print('saved conos_object.rds')

print('Done!')
