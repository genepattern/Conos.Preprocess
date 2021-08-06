print('==========================================================')
print("Loading libraries: optparse, Matrix, dplyr, pagoda2, conos")
library("optparse")
library("Matrix")
library("dplyr")
library("pagoda2")
library("conos")
library("tools")
print("Loaded libraries: optparse, Matrix, dplyr, pagoda2, conos")
print('==========================================================')

# # ====================================
# # PARAMETERS for pagoda2
# setk = 40
# setperplexity = 30
# setodgenes = 6000
# # ====================================
# # PARAMETERS for PCA/cPCA
# con_space='PCA' # option PCA or CPCA
# con_comps=100
# con_odgenes=6000
# # ====================================


# Parse input arguments
parser = OptionParser()
# ====================================
parser <- add_option(parser, c("--file_list"), help = "List of files to load.")
# ====================================
# PARAMETERS for pagoda2
parser <- add_option(parser, c("--k"), type='integer', default=40, help = "Pagoda2: default number of neighbors to use in kNN graph")
parser <- add_option(parser, c("--perplexity"),type='integer', default=50, help = "Pagoda2: perplexity to use in generating tSNE and largeVis embeddings (default=50)")
parser <- add_option(parser, c("--pagoda_odgenes"),type='integer', default=3000, help = "Pagoda2: number of top overdispersed genes to use (default=3e3)")
# ====================================
# PARAMETERS for PCA/cPCA
parser <- add_option(parser, c("--projection_method"),type='character',default='PCA', help = "Whether to perform PCA or CPCA")
parser <- add_option(parser, c("--ncomps"),type='integer', default=50, help = "How many PCA components to use (?) (default=50)")
parser <- add_option(parser, c("--conos_odgenes"),type='integer', default=1000, help = "Pagoda2: number of top overdispersed genes to use (default=1e3)")
# ====================================
print('==========================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('==========================================================')
# ====================================
# PARAMETERS for pagoda2
setk = args$k
setperplexity = args$perplexity
setodgenes = args$pagoda_odgenes
# ====================================
# PARAMETERS for PCA/cPCA
if((args$projection_method == 'PCA') || (args$projection_method == 'pca') ){
  con_space= 'PCA'
}else{
        if((args$projection_method == 'CPCA') || (args$projection_method == 'cpca') ){
        con_space= 'CPCA'
        }else{
        print('projection_method not recognized, only PCA or CPCA are acceptable; value provided was')
        print(args$projection_method)
        quit()
        }
      }
con_comps=args$ncomps
con_odgenes=args$conos_odgenes
# ====================================


print(args$file_list)
con <- file(args$file_list, open = "r")
lines = readLines(con)

panel = list(NA)
i = 1
for (line in lines) {
 print("About to read")
 print(line)
 if (file_ext(line) == "rds") {
  RDSinput = readRDS(line)
  name = tail(strsplit(line, "/")[[1]], 1)
  name <- gsub("\\.rds$", "", name)
  print(paste0("Using ", name, " as the name of the column"))
  panel[[i]] = RDSinput
  names(panel)[[i]] = as.character(name)
  i = i + 1
 } else {
  name = tail(strsplit(line, "/")[[1]], 1)
  print(paste0("Using ", name, " as the name of the column"))
  readTable = read.table(line, sep = "\t", header = TRUE, )

  # readTable = readTable[1:500,1:501] write.table(readTable, file =
  # paste('small_500x500',name,sep='_'), sep = '\t',row.names = F, quote=F)

  row.names(readTable) <- readTable$symbol
  readTable[1] <- NULL

  # print(readTable[1:6,1:6])

  panel[[i]] = as(as.matrix(readTable), "dgCMatrix")
  names(panel)[[i]] = as.character(name)
  i = i + 1
 }
}
print("finished with the panel")
print(str(panel, 1))

classes <- rep(NA, length(panel))
for (i in seq_len(length(panel))) {
 classes[i] <- class(panel[[i]])
}

if (all(classes == "dgCMatrix")) {
 print(head(colnames(panel[[1]])))
 print(any(duplicated(unlist(lapply(panel, colnames)))))
 print("About to save figures.")
 panel.preprocessed <- lapply(panel, basicP2proc, min.cells.per.gene = 1, k = setk, 
  perplexity = setperplexity, n.odgenes = setodgenes, get.largevis = FALSE, 
  make.geneknn = FALSE)
 mode <- "matrix"
} else if (all(classes == "Seurat")) {
 panel.preprocessed <- panel
 mode <- "seurat"
 # Check for which dimensionality reduction embeddings exist in the objects
 tsne_embeddings <- rep(NA, length(panel.preprocessed))
 for (i in seq_len(length(panel))) {
  tsne_embeddings[i] <- !is.null(panel.preprocessed[[i]]@reductions$tsne)
 }
 tsne_embeddings <- all(tsne_embeddings == TRUE)

 umap_embeddings <- rep(NA, length(panel.preprocessed))
 for (i in seq_len(length(panel))) {
  umap_embeddings[i] <- !is.null(panel.preprocessed[[i]]@reductions$umap)
 }
 umap_embeddings <- all(umap_embeddings == TRUE)

 if (tsne_embeddings == FALSE && umap_embeddings == FALSE) {
  stop(paste0("No valid common dimensionality reduction was detected in the supplied seurat objects. Please run Seurat.Clustering with the same dimensionality reduction parameters."))
 }

} else {
 stop(paste0("Mixed type datasets are not currently supported. Data types reported as: ", 
  paste0(unique(classes), collapse = ", ")))
}

con <- Conos$new(panel.preprocessed, n.cores = 1)  # n.cores=1 is just so TNSE is reproducible. This is okay for smaller datastets, I have not tested it on larger ones.
panel <- NULL
gc()
# print(con)
if (mode == "matrix") {
 png("Sample_Independent_Clusters_tSNE.png", width = 16, height = 9, units = "in", 
  res = 300)
 # plot(x, y) # Make plot
 print(con$plotPanel(clustering = "multilevel", use.local.clusters = T, title.size = 4))
 dev.off()
}
if (mode == "seurat") {
 if (umap_embeddings == TRUE) {
  png("Sample_Independent_Clusters_UMAP.png", width = 16, height = 9, units = "in", 
   res = 300)
  print(con$plotPanel(embedding = "umap", use.local.clusters = T, title.size = 4))
  dev.off()
 }
 if (tsne_embeddings == TRUE) {
  png("Sample_Independent_Clusters_tSNE.png", width = 16, height = 9, units = "in", 
   res = 300)
  print(con$plotPanel(embedding = "tsne", use.local.clusters = T, title.size = 4))
  dev.off()
 }
}

print("Done saving figures.")

panel.preprocessed <- NULL
gc()

con$buildGraph(k = 15, k.self = 5, space = con_space, ncomps = con_comps, n.odgenes = con_odgenes, 
 matching.method = "mNN", metric = "angular", score.component.variance = TRUE, 
 verbose = TRUE)
print(paste("Finished Projections:", Sys.time()))

## Write a table of Cell to Sample Relationships useful for later
write.table(con$getDatasetPerCell(), file = "Cell_to_Sample_Memberships.txt", sep = "\t", 
 quote = FALSE, col.names = FALSE)

print("About to save figures.")
## Capture CPCA communities variance
png(paste0(con_space, "_Variance.png"), width = 16, height = 9, units = "in", res = 300)
print(plotComponentVariance(con, space = con_space))
print("Done saving figures.")


# Save an object to a file
saveRDS(list(con = con, con_space = con_space, mode = mode), "conos_preprocess_output.rds")
print("saved conos_preprocess_output.rds")

# # Restore the object readRDS(file = 'conos_preprocess_output.rds') # reads in a
# varaible called 'con'

print("Done!")
