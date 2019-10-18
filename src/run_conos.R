print("Loading libraries: optparse,Matrix")
library("optparse")
library("Matrix")
library("dplyr")
library("pagoda2")
library("conos")

# Parse input arguments
parser = OptionParser()

parser <- add_option(parser, c("--file_list"), help = "List of files to load.")
args <- parse_args(parser)

print(args$file_list)
con <- file(args$file_list, open = "r")
lines = readLines(con)

panel = list(NA)
i = 1
for (line in lines){
  print("About to read")
  print(line)
  name = tail(strsplit(line,'/')[[1]],1)
  print(c('Using',name, 'as the name of the column'))
  readTable = read.table(line, sep='\t',header=TRUE,)

  # readTable = readTable[1:500,1:501]
  # write.table(readTable, file = paste("small_500x500",name,sep='_'), sep = "\t",row.names = F, quote=F)

  row.names(readTable) <- readTable$symbol
  readTable[1] <- NULL

  #print(readTable[1:6,1:6])

  panel[[i]] = as(as.matrix(readTable),"dgCMatrix")
  names(panel)[[i]] = as.character(name)
  i = i+1
}
print('finished with the panel')
print(str(panel,1))
print(head(colnames(panel[[1]])))
print(any(duplicated(unlist(lapply(panel,colnames)))))

##PARAMETERS
setk = 40
setperplexity = 30
setodgenes = 6000

panel.preprocessed <- lapply(panel, basicP2proc, min.cells.per.gene=1, k=setk, perplexity=setperplexity , n.odgenes=setodgenes, get.largevis=FALSE, make.geneknn=FALSE)
con <- Conos$new(panel.preprocessed, n.cores=1)
panel <- NULL
gc()
png("Sample_Independent_Clusters.png", width=16, height=9, units = 'in', res=300)
# plot(x, y) # Make plot

# print(con)
print(con$plotPanel(clustering="multilevel", use.local.clusters=T, title.size=4))
dev.off()

# PARAMETERS
con_space='PCA' # option PCA or CPCA
con_comps=100
con_odgenes=6000

panel.preprocessed <- NULL
gc()

con$buildGraph(k=15, k.self=5, space=con_space, ncomps=con_comps, n.odgenes=con_odgenes, matching.method='mNN', metric='angular', score.component.variance = TRUE, verbose = TRUE)
print(paste("Finished Projections:", Sys.time()))

## Write a table of Cell to Sample Relationships useful for later
write.table(con$getDatasetPerCell(), file = "Cell_to_Sample_Memberships.txt", sep="\t", quote= FALSE, col.names = FALSE)

## Capture CPCA communities variance
png("CPCA_Variance.png", width=16, height=9, units = 'in', res=300)
print(plotComponentVariance(con, space=con_space))

print('Done!')
