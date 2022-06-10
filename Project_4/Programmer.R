#installing and calling packages 
BiocManager::install("tximport")
BiocManager::install("SeqGSEA")
BiocManager::install('biomaRt')
BiocManager::install('fishpond')
BiocManager::install('EnsDb.Hsapiens.v79')
BiocManager::install("GenomicFeatures")
BiocManager::install("grimbough/biomaRt")
install.packages("Seurat")
install.packages("Matrix")
install.packages("tidyverse")
library('biomaRt')
library("GenomicFeatures")
library('tidyverse')
library('Matrix')
library('Seurat')
library('tximport')
library('SeqGSEA')
library('fishpond')
library('EnsDb.Hsapiens.v79')

x <- file.path("/projectnb/bf528/users/glass_bottom_boat/marina/data_curator/alevin_output/alevin/quants_mat.gz")  #defining file path
txi<- tximport(x, type="alevin")  #reading in alevin counts using txi
dim(txi)
panc_cells <- CreateSeuratObject(counts = txi$counts, project = "proj4_panc", min.cells = 3, min.features = 200)
dim(panc_cells)
panc_cells

#mapping gene symbols
gene_names <- sub("[.][0-9]*$", "", panc_cells@assays$RNA@counts@Dimnames[[1]])
ensembl <- convertEnsembl2Symbol(gene_names)
gene_original_names <- ensembl$hgnc_symbol

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= gene_names, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

gene_original_names <- geneIDs$SYMBOL

panc_cells@assays$RNA@counts@Dimnames[[1]] <- gene_original_names
panc_cells@assays$RNA@data@Dimnames[[1]] <- gene_original_names

head(panc_cells)
dim(panc_cells)

panc_cells[["percent.mt"]] <- PercentageFeatureSet(panc_cells, pattern = "^MT-") #calculates percentages of counts and does quality check
head(panc_cells@meta.data, 5)
dim(panc_cells)

VlnPlot(panc_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #visualises as violin plot

plot1 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filtering low quality cells 
panc_cells <- subset(panc_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
panc_cells

#normalize the data
panc_cells <- NormalizeData(panc_cells)
panc_cells[["RNA"]]@meta.features <- data.frame(row.names = rownames(panc_cells[["RNA"]]))

#identify highly variable features 
panc_cells <- FindVariableFeatures(panc_cells, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(panc_cells), 10)
#plot variable features with and without labels
plot1 <- VariableFeaturePlot(panc_cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 
plot2
dim(panc_cells)
#scaling the data
all.genes <- rownames(panc_cells)
panc_cells <- ScaleData(panc_cells, features = all.genes)

#perform linear dimensional reduction
panc_cells <- RunPCA(panc_cells, features = VariableFeatures(object = panc_cells))
dim(panc_cells)
print(panc_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(panc_cells, dims = 1:2, reduction = "pca")
DimPlot(panc_cells, reduction = "pca")

#testing for association between observed and latent variables
panc_cells <- JackStraw(panc_cells, num.replicate = 100)
panc_cells <- ScoreJackStraw(panc_cells, dims = 1:20)
JackStrawPlot(panc_cells, dims = 1:15)
ElbowPlot(panc_cells)

#clustering the cells
panc_cells <- FindNeighbors(panc_cells, dims = 1:10)
panc_cells <- FindClusters(panc_cells, resolution = 0.5)
cluster_assignments <- Idents(panc_cells)
plot(cluster_assignments, labels(TRUE))

#run non-linear reduction
panc_cells <- RunUMAP(panc_cells, dims = 1:10)
DimPlot(panc_cells, reduction = "umap")
dim(panc_cells)

save(panc_cells, file = "newpanc_cells.rda")
save(panc_cells, file = "newpanc_cells.Rdata")
saveRDS(panc_cells, file = "newpanc_cells.rds")

