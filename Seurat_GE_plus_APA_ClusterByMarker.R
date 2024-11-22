
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)

GE4 <- readRDS(paste0("./Seurat/GE.RData"))

##### Setup the Seurat Object

# Initialize the Seurat object with the raw (non-normalized data).
SeuratObject <- CreateSeuratObject(counts = GE4, project = "GSE116237", min.cells = 3, min.features = 200)
SeuratObject
SeuratObject@meta.data %>% head()

sra <- read.csv("./scRNASeqData_GSE116237.csv", header = T)

sra$Age <- gsub("before treatment","T00",sra$Age)
sra$Age <- gsub("4d on treatment","T04",sra$Age)
sra$Age <- gsub("28d on treatment","T28",sra$Age)
sra$Age <- gsub("57d on treatment","T57",sra$Age)

SeuratObject$timePoint<- sra$Age[match(SeuratObject@meta.data %>% row.names(),sra$Run)]

SeuratObject@meta.data %>% head()



# SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT")
# VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# SeuratObject <- subset(SeuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



##### Normalizing the data
SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
# SeuratObject <- NormalizeData(SeuratObject, normalization.method = "RC", scale.factor = 1e6)





##### Identification of highly variable features (feature selection)

SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top50 <- head(VariableFeatures(SeuratObject), 50)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SeuratObject)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
plot1
plot2


##### Scaling the data

all.genes <- rownames(SeuratObject)
SeuratObject <- ScaleData(SeuratObject, features = all.genes)


##### Perform linear dimensional reduction

SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(object = SeuratObject))

# Examine and visualize PCA results a few different ways
print(SeuratObject[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(SeuratObject, dims = 1:2, reduction = "pca")

DimPlot(SeuratObject, reduction = "pca")

DimHeatmap(SeuratObject, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(SeuratObject)

SeuratObject <- FindNeighbors(SeuratObject, dims = 1:15)
SeuratObject <- FindClusters(SeuratObject, resolution = 0.5)

SeuratObject <- RunUMAP(SeuratObject, dims = 1:15)

DimPlot(SeuratObject, reduction = "umap")
pA <- DimPlot(SeuratObject, reduction = "umap", cols = rep("lightgrey",6))

SeuratObject.markers <- FindAllMarkers(SeuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
x <- SeuratObject.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

p0 <- FeaturePlot(SeuratObject, features = x$gene)


id_convert <- read.table("id convert.tbl", header = F, sep = "\t")
names(id_convert) <- c("gene_id","gene_symbol")
id_convert2 <- id_convert[!duplicated(id_convert), ]
id_convert3 <- id_convert2 %>% separate(gene_id, sep = "\\.", into = c("gene_id",NA), remove = T)

signatures <- read.csv("/GSE116237/PIIS0092867418307931_signatures.csv", quote = "", na.strings = c("","NA"))




############################### add APA data


# library(APAlyzer)
# library(Rsamtools)
# library(easyRNASeq)
# 
# library(tidyverse)
# library(RColorBrewer)
# 
# library(dplyr)
# library(Seurat)
# library(patchwork)
# 
# 
# 
# 
# tbl_3UTR <- read.table("./APAlyzer/3UTR_hg19.tbl", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
# # tbl_IPA <- read.table("./APAlyzer/IPA_hg19.tbl", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
# 
# 
# RE2 <- tbl_3UTR[grepl("gene_symbol|_3UTR_RE", names(tbl_3UTR))]
# RE3 <- RE2 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
# # RE3 <- na.omit(RE3)
# 
# 
# RE4 <- data.frame(RE3[-c(1)],row.names=RE3$gene_symbol)
# 
# # saveRDS(RE4, file=paste0("./Seurat/RE.RData"))
# 
# # write.csv(RE4, paste0("./Seurat/RE.csv"), row.names = T)



library(tidyverse)
library(RColorBrewer)

library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)

RE4 <- readRDS(paste0("./Seurat/RE_3UTRAPA.RData"))
names(RE4) <- gsub("_3UTR_RE","",names(RE4))
RE5 <- RE4

RE5$median_RE <- apply(RE5,1,median,na.rm=T)
ncol <- ncol(RE5)
RE5[c(1:(ncol-1))] <- RE5[c(1:(ncol-1))]-RE5$median_RE

RE6 <- RE5[c(1:(ncol-1))]


# create a new assay to store APA information
APA_assay  <- CreateAssayObject(counts = RE6)
APA_assay


# add this assay to the previously created Seurat object
SeuratObject[["APA"]] <- APA_assay

# Validate that the object now contains multiple assays
Assays(SeuratObject)

# Extract a list of features measured in the APA assay
rownames(SeuratObject[["APA"]])

# Note that we can easily switch back and forth between the two assays to specify the default
# for visualization and analysis

# List the current default assay
DefaultAssay(SeuratObject)

# Switch the default to APA
DefaultAssay(SeuratObject) <- "APA"
DefaultAssay(SeuratObject)

# Now, we will visualize GE and APA By setting the default assay, we can
# visualize one or the other
DefaultAssay(SeuratObject) <- "APA"
p1 <- FeaturePlot(SeuratObject, features = x$gene, cols = c("blue", "red")) + ggtitle("APA")
DefaultAssay(SeuratObject) <- "RNA"
p2 <- FeaturePlot(SeuratObject, features = x$gene, cols = c("blue", "red")) + ggtitle("GE")

# # place plots side-by-side
# p0 | p1 | p2




RE7 <- RE6
RE7.median0 <- RE7.median <- data.frame(RE7.median=apply(RE7,2,median,na.rm=T))
names(RE7.median) <- "medianRED"
RE7.median <- t(RE7.median)

assay  <- CreateAssayObject(counts = RE7.median)
SeuratObject[["medianRED"]] <- assay



DefaultAssay(SeuratObject) <- "medianRED"
FeaturePlot_scCustom(seurat_object = SeuratObject, features = "medianRED", na_cutoff = NA) & scale_color_gradientn(colours = c("blue", "white","red"), limits = c(-0.5, 0.5), oob = scales::squish)

# saveRDS(SeuratObject, file = "./Seurat/SeuratObject_GE_APA.rds")

RE7.median1 <- RE7.median0
names(RE7.median1) <- "median_RED"
RE7.median1$Run <- row.names(RE7.median1)

sra_heatmap2 <- merge(RE7.median1,sra[c(1,2)])
sra_heatmap2$medianRED.T00 <- sra_heatmap2$median_RED
sra_heatmap2$medianRED.T04 <- sra_heatmap2$median_RED
sra_heatmap2$medianRED.T28 <- sra_heatmap2$median_RED
sra_heatmap2$medianRED.T57 <- sra_heatmap2$median_RED
sra_heatmap2[sra_heatmap2$Age!="T00",]$medianRED.T00 <- NA
sra_heatmap2[sra_heatmap2$Age!="T04",]$medianRED.T04 <- NA
sra_heatmap2[sra_heatmap2$Age!="T28",]$medianRED.T28 <- NA
sra_heatmap2[sra_heatmap2$Age!="T57",]$medianRED.T57 <- NA

REDbyT <- t(data.frame(sra_heatmap2[-c(1:3)],row.names=sra_heatmap2$Run))
assay  <- CreateAssayObject(counts = REDbyT)
SeuratObject[["medianREDbyT"]] <- assay
rownames(SeuratObject[["medianREDbyT"]])



DefaultAssay(SeuratObject) <- "medianREDbyT"
FeaturePlot_scCustom(seurat_object = SeuratObject, features = "medianRED.T00", na_cutoff = NA, alpha_exp = 0.75) & scale_color_gradientn(na.value = "grey90", colours = c("blue", "white","red"), limits = c(-0.5, 0.5), oob = scales::squish)|
FeaturePlot_scCustom(seurat_object = SeuratObject, features = "medianRED.T04", na_cutoff = NA, alpha_exp = 0.75) & scale_color_gradientn(na.value = "grey90", colours = c("blue", "white","red"), limits = c(-0.5, 0.5), oob = scales::squish)|
FeaturePlot_scCustom(seurat_object = SeuratObject, features = "medianRED.T28", na_cutoff = NA, alpha_exp = 0.75) & scale_color_gradientn(na.value = "grey90", colours = c("blue", "white","red"), limits = c(-0.5, 0.5), oob = scales::squish)|
FeaturePlot_scCustom(seurat_object = SeuratObject, features = "medianRED.T57", na_cutoff = NA, alpha_exp = 0.75) & scale_color_gradientn(na.value = "grey90", colours = c("blue", "white","red"), limits = c(-0.5, 0.5), oob = scales::squish)


# SeuratObject <- readRDS(paste0("./Seurat/SeuratObject_GE_APA.rds"))


# SeuratObject0 <- SeuratObject
# SeuratObject <- SeuratObject0


DefaultAssay(SeuratObject) <- "RNA"
GE4.counts <- data.frame(GetAssayData(object = SeuratObject, slot = "counts")) ##get raw read counts
GE4.norm <- data.frame(GetAssayData(object = SeuratObject, slot = "data")) ##get normalized data
GE4.scaled <- data.frame(GetAssayData(object = SeuratObject, slot = "scale.data")) ##get normalized + scaled data



plotListRNA1 <- list()
plotListRNA2 <- list()
plotListRNA3 <- list()
plotListRNA4 <- list()

plotListAPA1 <- list()
plotListAPA2 <- list()
plotListAPA3 <- list()
plotListAPA4 <- list()

for (col in colnames(signatures)){
  
  # col <- "mitosis"
  
  genes <- signatures[[col]][!is.na(signatures[[col]])]
  genes2 <- id_convert3$gene_symbol[match(genes, id_convert3$gene_id)]
  
  # gene_list <- list(genes2)
  
  DefaultAssay(SeuratObject) <- "RNA"
  
  geneSet <- subset(GE4.scaled, row.names(GE4.scaled) %in% genes2)
  geneSet.median <- data.frame(geneSet.median=apply(geneSet,2,median,na.rm=T))
  names(geneSet.median) <- col
  geneSet.median2 <- geneSet.median %>% mutate_at(c(col), ~(scale(.) %>% as.vector))
  
  Assay.geneSet.median <- t(geneSet.median2)
  assay  <- CreateAssayObject(counts = Assay.geneSet.median)
  SeuratObject[[paste0("GE",col)]] <- assay
  
  geneSet.median$cell <- row.names(geneSet.median)
  geneSet.median <- geneSet.median[order(geneSet.median[[col]],decreasing=F),]
  bin1=bin2=round(nrow(geneSet.median)/3)
  bin3=nrow(geneSet.median)-(bin1+bin2)
  geneSet.median$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3))
  SeuratObject[[paste0(col,"_explvl")]] <- geneSet.median$bin[match(SeuratObject@meta.data %>% row.names(),geneSet.median$cell)]
  
  sra[[paste0(col,"_explvl")]] <- geneSet.median$bin[match(sra$Run,geneSet.median$cell)]
  sra[[paste0(col,"_explvl_timePoint")]] <- paste0(sra$Age,sra[[paste0(col,"_explvl")]])
  SeuratObject[[paste0(col,"_explvl_timePoint")]] <- sra[[paste0(col,"_explvl_timePoint")]][match(SeuratObject@meta.data %>% row.names(),sra$Run)]
  
  
  
  DefaultAssay(SeuratObject) <- paste0("GE",col)
  
  # SeuratObject <- AddModuleScore(object = SeuratObject, features = gene_list, name = col)
  plotListRNA1[[paste0("GE",col,1)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), na_cutoff = NA) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(Assay.geneSet.median, na.rm = T))), max(abs(range(Assay.geneSet.median, na.rm = T)))))
  plotListRNA2[[paste0("GE",col,2)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), split.by = "timePoint", na_cutoff = NA) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(Assay.geneSet.median, na.rm = T))), max(abs(range(Assay.geneSet.median, na.rm = T)))))
  plotListRNA3[[paste0("GE",col,3)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), split.by = paste0(col,"_explvl"), na_cutoff = NA) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(Assay.geneSet.median, na.rm = T))), max(abs(range(Assay.geneSet.median, na.rm = T)))))
  plotListRNA4[[paste0("GE",col,4)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), split.by = paste0(col,"_explvl_timePoint"), na_cutoff = NA, num_columns = 3) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(Assay.geneSet.median, na.rm = T))), max(abs(range(Assay.geneSet.median, na.rm = T)))))
  

  
  RE.signature <- subset(RE6, row.names(RE6) %in% genes2)
  RE.signature.median <- data.frame(RE.signature.median=apply(RE.signature,2,median,na.rm=T))
  names(RE.signature.median) <- col
  RE.signature.median <- RE.signature.median %>% mutate_at(c(col), ~(scale(.) %>% as.vector))
  
  RE.signature.median <- t(RE.signature.median)
  assay  <- CreateAssayObject(counts = RE.signature.median)
  SeuratObject[[paste0("RED",col)]] <- assay
  
  DefaultAssay(SeuratObject) <- paste0("RED",col)
  plotListAPA1[[paste0("RED",col,1)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), na_cutoff = NA) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(RE.signature.median, na.rm = T))), max(abs(range(RE.signature.median, na.rm = T)))))
  plotListAPA2[[paste0("RED",col,2)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), split.by = "timePoint", na_cutoff = NA) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(RE.signature.median, na.rm = T))), max(abs(range(RE.signature.median, na.rm = T)))))
  plotListAPA3[[paste0("RED",col,3)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), split.by = paste0(col,"_explvl"), na_cutoff = NA) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(RE.signature.median, na.rm = T))), max(abs(range(RE.signature.median, na.rm = T)))))
  plotListAPA4[[paste0("RED",col,4)]] <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), split.by = paste0(col,"_explvl_timePoint"), na_cutoff = NA, num_columns = 3) &
    scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(RE.signature.median, na.rm = T))), max(abs(range(RE.signature.median, na.rm = T)))))
  
  

}

# write.csv(sra,"./scRNASeqData_GSE116237_by_markerGeneExpLvl_scaled.csv", row.names = F)

library(cowplot)
plot_grid(plotlist=c(plotListRNA1),ncol = 4)
plot_grid(plotlist=c(plotListRNA2),ncol = 2,labels = names(plotListRNA2))
plot_grid(plotlist=c(plotListRNA3),ncol = 2,labels = names(plotListRNA3))
plot_grid(plotlist=c(plotListRNA4[c(1:2)]),ncol = 2,labels = names(plotListRNA4[c(1:2)]))
plot_grid(plotlist=c(plotListRNA4[c(3:4)]),ncol = 2,labels = names(plotListRNA4[c(3:4)]))
plot_grid(plotlist=c(plotListRNA4[c(5:6)]),ncol = 2,labels = names(plotListRNA4[c(5:6)]))
plot_grid(plotlist=c(plotListRNA4[c(7:8)]),ncol = 2,labels = names(plotListRNA4[c(7:8)]))

plot_grid(plotlist=c(plotListAPA1),ncol = 4)
plot_grid(plotlist=c(plotListAPA2),ncol = 2,labels = names(plotListAPA2))
plot_grid(plotlist=c(plotListAPA3),ncol = 2,labels = names(plotListAPA3))
plot_grid(plotlist=c(plotListAPA4[c(1:2)]),ncol = 2,labels = names(plotListAPA4[c(1:2)]))
plot_grid(plotlist=c(plotListAPA4[c(3:4)]),ncol = 2,labels = names(plotListAPA4[c(3:4)]))
plot_grid(plotlist=c(plotListAPA4[c(5:6)]),ncol = 2,labels = names(plotListAPA4[c(5:6)]))
plot_grid(plotlist=c(plotListAPA4[c(7:8)]),ncol = 2,labels = names(plotListAPA4[c(7:8)]))


SeuratObject@meta.data %>% head()


p <- FeaturePlot_scCustom(seurat_object = SeuratObject, features = paste0(col), split.by = "timePoint", na_cutoff = NA) &
  scale_color_gradientn( colours = c("blue", "yellow","red"), limits = c(-max(abs(range(Assay.geneSet.median, na.rm = T))), max(abs(range(Assay.geneSet.median, na.rm = T)))))
