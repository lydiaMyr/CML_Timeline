.libPaths("D:/Program Files/R/R-3.5.3/library")
.libPaths()
library(Seurat)
library(SeuratData)
library(cowplot)
library(scater)
library(Matrix)
library(sva)
library(scran)
library(ggplot2)
library(dplyr)


# PART one: Load data and filter cells
# step1: Load data
setwd("E:/YXP_result")
YXP1.data <- Read10X("./filtered_feature_bc_matrix_YXP1")
YXP2.data <- Read10X("./filtered_feature_bc_matrix_YXP2")

dim(YXP1.data) # 31053  4400
dim(YXP2.data) # 31053  2780

#setwd("E:/YXP_result/result")
# write.table(as.matrix(raw.data),"raw_data.txt",sep="\t",col.names=T,row.names=T, quote=F)

# step2: Seurat CreateSeuratObject and Preliminary filter
# First, create Seurat objects for each of the datasets, and then merge into one large seurat object.
sdata.YXP1 <- CreateSeuratObject(YXP1.data, project = "YXP1", min.cells = 3, min.features = 200)
sdata.YXP2 <- CreateSeuratObject(YXP2.data, project = "YXP2", min.cells = 3, min.features = 200)

# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the datasets.
alldata <- merge(sdata.YXP1,sdata.YXP2, add.cell.ids=c("YXP1", "YXP2"))
alldata

# check number of cells from each sample, is stored in the orig.ident slot of metadata and is autmatically set as active ident.
table(Idents(alldata))
# KL1  KL2 
# 4316 2642

# step3: Calculate mitochondrial proportion and ribosomal proportion
alldata <- PercentageFeatureSet(alldata, pattern = "^MT-", col.name = "percent.mito")
alldata <- PercentageFeatureSet(alldata, pattern = "^Rp[Sl]", col.name = "percent.ribo")
head(alldata@meta.data)

# step4: Secondary filter

VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA","percent.mito"), pt.size = 0.1) + NoLegend()

# step4.1: select cells with percent.mito < 10
dim(alldata) # 17925  6958
# dim(subset(alldata, percent.mito >= 25)) # 17925   266
# dim(subset(alldata, percent.mito < 25)) # 17925  6692
# dim(subset(alldata, percent.mito <= 10)) # 17925  6038
# dim(subset(alldata, percent.mito < 5)) # 17925  3544
data.filt <- subset(alldata, percent.mito < 25)

# check number of cells
ncol(data.filt) # 6038
table(Idents(data.filt))
# KL1  KL2 
# 3858 2180

# step5: Original dataset in Seurat class
save(data.filt,file="original_seurat_object_YXP1-YXP2.RData")

# write.table(as.matrix(raw.data),"raw_data.txt",sep="\t",col.names=T,row.names=T, quote=F)
write.table(as.matrix(GetAssayData(object = data.filt, slot = "counts")),"data.filt.txt",sep="\t", col.names=T, row.names=T, quote=F)


######################################################################################
######################################################################################
# PART two: cluster cells
# step1: Load data
load("original_seurat_object_YXP1-YXP2.RData")
data.filt
dim(data.filt) # 18792 13677
head(data.filt@meta.data)

# bak.filt <- data.filt
# data.filt <- bak.filt 
# step2: normal workflow

## normalization respectively
YXP1_data <- subset(data.filt, orig.ident=="YXP1")
YXP2_data <- subset(data.filt, orig.ident=="YXP2")

# FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# length(VariableFeatures(immune.combined)) # 2000
YXP1_data <- NormalizeData(object = YXP1_data)
YXP1_data <- FindVariableFeatures(object = YXP1_data, selection.method = "vst")

YXP2_data <- NormalizeData(object = YXP2_data)
YXP2_data <- FindVariableFeatures(object = YXP2_data, selection.method = "vst") 

## Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = c(YXP1_data, YXP2_data), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, features = rownames(immune.combined)) # vars.to.regress = c("percent.mito", "percent.ribo", "orig.ident")
# immune.combined <- RunPCA(immune.combined, features = VariableFeatures(immune.combined), verbose = FALSE)
immune.combined <- ScaleData(immune.combined, features = rownames(immune.combined)) # vars.to.regress = c("percent.mito", "percent.ribo", "orig.ident")
immune.combined <- RunPCA(immune.combined, features = rownames(immune.combined), verbose = FALSE)

# # Choose the PC numbers
# We will start by calculating the first metric:
seurat_control <- immune.combined # from 15 to 40
# Determine percent of variation associated with each PC
pct <- seurat_control[["pca"]]@stdev / sum(seurat_control[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2

# Printing out the most variable genes driving PCs
print(x = seurat_control[["pca"]], dims = 1:50, nfeatures = 5)

ElbowPlot(seurat_control, ndims = 50)
use.pcs = 1:41
#use.pcs = 1:20
immune.combined <- FindNeighbors(immune.combined, reduction="pca", dims = use.pcs)
immune.combined <- FindClusters(object = immune.combined, resolution = seq(0.1,1.5,0.05), verbose = FALSE)
sapply(grep("res",colnames(immune.combined@meta.data),value = TRUE),function(x) length(unique(immune.combined@meta.data[,x])))
# integrated_snn_res.0.1 integrated_snn_res.0.15  integrated_snn_res.0.2 integrated_snn_res.0.25  integrated_snn_res.0.3 integrated_snn_res.0.35 
# 10                      10                      11                      12                      13                      13 
# integrated_snn_res.0.4 integrated_snn_res.0.45  integrated_snn_res.0.5 integrated_snn_res.0.55  integrated_snn_res.0.6 integrated_snn_res.0.65 
# 13                      15                      15                      17                      17                      18 
# integrated_snn_res.0.7 integrated_snn_res.0.75  integrated_snn_res.0.8 integrated_snn_res.0.85  integrated_snn_res.0.9 integrated_snn_res.0.95 
# 18                      18                      18                      20                      20                      20 
# integrated_snn_res.1 integrated_snn_res.1.05  integrated_snn_res.1.1 integrated_snn_res.1.15  integrated_snn_res.1.2 integrated_snn_res.1.25 
# 20                      20                      23                      23                      23                      24 
# integrated_snn_res.1.3 integrated_snn_res.1.35  integrated_snn_res.1.4 integrated_snn_res.1.45  integrated_snn_res.1.5 
# 24                      24                      24                      25                      25 


# t-SNE and Clustering
immune.combined <- RunTSNE(object = immune.combined, reduction = "pca", dims = use.pcs, do.fast = TRUE)
immune.combined <- RunUMAP(object = immune.combined, reduction = "pca", dims = use.pcs)

# table(Idents(data.filt),data.filt$orig.ident)
Idents(immune.combined) <- "integrated_snn_res.0.6"
table(Idents(immune.combined))
table(Idents(immune.combined),immune.combined@meta.data$orig.ident)
orig=table(Idents(immune.combined),immune.combined@meta.data$orig.ident)
cell_orig=data.frame(Count=as.vector(unlist(c(orig[,1],orig[,2]))),cell=rep(row.names(orig),times=2),orig=rep(c("T3","T4"),each=18))
cell_orig$cell <- factor(cell_orig$cell,levels=names(sort(orig[,1],decreasing = T)), ordered=TRUE)

pdf("single_cell_orig_distribution.pdf",15,10)
ggplot(cell_orig, aes(cell, Count, fill=orig,group=orig)) +
  geom_bar(stat="identity",position = "dodge")+theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

dev.off()

ggplot(cell_orig, aes(cell, Count, fill=orig,group=orig)) +
  geom_bar(stat="identity",position = "fill")+theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


data.filt <- immune.combined
# Visualization
DimPlot(data.filt, reduction = "pca") # split.by = "seurat_clusters"
DimPlot(data.filt, reduction = "tsne",label=TRUE)
DimPlot(data.filt, reduction = "umap",label=TRUE)

# Save plot
# # compare the distribution of different samples among umap
p3 <- DimPlot(object = data.filt, reduction = "umap", split.by = "orig.ident",label=T, label.size=4, pt.size=1)
p1 <- DimPlot(object = data.filt, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = data.filt, reduction = "umap", label = TRUE, pt.size=1)
p <- plot_grid(p1, p2)
# 
save_plot("Umap_compare_samples.png", p, base_height = 6, base_aspect_ratio = 2.5, base_width = NULL, dpi=600)
save_plot("Umap_compare_samples.pdf", p, base_height = 6, base_aspect_ratio = 2.5, base_width = NULL)

save_plot("Umap_split_samples.png", p3, base_height = 6, base_aspect_ratio = 2.5, base_width = NULL, dpi=600)
save_plot("Umap_split_samples.pdf", p3, base_height = 6, base_aspect_ratio = 2.5, base_width = NULL)

save(data.filt, file="Normalizeation_respectively_cluster_YXP1-YXP2.RData")

######################################################################################
######################################################################################
# PART three: Identifying markers and visualize features
# step1: Load data
load("Normalizeation_respectively_cluster_YXP1-YXP2.RData")
data.filt
head(data.filt@meta.data)
DefaultAssay(data.filt) <- "RNA"

# FindAllMarkers between clusters
markers_all <- FindAllMarkers(object = data.filt, only.pos = FALSE,min.cells.group=10)
dim(markers_all) 
head(markers_all)
table(table(markers_all$gene))
table(markers_all$cluster)
write.table(markers_all,"markers_all_DEGs_among_clusters_YXP1-YXP2.csv", sep=",", col.names=NA)
save(markers_all,file="markers_of_YXP1-YXP2_among_clusters.RData")

immune.combined=data.filt
immune.combined <- RenameIdents(immune.combined, `0` = "C1QTNF4+pro_B", `1` = "IGLC6+pro_B", `2` = "NK", `3` = "MT1X+NRIP1+pro_B", `4` = "CCL17+JCHAIN+pro_B", `5` = "Neural progenitor cell", `6` = "CD8+ T", `7` = "Erythrocytes", `8` = "CD4 central memory", `9`="Granulocyte-monocyte progenitor",`10`="Myelocyte",`11` = "Neutrophil", 
                                `12` = "IGHM+pro_B", `13` = "Monocyte",`14`="Memory B cell",`15`="NKT",`16`="plasma",`17`="Megakaryocyte") 
pdf("result/single_cell/singleCellCluster.pdf",25,10)
DimPlot(immune.combined, label = TRUE,split.by = "orig.ident", label.size = 6 )+NoLegend()
dev.off()
length(names(table(markers_all$gene))) 
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]
dim(markers_all_single) # 2018    7
head(markers_all_single)
table(markers_all_single$cluster)
write.table(markers_all_single,"markers_single_YXP1-YXP2.csv", sep=",", col.names=NA)

FeaturePlot(data.filt, features =c("CD79A","CD79B") ,  max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
FeaturePlot(data.filt, features =c("CD3D","CD3E","CD3G") ,  max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
FeaturePlot(data.filt, features =c("CD4","CD8A","CD8B") , max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
#naive effector
FeaturePlot(data.filt, features =c("PTPRC","CCR7") ,  max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
#NK
FeaturePlot(data.filt, features =c("FCGR3A","FCGR3B") ,  max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
#neutrophil
#MPO and Ly6G
pdf("result/single_cell/T3_up_gene.pdf",15,10)
for(gene in c("BTK","CD72","BLNK","CD19","CD79A","CD79B","CD22","LILRB4")){
  p=FeaturePlot(data.filt, features =gene ,  max.cutoff = 3, cols = c("grey", "red")
              ,min.cutoff = "q4",label=T, label.size=4, pt.size=1,split.by = "orig.ident")
  print(p)
}
dev.off()
#monocyte
FeaturePlot(data.filt, features =c("KLF2","KLF4") , max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
#macrophage  meiyou
FeaturePlot(data.filt, features =c("PDCD4") , split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
#DC
FeaturePlot(data.filt, features =c("TYMS","MPO","TPX2") , max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)
FeaturePlot(data.filt, features =t_inter[1:3] , max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)

FeaturePlot(data.filt, features =markers_all_single[which(markers_all_single$cluster==16)[1:3],"gene"] , max.cutoff = 3, cols = c("grey", "red")
            ,min.cutoff = "q4",label=T, label.size=4, pt.size=1)

markers.to.plot <- markers_all[which(markers_all$cluster==3),]
markers.to.plot=markers.to.plot[order(markers.to.plot$avg_logFC,decreasing = T),"gene"][1:5]

DotPlot(immune.combined, features =c("CD79A","CD79B"), cols = 
          c("blue", "red"), dot.scale = 8) + RotatedAxis()

plots <- VlnPlot(data.filt, features = c("CD79A","CD79B"), 
                 split.by = "orig.ident", group.by = "seurat_clusters", pt.size = 0, combine = 
                   FALSE)

save.image(file="single_cell.rda",ascii=FALSE,compress=TRUE)

cluster345_split_markers=list("YXP1_cluster3"=c(),"YXP1_cluster4"=c(),"YXP1_cluster5"=c(),
                              "YXP2_cluster3"=c(),"YXP2_cluster4"=c(),"YXP2_cluster5"=c())
for(cluster in c("YXP1","YXP2")){
  for(i in c(3,4,5)){
    print(paste("cluster",i,sep=""))
    markers <- FindMarkers(object = data.filt, group.by="orig.ident" ,subset.ident =i,ident.1 = cluster)
   # print(head(markers[order(markers$avg_logFC,decreasing = T),]))
    cluster345_split_markers[[paste(cluster,"_cluster",i,sep="")]]=markers
  }
  
}
