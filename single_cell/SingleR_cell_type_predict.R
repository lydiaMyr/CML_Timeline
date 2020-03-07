singler = CreateSinglerSeuratObject(data.filt@assays$RNA@counts, "YXP",
                                    min.genes = 500, technology, species = "Human",
                                    normalize.gene.length = F, min.cells = 2, npca = 17
                                    regress.out = "nUMI", reduce.seurat.object = T)
singler=CreateSinglerSeuratObject(data.filt@assays$RNA@counts, annot = NULL,  "YXP",
                          min.genes = 200, technology = "10X", species = "Human",
                          citation = "", ref.list = list(), normalize.gene.length = F,
                          variable.genes = "de", fine.tune = T, reduce.file.size = T,
                          do.signatures = F, min.cells = 2, npca = 17,
                          regress.out = "nUMI", do.main.types = T, reduce.seurat.object = T,
                          temp.dir = NULL, numCores = SingleR.numCores)
sample_lable=singler$singler[[1]]$SingleR.clusters.main$labels
for(i in seq(3,5)){
  print(paste("cluster",i,sep=""))
  cluster0=intersect(names(data.filt@active.ident[which(data.filt@active.ident==i)]),names(which(data.filt$orig.ident=="YXP1")))
  cluster0_lable=sample_lable[cluster0,]
  print(paste("Cell count:",length(cluster0),sep=""))
  a=which.max(table(cluster0_lable))
  print(a)
  
}
cluster0=names(data.filt$seurat_clusters[which(data.filt$seurat_clusters==0)])
cluster0_lable=sample_lable[cluster0,]
which.max(table(cluster0_lable))
cluster1=names(data.filt$seurat_clusters[which(data.filt$seurat_clusters==1)])
cluster1_lable=sample_lable[cluster1,]
which.max(table(cluster1_lable))



singler = SingleR(method = "single", sc_data, ref_data, types, clusters = NULL,
                  genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                  fine.tune = TRUE, fine.tune.thres = 0.05, sd.thres = 1,
                  do.pvals = T, numCores = SingleR.numCores)