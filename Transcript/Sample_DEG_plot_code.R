library("NOISeq")
data_raw=read.table("Timeline_sample_count_matrix_final",header=T,row.names=1,sep="\t")
BC_sample=read.table("BC_sample_count_matrix",header=T,row.names=1,sep="\t")
BC_sample1=read.table("/home/miaoyr/rna/BC_sample_new_count/BC_new_sample_count",header=T,row.names=1,sep="\t")
#data=data_raw[,-1]
tmp_g=intersect(row.names(data_raw),row.names(BC_sample))
com_gene=intersect(tmp_g,row.names(BC_sample1))
#BC_exp=as.matrix(apply(BC_sample,1,mean))
#BC sample 
data=cbind(data_raw[com_gene,],BC_sample[com_gene,],BC_sample1[com_gene,])
#colnames(data)[5]="BlastCrisis"
tissue=c("T1","T2","T3","T4",paste("BC",seq(1,8),sep="_"),paste("BC1",seq(1,6),sep="_"))
tissue_run=c("t1_1","t2_1","t3_1","t4_1",paste(paste("BC",seq(1,8),sep=""),1,sep="_"),paste(paste("BC1",seq(1,6),sep=""),1,sep="_"))
myfactors = data.frame(Tissue=tissue, TissueRun=tissue_run)
mydata=readData(data=data,factors=myfactors)

mynoiseq_T1_T2 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","T2"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T2_T3 <- noiseq(mydata, factor = "Tissue",conditions = c("T2","T3"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T3_T4 <- noiseq(mydata, factor = "Tissue",conditions = c("T3","T4"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")

mynoiseq_T1_T3 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","T3"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")

deg12=DEGs_nine(mynoiseq_T1_T2,"T1_VS_T2_count_diff")
deg23=DEGs_nine(mynoiseq_T2_T3,"T2_VS_T3_count_diff")
deg34=DEGs_nine(mynoiseq_T3_T4,"T3_VS_T4_count_diff")
deg13=DEGs_nine(mynoiseq_T1_T3,"T1_VS_T3_count_diff")
pdf("result/Transcript/T3_T4_deg_pathway.pdf",8,6)
enrich_fun(deg34$up)
enrich_fun(deg34$down)
dev.off()

pdf("result/Transcript/T1_T2_deg_pathway.pdf",8,6)
enrich_fun(deg12$up)
enrich_fun(deg12$down)
dev.off()

pdf("result/Transcript/T2_T3_deg_pathway.pdf",8,6)
enrich_fun(deg23$up)
enrich_fun(deg23$down)
dev.off()

pdf("result/Transcript/T1_T3_deg_pathway.pdf",8,6)
enrich_fun(deg13$up)
enrich_fun(deg13$down)
dev.off()


pheatmap_fun=function(gene){
  EG2Ensembl=toTable(org.Hs.egENSEMBL) #将ENTREZID和ENSEMBL对应的数据存入该变量
  geneLists=data.frame(ensembl_id=gene)
  results=merge(geneLists,EG2Ensembl,by='ensembl_id',all.x=T)
  id=na.omit(results$gene_id)  #提取出非NA的ENTREZID
  entrez_id=as.character(id)
  # ego <- enrichGO(
  #   OrgDb="org.Hs.eg.db", gene = entrez_id, ont = "CC", pvalueCutoff = 0.15,pAdjustMethod = "BH",
  #   qvalueCutoff = 0.15, readable= TRUE)
  # GO_plot=dotplot(ego,showCategory=20,title="Enrichment GO Top20"
  ekegg<-enrichKEGG(gene = entrez_id,organism = "hsa",keyType = "kegg",  pvalueCutoff = 0.15,pAdjustMethod = "BH",
                    qvalueCutoff = 0.15,minGSSize = 5,
                    use_internal_data = FALSE)
  kegg_plot=dotplot(ekegg,showCategory=20,title="Enrichment KEGG Top20")
  kegg_pathway=data.frame(ekegg)
  return(kegg_pathway)
}
kegg_pathway=pheatmap_fun(up_final1)
gene=as.vector(unlist(strsplit(kegg_pathway["hsa04068","geneID"],"\\/")))
sample_pheatmap(gene,-c(2,3,4))
sample_pheatmap=function(FOXO_gene,sample){
  FOXO_symbol = bitr(FOXO_gene, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  FOXO_ensembl=bitr(FOXO_gene, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
  #all_exp[FOXO_ensembl$ENSEMBL,][-4,]
  g=intersect(row.names(all_exp),FOXO_ensembl$ENSEMBL)
  exp_symbol=all_exp[g,]
  row.names(exp_symbol)=exp_symbol$Symbol
  exp_symbol=exp_symbol[-11,-1]
  p=pheatmap(exp_symbol[,sample],scale="row", cluster_rows=F, cluster_cols=F )
  print(p)
}

# enrich_fun(deg12$up)
# enrich_fun(deg23$up)
# enrich_fun(deg34$up)
# 
# 
# enrich_fun(deg12$down)
# enrich_fun(deg23$down)
# enrich_fun(deg34$down)


a=intersect(deg12$up,deg23$down)
bad_factor=intersect(a,deg34$up)
b=intersect(deg12$down,deg23$up)
good_factor=intersect(b,deg34$down)

# enrich_fun(good_factor)
# enrich_fun(good_factor)

list=select(org.Hs.eg.db,keys=row.names(noiseq_up),columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
gene=row.names(noiseq_up[order(noiseq_up$log2FC,decreasing = T),])[1:10]
list[which(list$ENSEMBL==gene[8]),]
DEGs_nine=function(mynoiseq,file_name){
  noiseq_result<-mynoiseq@results[[1]]
  noiseq_degene<-subset(noiseq_result,noiseq_result$prob>=0.99)
  log2FC=log(noiseq_degene[,1]/noiseq_degene[,2],2)
  #Symbol=as.vector(unlist(data_raw[row.names(noiseq_degene),"Gene.Name"]))
  noiseq_degene_new=cbind(noiseq_degene,log2FC)#,Symbol)
  noiseq_up<-subset(noiseq_degene_new,noiseq_degene_new$log2FC>=1)
  noiseq_down<-subset(noiseq_degene_new,noiseq_degene_new$log2FC<=(-1))
  deg=list(up=row.names(noiseq_up),down=row.names(noiseq_down))
  return(deg)
  # write.table(noiseq_up,file=paste(file_name,"up_genes.txt",sep="_"),sep="\t",quote=F)
  # write.table(noiseq_down,file=paste(file_name,"down_genes.txt",sep="_"),sep="\t",quote=F)
}

mynoiseq_T1_BC1 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_1"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC2 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_2"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC3 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_3"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC4 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_4"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC5 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_5"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC6 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_6"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC7 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_7"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC8 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC_8"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")

mynoiseq_T1_BC1_1 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC1_1"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC1_2 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC1_2"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC1_3 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC1_3"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC1_4 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC1_4"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC1_5 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC1_5"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")
mynoiseq_T1_BC1_6 <- noiseq(mydata, factor = "Tissue",conditions = c("T1","BC1_6"), k = NULL, norm = "rpkm", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "no")




deg1=DEGs(mynoiseq_T1_BC1,"T1_VS_BC1_count_diff")
deg2=DEGs(mynoiseq_T1_BC2,"T1_VS_BC2_count_diff")
deg3=DEGs(mynoiseq_T1_BC3,"T1_VS_BC3_count_diff")
deg4=DEGs(mynoiseq_T1_BC4,"T1_VS_BC4_count_diff")
deg5=DEGs(mynoiseq_T1_BC5,"T1_VS_BC5_count_diff")
deg6=DEGs(mynoiseq_T1_BC6,"T1_VS_BC6_count_diff")
deg7=DEGs(mynoiseq_T1_BC7,"T1_VS_BC7_count_diff")
deg8=DEGs(mynoiseq_T1_BC8,"T1_VS_BC8_count_diff")

#new bc_sample
deg1_1=DEGs(mynoiseq_T1_BC1_1,"T1_VS_BC1_1_count_diff")
deg1_2=DEGs(mynoiseq_T1_BC1_2,"T1_VS_BC1_2_count_diff")
deg1_3=DEGs(mynoiseq_T1_BC1_3,"T1_VS_BC1_3_count_diff")
deg1_4=DEGs(mynoiseq_T1_BC1_4,"T1_VS_BC1_4_count_diff")
deg1_5=DEGs(mynoiseq_T1_BC1_5,"T1_VS_BC1_5_count_diff")
deg1_6=DEGs(mynoiseq_T1_BC1_6,"T1_VS_BC1_6_count_diff")


deg_ls=list(deg1=deg1,deg1_1=deg1_1,deg2=deg2,deg1_2=deg1_2,deg1_3=deg1_3,deg1_4=deg1_4,deg1_5=deg1_5,deg1_6=deg1_6,deg3=deg3,deg4=deg4,deg5=deg5,deg6=deg6,deg7=deg7,deg8=deg8)
up_gene=unique(c(deg1$up,deg1_1$up,deg2$up,deg1_2$up,deg3$up,deg1_3$up,deg4$up,deg1_4$up,deg5$up,deg1_5$up,deg6$up,deg1_6$up,deg7$up,deg8$up))
down_gene=unique(c(deg1$down,deg1_1$down,deg2$down,deg1_2$down,deg3$down,deg1_3$down,deg4$down,deg1_4$down,deg5$down,deg1_5$down,deg6$down,deg1_6$down,deg7$down,deg8$down))

up_gene_fre=deg_intersect(up_gene,"up")
down_gene_fre=deg_intersect(down_gene,"down")
up_final=names(which(up_gene_fre>5))
down_final=names(which(down_gene_fre>5))

up_kegg=enrich_fun(up_final)
down_kegg=enrich_fun(down_final)

#8组交集
up_final1=names(which(up_gene_fre>13))
down_final1=names(which(down_gene_fre>13))

# up_kegg_8=enrich_fun(up_final1)
# down_kegg_8=enrich_fun(down_final1)
# pdf("up381_down233_enrichment_pathway.pdf",8,6)
# print(up_kegg_8)
# print(down_kegg_8)
# dev.off()
up_kegg_8_717=enrich_fun(up_final1)
down_kegg_8_574=enrich_fun(down_final1)
pdf("result/Transcript/T1_NCBI_BC_enrichment_pathway_2_24.pdf",7,5)
print(up_kegg_8_717)
print(down_kegg_8_574)
dev.off()

CML_tmp=read.table("Timeline_sample_TPM_all_sample",header=T,row.names=1,sep="\t")
BC_sample_TPM=read.table("BC_sample_TPM",header=T,row.names=1,sep="\t")
BC_new_sample_TPM=read.table("BC_sample_new_exp_TPM",header=T,row.names=1,sep="\t")
com_gene=intersect(row.names(BC_sample),row.names(CML_tmp))
com_gene=intersect(com_gene,row.names(BC_new_sample_TPM))
all_exp=cbind(CML_tmp[com_gene,],BC_sample_TPM[com_gene,],BC_new_sample_TPM[com_gene,])[,-6]
down_gene_exp=all_exp[down_gene,-c(3,4,5)]



deg_intersect=function(up_gene,deg_type){
  gene_count=c()
  for(gene in up_gene){
    count=0
    for(de in names(deg_ls)){
      deg_tmp=deg_ls[[de]]
      up_tmp=deg_tmp[[deg_type]]
      if(gene %in% up_tmp){
        count=count+1
      }
    }
    gene_count=c(gene_count,count)
  }
  names(gene_count)=up_gene
  return(gene_count)
}




library("org.Hs.eg.db")
enrich_fun=function(gene){
  EG2Ensembl=toTable(org.Hs.egENSEMBL) #将ENTREZID和ENSEMBL对应的数据存入该变量
  geneLists=data.frame(ensembl_id=gene)
  results=merge(geneLists,EG2Ensembl,by='ensembl_id',all.x=T)
  id=na.omit(results$gene_id)  #提取出非NA的ENTREZID
  entrez_id=as.character(id)
  ego <- enrichGO(
    OrgDb="org.Hs.eg.db", gene = entrez_id, ont = "CC", pvalueCutoff = 0.15,pAdjustMethod = "BH",
    qvalueCutoff = 0.15, readable= TRUE)
  GO_plot=dotplot(ego,showCategory=20,title="Enrichment GO Top20"
  ekegg<-enrichKEGG(gene = entrez_id,organism = "hsa",keyType = "kegg",  pvalueCutoff = 0.15,pAdjustMethod = "BH",
                    qvalueCutoff = 0.15,minGSSize = 5,
                    use_internal_data = FALSE)
  kegg_plot=dotplot(ekegg,showCategory=20,title="Enrichment KEGG Top20")
  return(kegg_plot)
}
#FOXO pathway
library(pheatmap)
kegg_pathway=data.frame(ekegg)
#FOXO_gene_raw=c("1027","604","10000","6794","10365","4193","26260","10018","5896","3575","5897","4303","8660")
#FOXO_gene_raw=c("208","1027","1647","10365","4193","5295","472","5896","3575","11337","5897","8660")
hema_gene=c("928","3655","1791","2322","3559","947","930","4311","100133941")
# FOXO_gene_raw=c("208","1027","604","1647","10365","4193","5295","472","5896","3575","11337","5897","8660")
exp_heatmap_func(hema_gene)
Immuno_deficiency=c("115650","5896","3575","5897","916","6891")
#Transcriptional_misregulation=c("7185","5218","6688","6929","3560","3002","2321","84444","1027","604","4193","55589","5079","2119","8148")
Transcriptional_misregulation=c("6688","6929","2521","1027","604","1647","4193","472","2078","4299","5079","8148")
Transcriptional_misregulation_new=c("4353","330","3560","3002","4318","860","64919","597","2209","3576","3684","929","4300","4094","6256","1991","1668")
exp_heatmap_func(Transcriptional_misregulation_new)
platinum_drug_resistance=c("7155","208","2956","4193","5295","472")
apoptosis=c("1522","208","1647","8739","4170","5295","472","3725","5551")
apoptosis_new=c("330","355","3710","1439","3002","1512","8797","1075","356","1509","597","1521","5551","4217")
exp_heatmap_func(apoptosis)
exp_heatmap_func=function(FOXO_gene){
  FOXO_symbol = bitr(FOXO_gene, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  FOXO_ensembl=bitr(FOXO_gene, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
  #all_exp[FOXO_ensembl$ENSEMBL,][-4,]
  g=intersect(row.names(all_exp),FOXO_ensembl$ENSEMBL)
  exp_symbol=all_exp[g,]
  row.names(exp_symbol)=exp_symbol$Symbol
  exp_symbol=exp_symbol[,-1]
  pheatmap(exp_symbol,scale="row")
}
library(pheatmap)
pheatmap(exp_symbol,scale="row")
DEGs=function(mynoiseq,file_name){
  noiseq_result<-mynoiseq@results[[1]]
  noiseq_degene<-subset(noiseq_result,noiseq_result$prob>=0.9)
  log2FC=log(noiseq_degene[,1]/noiseq_degene[,2],2)
  #Symbol=as.vector(unlist(data_raw[row.names(noiseq_degene),"Gene.Name"]))
  noiseq_degene_new=cbind(noiseq_degene,log2FC)#,Symbol)
  noiseq_up<-subset(noiseq_degene_new,noiseq_degene_new$log2FC>=1)
  noiseq_down<-subset(noiseq_degene_new,noiseq_degene_new$log2FC<=(-1))
  deg=list(up=row.names(noiseq_up),down=row.names(noiseq_down))
  return(deg)
  # write.table(noiseq_up,file=paste(file_name,"up_genes.txt",sep="_"),sep="\t",quote=F)
  # write.table(noiseq_down,file=paste(file_name,"down_genes.txt",sep="_"),sep="\t",quote=F)
}


exp_symbol=all_exp[g,]
row.names(exp_symbol)=exp_symbol$Symbol
exp_symbol=exp_symbol[,-1]
#exp_symbol=exp_symbol[-10,]
pheatmap(exp_symbol,scale="row")
