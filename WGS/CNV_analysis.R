varscan_result=read.table("/home/miaoyr/CML_sample_RNA_DNA_scRNA/varscan2/gistic/varscan_segmentation.txt",header=T,sep="\t")
result_filter=varscan_result[intersect(which(varscan_result$num.mark>10),which(abs(varscan_result$seg.mean)>0.25)),]
write.table(result_filter,file="/home/miaoyr/CML_sample_RNA_DNA_scRNA/varscan2/gistic/varscan_result_filter.txt",sep="\t",quote=F)
#CNV gene file
cnv_gene=read.table("/home/miaoyr/CML_sample_RNA_DNA_scRNA/varscan2/gistic/CNV_gene_info.txt",header=T,row.names = 1,sep="\t")
list=select(org.Hs.eg.db,keys=row.names(cnv_gene),columns = c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
#cnv_gene=list$SYMBOL
#list[which(list$SYMBOL=="IKZF1"),]
list[which(list$SYMBOL=="CDKN2A"),]
list[which(list$SYMBOL=="PAX5"),]
list[which(list$SYMBOL=="ETV6"),]
sig_cnv=cnv_gene[c("ENSG00000147889","ENSG00000196092","ENSG00000139083"),]
sig_cnv=cbind(c("CDKN2A","PAX5","ETV6"),sig_cnv)
colnames(sig_cnv)[1]="Symbol"
