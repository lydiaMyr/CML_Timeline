exp=read.table("SRR944626_count",header=T,sep="\t",row.names=1)
count_mat=as.matrix(exp[,6])
row.names(count_mat)=row.names(exp)
for(i in c(27,28,29,40,41)){
  count_tmp=read.table(paste("SRR9446",i,"_count",sep=""),row.names=1)
  gg=intersect(row.names(count_mat),row.names(count_tmp))
  count_mat=cbind(count_mat[gg,],count_tmp[gg,6])
}
count_mat_new=apply(count_mat,2,as.numeric)
row.names(count_mat_new)=row.names(count_mat)
colnames(count_mat_new)=paste("BC_new",seq(1,6),sep="")



exp=read.table("/home/miaoyr/rna/BC_sample6_new_TMP/SRR944626.gene_abund.tab",header=T,sep="\t")
exp_new=exp[!duplicated(exp$Gene.ID),]
row.names(exp_new)=exp_new$Gene.ID
exp_new=exp_new[,-1]
tpm_mat=as.matrix(exp_new$TPM)
row.names(tpm_mat)=row.names(exp_new)
for(i in c(27,28,29,40,41)){
  tpm_tmp=data.frame(read.table(paste("/home/miaoyr/rna/BC_sample6_new_TMP/SRR9446",i,".gene_abund.tab",sep=""),header=T,sep="\t"))
  tpm_tmp_new=tpm_tmp[!duplicated(tpm_tmp$Gene.ID),]
  row.names(tpm_tmp_new)=tpm_tmp_new$Gene.ID
  tpm_tmp_new=tpm_tmp_new[,-1]
  gg=intersect(row.names(tpm_mat),row.names(tpm_tmp_new))
  tpm_mat=cbind(tpm_mat[gg,],tpm_tmp_new[gg,"TPM"])
}
tpm_mat_new=apply(tpm_mat,2,as.numeric)
row.names(tpm_mat_new)=row.names(tpm_mat)
colnames(tpm_mat_new)=paste("BC_new",seq(1,6),sep="")
write.table(tpm_mat_new,file="BC_sample_new_exp_TPM",sep="\t",quote=F)





