T1_mut_gene=read.table("T1.vep.vcf.split.sig.mut.gene")
T2_mut_gene=read.table("T2.vep.vcf.split.sig.mut.gene")
T3_mut_gene=read.table("T3.vep.vcf.split.sig.mut.gene")
T4_mut_gene=read.table("T4.vep.vcf.split.sig.mut.gene")

a=setdiff(T2_mut_gene$V1,T1_mut_gene$V1)
b=setdiff(a,T3_mut_gene$V1)

TMB_df=data.frame(Sample=c("T1","T2","T3","T4"),Count=c(nrow(T1_mut_gene),nrow(T2_mut_gene),
                                                        nrow(T3_mut_gene),nrow(T4_mut_gene)))
#TMB_df$Count=TMB_df$Count
pdf("result/WGS/Mutation_count.pdf",4,3)
ggplot(TMB_df, aes(Sample, Count, fill=Sample)) +
  geom_bar(stat="identity",width=0.3)+theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  geom_text(aes(x=Sample,y=Count+10,label=Count))
dev.off()
drive_mut=read.table("/home/miaoyr/CML_sample_RNA_DNA_scRNA/mutation_download_tab.txt",header=T,sep="\t")
AML_driver=drive_mut[which(drive_mut$cancer_project=="Acute_myeloid_leukemia(TCGA,US)"),]

AML_driver_gene=c()
for(i in AML_driver$driver_gene){
  AML_driver_gene=c(AML_driver_gene,as.vector(unlist(strsplit(i,", "))))
}
AML_driver_gene=unique(AML_driver_gene)
T1_driver=intersect(T1_mut_gene$V1,AML_driver_gene)
T2_driver=intersect(T2_mut_gene$V1,AML_driver_gene)
T3_driver=intersect(T3_mut_gene$V1,AML_driver_gene)
T4_driver=intersect(T4_mut_gene$V1,AML_driver_gene)

#T1 single mut
a=setdiff(T1_mut_gene$V1,T2_mut_gene$V1)
#b=intersect(a, setdiff(T1_mut_gene$V1,T3_mut_gene$V1))
c=intersect(a, setdiff(T1_mut_gene$V1,T4_mut_gene$V1))

#T3 single mut
a=setdiff(T3_mut_gene$V1,T1_mut_gene$V1)
b=intersect(a, setdiff(T3_mut_gene$V1,T2_mut_gene$V1))
c=intersect(b, setdiff(T3_mut_gene$V1,T4_mut_gene$V1))


onco_mut=c("RUNX1","ASXL1","IKZF1","WT1","TET2"," IDH1"," NRAS"," KRAS"," CBL","TP53")
a=c("CYP1A2","CYP2A6","CYP2C19","CYP2C9","CYP2D6","CYP2E1","CYP3A4","CYP3A5","GSTM1","GSTP1","GSTT1")