#BCR result
T1_BCR=read.table("/home/miaoyr/rna/CATT_BCR_result/T1_bcr.IGH.CDR3.CATT.csv",header=T,sep=",")
T2_BCR=read.table("/home/miaoyr/rna/CATT_BCR_result/T2_bcr.IGH.CDR3.CATT.csv",header=T,sep=",")
T3_BCR=read.table("/home/miaoyr/rna/CATT_BCR_result/T3_bcr.IGH.CDR3.CATT.csv",header=T,sep=",")
T4_BCR=read.table("/home/miaoyr/rna/CATT_BCR_result/T4_bcr.IGH.CDR3.CATT.csv",header=T,sep=",")

region_data=data.frame(Sample=c("T1","T2","T3","T4","T1","T2","T3","T4"),group=rep(c("Jregion","Vregion"),each=4),Count=c(length(which(T1_BCR$Jregion!="None")),
                                                         length(which(T2_BCR$Jregion!="None")),
                                                         length(which(T3_BCR$Jregion!="None")),
                                                         length(which(T4_BCR$Jregion!="None")),
                                                         length(which(T1_BCR$Vregion!="None")),
                                                         length(which(T2_BCR$Vregion!="None")),
                                                         length(which(T3_BCR$Vregion!="None")),
                                                         length(which(T4_BCR$Vregion!="None"))))
pdf("result/Transcript/BCR_distribution.pdf",4,3)
ggplot(region_data, aes(Sample, Count,fill=group,group=group)) +
  geom_bar(stat="identity",width=0.5,position = "dodge")+theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) 
dev.off()
all_Vregion=unique(c(T1_BCR$Vregion,T2_BCR$Vregion,T2_BCR$Vregion,T3_BCR$Vregion,T4_BCR$Vregion))
all_Vregion=all_Vregion[-which(all_Vregion=="None")]


v_df=data.frame()
for(i in all_Vregion){
  T1_count=length(which(T1_BCR$Vregion==i))
  T2_count=length(which(T2_BCR$Vregion==i))
  T3_count=length(which(T3_BCR$Vregion==i))
  T4_count=length(which(T4_BCR$Vregion==i))
  v_df=rbind(v_df,c("T1",T1_count,i))
  v_df=rbind(v_df,c("T2",T2_count,i))
  v_df=rbind(v_df,c("T3",T3_count,i))
  v_df=rbind(v_df,c("T4",T4_count,i))
}
colnames(v_df)=c("Sample","Count","Vchain")
v_df$Count=as.numeric(v_df$Count)
Top10=v_df[which(v_df$Count>30),]


j_df=data.frame()
for(i in all_Jregion){
  T1_count=length(which(T1_BCR$Jregion==i))
  T2_count=length(which(T2_BCR$Jregion==i))
  T3_count=length(which(T3_BCR$Jregion==i))
  T4_count=length(which(T4_BCR$Jregion==i))
  j_df=rbind(j_df,c("T1",T1_count,i))
  j_df=rbind(j_df,c("T2",T2_count,i))
  j_df=rbind(j_df,c("T3",T3_count,i))
  j_df=rbind(j_df,c("T4",T4_count,i))
}
colnames(j_df)=c("Sample","Count","Jchain")
j_df$Count=as.numeric(j_df$Count)

pdf("result/Transcript/BCR_type.pdf",6,4)
ggplot(Top10, aes(Vchain, Count,fill=Sample,group=Sample)) +
  geom_bar(stat="identity",width=0.5,position = "fill")+theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

ggplot(j_df, aes(Jchain, Count,fill=Sample,group=Sample)) +
  geom_bar(stat="identity",width=0.5,position = "fill")+theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

all_Jregion=unique(c(T1_BCR$Jregion,T2_BCR$Jregion,T2_BCR$Jregion,T3_BCR$Jregion,T4_BCR$Jregion))
all_Jregion=all_Jregion[-which(all_Jregion=="None")]
