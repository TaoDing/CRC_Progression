#########################interaction of genera##################################################
library(ggplot2)
df=read.table("genus.N-RLRE.txt",sep="\t",header=F)
pdf("Normal.point.pdf",width=11,height=7)
ggplot(data = df, mapping = aes(x = V1, y = V3,color=V5)) + 
geom_jitter(alpha=0.6,width=0.3)+ geom_hline(aes(yintercept=0.5), colour="#990000", linetype="dashed")+
geom_hline(aes(yintercept=-0.5), colour="#990000", linetype="dashed")+
theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))+scale_color_manual(values=c("#5ab4ac","#01665e","#8c510a"))+cuiying_theme+xlab("")+ylab("correlation coefficient")
dev.off()

df2=read.table("N.bar.plot",sep="\t")
pdf("Normal.barplot.pdf",width=11,height=5)
ggplot(data = df2, mapping = aes(x = V2, y = V1,color=V3,fill=V3)) +
geom_bar(stat="identity", position = 'dodge')+theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))+
scale_color_manual(values=c("#5ab4ac","#01665e","#8c510a"))+cuiying_theme+xlab("")+
ylab("correlation coefficient")+scale_fill_manual(values=c("#5ab4ac","#01665e","#8c510a"))
dev.off()


#########################biofilm#########################
library(ggplot2)
df<-read.table("Kegg.tumor.rst.type1.txt",sep="\t",header=T)
df$qvalue<0.01
df$Description<-factor(df$Description,df$Description)
cuiying_theme<- theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',size=1))+theme(axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black")) + theme(axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))+theme(plot.margin = margin(t = 3, r = 3,  b = 3, l = 3,unit = "cm"))+theme(plot.title = element_text(size = 16, face = "bold"),legend.title=element_text(size=16), legend.text=element_text(size=16))

library(ggplot2)
ggplot(data = mtcars,
       aes(x = Foldchange, y = -log(pavlue), 
       colour = factor(cyl))) +
       geom_point() +
       coord_cartesian() +
       facet_wrap(~ cyl) +
       theme_bw()
pdf("fujitu0606-3.pdf",height=7,width=11)
ggplot(df, aes(x=-log(p.adjust),y=Description,fill=-log(pvalue),na.rm = FALSE))+geom_bar(stat="identity",na.rm = FALSE)+geom_text(aes(x=-log(p.adjust),y=Description,label=Count))+scale_fill_distiller(palette="BuPu",direction = 2)+cuiying_theme
dev.off()
E:\博士学习\崔莹-sysu-2022假期工作\03-CRC\02-文章撰写\CRC文章提交\CRC绘图汇总\F3
#########################Genus boxplot#########################
library(ggplot2)
library(ggpubr)
df=read.table("all.genus.Type.R_L_Site.txt",sep="\t",header=T)
head(df)
dfN=subset(df,Type=="Normal")
dfT=subset(df,Type=="Tumor")
dim(dfN)
dim(dfT)
my_comparisons1 <- list( c("Rectum", "Left colon"), c("Rectum", "Right colon"), c("Left colon", "Right colon"))
cuiying_theme<- theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',size=1))+theme(axis.text.x = element_text(size = 12,color="black"),axis.text.y = element_text(size = 12,color="black")) + theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14))+theme(plot.margin = margin(t = 3, r = 3,  b = 3, l = 3,unit = "cm"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.title=element_text(size=10), legend.text=element_text(size=9))
library(ggplot2)
library(reshape2)
library(ggpubr)
pdf("Carcinoma-g__Bacteroides.pdf",width=8,height=7)
ggplot(dfT,aes(y=g__Bacteroides,x=R_L_Site,fill=R_L_Site))+stat_boxplot(geom='errorbar', linetype=1, width=0.5)+geom_boxplot()+stat_compare_means(comparisons = my_comparisons1,method = "wilcox.test")+cuiying_theme+scale_fill_manual(values=c("#8c510a","#5ab4ac","#01665e"))+ labs(title = "Carcinoma:g__Bacteroides")
dev.off()
pdf("Carcinoma-g__Bacteroides.pdf",width=9,height=7)
ggplot(dfT,aes(y=g__Bacteroides,x=R_L_Site,fill=R_L_Site))+stat_boxplot(geom='errorbar', linetype=1, width=0.5)+geom_boxplot()+stat_compare_means(comparisons = my_comparisons1,method = "wilcox.test")+cuiying_theme+scale_fill_manual(values=c("#8c510a","#5ab4ac","#01665e"))+ labs(title = "Carcinoma:g__Bacteroides")
dev.off()
pdf("Carcinoma-g__Bacteroides.pdf",width=8.7,height=7)
ggplot(dfT,aes(y=g__Bacteroides,x=R_L_Site,fill=R_L_Site))+stat_boxplot(geom='errorbar', linetype=1, width=0.5)+geom_boxplot()+stat_compare_means(comparisons = my_comparisons1,method = "wilcox.test")+cuiying_theme+scale_fill_manual(values=c("#8c510a","#5ab4ac","#01665e"))+ labs(title = "Carcinoma:g__Bacteroides")
dev.off()
pdf("Carcinoma-g__Leptotrichia.pdf",width=9,height=7)
ggplot(dfT,aes(y=g__Leptotrichia,x=R_L_Site,fill=R_L_Site))+stat_boxplot(geom='errorbar', linetype=1, width=0.5)+geom_boxplot()+stat_compare_means(comparisons = my_comparisons1,method = "wilcox.test")+cuiying_theme+scale_fill_manual(values=c("#8c510a","#5ab4ac","#01665e"))+ labs(title = "Carcinoma:g__Leptotrichia")
dev.off()
pdf("CarcinomaA-g__Bacteroides.pdf",width=9,height=7)
ggplot(dfN,aes(y=g__Bacteroides,x=R_L_Site,fill=R_L_Site))+stat_boxplot(geom='errorbar', linetype=1, width=0.5)+geom_boxplot()+stat_compare_means(comparisons = my_comparisons1,method = "wilcox.test")+cuiying_theme+scale_fill_manual(values=c("#8c510a","#5ab4ac","#01665e"))+ labs(title = "Carcinoma adjacent
:g__Bacteroides")
dev.off()
pdf("Carcinoma-g__Leptotrichia.pdf",width=9,height=7)
ggplot(dfT,aes(y=g__Leptotrichia,x=R_L_Site,fill=R_L_Site))+stat_boxplot(geom='errorbar', linetype=1, width=0.5)+geom_boxplot()+stat_compare_means(comparisons = my_comparisons1,method = "t.test")+cuiying_theme+scale_fill_manual(values=c("#8c510a","#5ab4ac","#01665e"))+ labs(title = "Carcinoma:g__Leptotrichia")
dev.off()

