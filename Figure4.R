#########################PICRUSt2#########################
conda create -n picrust2 -c bioconda -c conda-forge picrust2
source activate picrust2
#基于我们常用的代表序列otus.fa和特征表otutab.txt
picrust2_pipeline.py -s otus.fa -i otutab.txt -o picrust2_out 


##https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.4.1)
metagenome_pipeline.py -i ../table.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz \
-o EC_metagenome_out --strat_out

#添加EC的注释
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
#KO添加注释
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
#pathway添加注释
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#########################Differentlly expression analysis#########################
source("/public5/pulic4_backups/cuiying/03project_public4/yuchen_server_220112/yuchen_severe_20210525/cuiying/test/0116/0119/bayesreg.R")
library("multtest")
data<-read.table("KO.gene.FPKM.txt", header=TRUE,sep="\t",row.names=1,check.names=F)
data1<-data[-1,]
data2 <- as.data.frame(lapply(data1,function(x) as.numeric(as.character(x))))
rownames(data2)=rownames(data1)
rst <- bayesT(data2,25,14,doMulttest=T)
ig.edger <-rst[((rst$pVal < 0.05) &( rst$fold >2 | rst$fold < -2)),]
ig.edger1<-cbind(as.data.frame(ig.edger$fold),as.data.frame(ig.edger$BH))
row.names(ig.edger1)<-row.names(ig.edger)
write.table(rst,"T_RL.KO_25_14.DE.rst.txt",quote=F,sep="\t")
write.table(ig.edger1,"KO.gene.DE.txt",quote=F,sep="\t")
#########################KEGG pathway enrichment analysis#########################
library(MicrobiomeProfiler)
gene<-read.table("KO.gene.DE.txt",sep="\t")
ko <- enrichKO(gene$V1,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500)
df1=read.table("ko_T1_all1.txt",sep="\t",header=T)
df2=read.table("ko_T2_all1.txt",sep="\t",header=T)
df3=read.table("ko_T3_all1.txt",sep="\t",header=T)
df4=read.table("ko_T4_all1.txt",sep="\t",header=T)
pdf("T1.gene.enrichment.pdf")
ggplot(df1,aes(x = parse_ratio(GeneRatio),y = -log10(pvalue),size = Count,color = Type),alpha = 0.8) + geom_point() +theme_classic(base_size = 20) +theme(aspect.ratio = 6/6) + geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 1) +xlab('GeneRatio')+scale_color_manual(values = c('up-sig' = '#35978f','down-sig' = '#bf812d','up-none'='#969696','down-none'='#969696'))
dev.off()
pdf("T2.gene.enrichment.pdf")
ggplot(df2,aes(x = parse_ratio(GeneRatio),y = -log10(pvalue),size = Count,color = Type),alpha = 0.8) + geom_point() +theme_classic(base_size = 20) +theme(aspect.ratio = 6/6) + geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 1) +xlab('GeneRatio')+scale_color_manual(values = c('up-sig' = '#7fbc41','down-sig' = '#de77ae','up-none'='#969696','down-none'='#969696'))
dev.off()
pdf("T3.gene.enrichment.pdf")
ggplot(df3,aes(x = parse_ratio(GeneRatio),y = -log10(pvalue),size = Count,color = Type),alpha = 0.8) + geom_point() +theme_classic(base_size = 20) +theme(aspect.ratio = 6/6) + geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 1) +xlab('GeneRatio')+scale_color_manual(values = c('up-sig' = "brown",'down-sig' = '#e08214','up-none'='#969696','down-none'='#969696'))
dev.off()
pdf("T4.gene.enrichment.pdf")
ggplot(df4,aes(x = parse_ratio(GeneRatio),y = -log10(pvalue),size = Count,color = Type),alpha = 0.8) + geom_point() +theme_classic(base_size = 20) +theme(aspect.ratio = 6/6) + geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 1) +xlab('GeneRatio')+scale_color_manual(values = c('ud-sig' = '#8073ac' ,'ud-sig' ='#969696'))
dev.off()

#########################correlation coefficient plot#########################
df=read.table("cor.plot.rst",sep="\t",header=T)
library(reshape2)
library(ggsci)
dfmelt<-melt(df,id=c("X","GLYCOCAT.PWY"))
cuiying_theme<- theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',size=1))+
theme(axis.text.x = element_text(size = 12,color="black"),axis.text.y = element_text(size = 12,color="black")) +
theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14))+
theme(plot.margin = margin(t = 3, r = 3,  b = 3, l = 3,unit = "cm"))

pdf("cor.pdf",height=8,width=10)
ggplot(dfmelt, aes(x = GLYCOCAT.PWY,y =value,group=variable))+geom_point(aes(color =variable))+geom_smooth(aes(color = variable, fill = variable),method = "lm", fill = "lightgray")+scale_color_lancet()+scale_fill_lancet()+ggpubr::stat_cor(aes(color = variable),method="spearman")+cuiying_theme+xlab("Pathway")+ylab("Genus")+ylim(c(0,0.01))+xlim(c(0,0.01))
dev.off()


