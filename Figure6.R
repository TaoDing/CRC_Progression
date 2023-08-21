#########################Distribution of adenoma types#########################
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pheatmap)
df<-read.table("a.txt",sep="\t",header=T,row.names=1)
pdf("p_RLRE.heat.pdf")
pheatmap(df,cluster_rows = F,cluster_cols = F,cellwidth = 50, 
cellheight = 50, color = colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(100),
display_numbers=T)
dev.off()
#########################Distribution of adenoma types#########################


library(ggplot2)
pdf("RLRE.same.diff.pdf",height=3.5,width=5)
ggplot(dfs1_same, aes(Var1, value, fill=Var2)) +geom_bar(position="fill",stat="identity") + guides(fill=guide_legend(reverse=F)) +scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#EDF8B1","#7FCDBB","#2C7FB8"))+theme_pander()  +theme(panel.background = element_rect(colour = "black",size = 1))
dfs1_same<-melt(as.matrix(table(dfs1$NT_DMM,dfs1$diff)))
as.matrix(table(dfs$NT_DMM,dfs$T_stage)
)
dev.off()

pdf("RLRE.same.diff.pdf",height=3.5,width=5)
ggplot(dfs1_same, aes(Var1, value, fill=Var2)) +geom_bar(position="fill",stat="identity") + guides(fill=guide_legend(reverse=F)) +scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#EDF8B1","#7FCDBB","#2C7FB8"))+theme_pander()  +theme(panel.background = element_rect(colour = "black",size = 1))
dfs1_same<-melt(as.matrix(table(dfs1$NT_DMM,dfs1$diff)))
as.matrix(table(dfs$NT_DMM,dfs$T_stage)
)
dev.off()


pdf("RLRE.same.diff.pdf",height=3.5,width=5)
ggplot(dfs1_same, aes(Var1, value, fill=Var2)) +geom_bar(position="fill",stat="identity") + guides(fill=guide_legend(reverse=F)) +scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#EDF8B1","#7FCDBB","#2C7FB8"))+theme_pander()  +theme(panel.background = element_rect(colour = "black",size = 1))
dfs1_same<-melt(as.matrix(table(dfs1$NT_DMM,dfs1$diff)))
as.matrix(table(dfs$NT_DMM,dfs$T_stage)
)
dev.off()


#########################corrplot#########################
df<-read.table("a.txt",sep="\t",header=T,row.names=1)
library(RColorBrewer)
 
corrplot(M, type = "upper", order = "hclust",col = brewer.pal(n = 8, name = "RdBu"))
corrplot(M, type = "upper", order = "hclust",col = brewer.pal(n = 8, name = "RdBu"))
corrplot(M, type = "upper", order = "hclust",col = brewer.pal(n = 8, name = "RdBu"))

######################correlation coefficient######################


import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
df1=pd.read_csv("df.rst.merge.mfuzz.txt",sep="\t",index_col=0)
df2=pd.read_csv("lefse.pathway.mfuzz.RA.rst.txt",sep="\t",index_col=0)
df1
s=df1.columns
s1=df1.columns
s2=df2.columns
s2
with open("cor.txt", 'w') as fout:
    for n1 in s1:
        for n2 in s2:
            r, p = spearmanr(df1[n1], df1[n2])
            fout.write(s1 + '\t' + s2 + '\t' + str(r) + '\t' + str(p)+ '\n')
