#########################MFuzz#########################
setwd("/public4/cuiying/03project_public4/yuchen_severe_20210525/cuiying/cuiying2/0531/01Diff_type_genus_mfuzz/RLRE-NPT/")

library('Mfuzz')

df3=read.table("all.NPT.rst.eset.df2.txt",sep="\t",header=T,row.names=1 )
df4<-as.matrix(t(df3))
df4Ex<- ExpressionSet(assayData = df4)
df4F <- standardise(df4Ex)
set.seed(2021)
cl <- mfuzz(df4F,c=8,m=1.25)
pdf("R_NPT.c8.pdf")
mfuzz.plot2(df4F, cl=cl,mfrow=c(4,4),centre=TRUE,x11=F,centre.lwd=0.2)
dev.off()
dir.create(path="mfuzz.R",recursive = TRUE)
for(i in 1:8){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("mfuzz.R","/mfuzz.R",i,".csv"))
}



dir.create(path="mfuzz.pathway",recursive = TRUE)
for(i in 1:12){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("mfuzz.pathway","/mfuzz.pathway",i,".csv"))
}


#########################Enrichment analysis of disease related species #########################
/publicgp/cuiying/03project_public6/01-project03_CRC_HZ/HOM_database/

dfPG=read.delim("pathway.Genus.txt",sep="\t",row.names=1,header=F)
dfPGT<-t(dfPG)
test2<-as.list(as.data.frame(dfPGT))
for(i in names(test2)){
if(length(which(test2[[i]]==""))>0){
d = which(test2[[i]]=="")
test2[[i]] = test2[[i]][-d]
}}
dfb=read.table("all.genus.backgroud6.txt",sep="\t".header=F)
dfb=read.table("all.genus.backgroud6.txt",sep="\t",header=F)
dfm=read.table("Tra.plot.melt3",sep="\t",header=T)
results <- hyperGtest_jimmy(as.array(test1),as.array(test2),dfm$subject,dfb$V1)
hyperGtest_jimmy <- function(GeneID2Path = GeneID2kegg_list, Path2GeneID = kegg2GeneID_list, diff_gene = sample(unique(hgu95av2_id$gene_id),
    500), universeGeneIds = unique(hgu95av2_id$gene_id)) {
    diff_gene_has_path = intersect(diff_gene, names(GeneID2Path))
    n = length(diff_gene)  #306
    N = length(universeGeneIds)  #5870
    results = c()
    for (i in names(Path2GeneID)) {
        M = length(intersect(Path2GeneID[[i]], universeGeneIds))
        # print(M)
        if (M < 5)
            next
        exp_count = n * M/N
        # print(paste(n,N,M,sep='\t'))
        k = 0
        for (j in diff_gene_has_path) {
            if (i %in% GeneID2Path[[j]])
                k = k + 1
        }
        OddsRatio = k/exp_count
        if (k == 0)
            next
        p = phyper(k - 1, M, N - M, n, lower.tail = F)
        # print(paste(i,p,OddsRatio,exp_count,k,M,sep=' '))
        results = rbind(results, c(i, p, OddsRatio, exp_count, k, M))
    }
    colnames(results) = c("PathwayID", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size")
    results = as.data.frame(results, stringsAsFactors = F)
    results$p.adjust = p.adjust(results$Pvalue, method = "BH")
    results = results[order(results$Pvalue), ]
    rownames(results) = 1:nrow(results)
    return(results)
}
results <- hyperGtest_jimmy(as.array(test1),as.array(test2),dfm$subject,dfb$V1)




library(ggplot2)
df=read.table("all.p2.txt",sep="\t",header=T)
dfT1=subset(df,Type=="T1")
dfT2=subset(df,Type=="T2")
dfT3=subset(df,Type=="T3")
dfT4=subset(df,Type=="T4")
pdf("anno.pdf",height=4.5,width=9)
ggplot(data=dfT4) +
  geom_bar(mapping = aes(x = PathwayID, y = Count/3,  fill = as.factor(PathwayID1)), stat = "identity") +
  geom_line(mapping = aes(x = PathwayID, y = -log(Pvalue,10)),color = "#4d4d4d") +
  geom_point(mapping = aes(x = PathwayID, y = -log(Pvalue,10)),color = "#4d4d4d",fill="#4d4d4d", size = 2, shape = 21) +
  scale_y_continuous(
    name = expression("-log(pvalue)"),
    sec.axis = sec_axis(~.*3 , name = "Count")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)+scale_fill_manual(values=c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#8073ac","#b2abd2","#de77ae","#7fbc41"))+theme(axis.text=element_text(angle=25,vjust=1, hjust=1))
dev.off()
