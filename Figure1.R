#########################permutational multivariate analysis of variance (PERMANOVA)#########################
library(vegan)
explain <- function(a,b){a/b}
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
dis <- vegdist(otu, method = 'bray')
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE)
adonis_result_dis <- adonis(dis~site, group, permutations = 999) #根据 group$site 这一列样本分组信息进行 PERMANOVA 分析，随机置换检验 999 次
otuput <- data.frame(adonis_result_dis$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'PERMANOVA.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

pdf("PERMANOVA.rst.plot.pdf",height=10,width=9.5)
ggplot(data = dfm, aes(x = X, y = value, fill = X))+geom_bar(stat = "identity", position = "stack")+facet_grid(variable~.)+ scale_fill_manual(values = c('#941319','#4a8a53','#8B864E','#E77D13','#FFD700','#EEB4B4','LemonChiffon3','#5D478B','#696969'))+labs(x=NULL,y='Pseudotype Entry')+
  theme_test(base_size = 15)+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 60,vjust = 0.6),
        axis.text = element_text(color='black'),
        strip.background = element_rect(color = 'transparent', fill = 'transparent'))
dev.off()


#########################Genera average relative abundance and prevalent plot#########################
library(ggplot2)
df=read.table("RA.prevlance.plot",sep="\t",header=F)
df$V1<- factor(df$V1,levels=rev(unique(df$V1)))
df$V4<- factor(df$V4,levels=c("ALL","Carcinoma adjacent","Adenoma","Carcinoma"))
cuiying_theme<- theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',size=1))+
theme(axis.text.x = element_text(size = 12,color="black"),axis.text.y = element_text(size = 12,color="black")) +
 theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14))+
 theme(plot.margin = margin(t = 3, r = 3,  b = 3, l = 3,unit = "cm"))+
 theme(plot.title = element_text(size = 12, face = "bold"),legend.title=element_text(size=10), legend.text=element_text(size=9))
df$V4<- factor(df$V4,levels=c("ALL","Carcinoma adjacent","Adenoma","Carcinoma"))
pdf("RA.prevlance.pdf",height=11,width=10)
ggplot(df, aes(x = V4, y = V1)) +  
geom_point(aes(size = V3, fill = V2), alpha = 0.85, shape = 21,width=10,stroke = 1) + 
scale_fill_material("red")+scale_size_continuous(range=c(0.1,8.5))+
cuiying_theme+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))+xlab("")+ylab("")

dev.off()
