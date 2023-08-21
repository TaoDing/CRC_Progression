#########################Principal coordinate analysis (PCoA) on Bray–Curtis distances#########################
##########β多样性
setwd("/data/yuchen_data/cuiying/cuiying2/flat.feature.table")
library(vegan)
library(ggplot2)
#颜色选择
color=c( "#3C5488B2","#00A087B2",          
"#F39B7FB2","#91D1C2B2",         
 "#8491B4B2", "#DC0000B2",         
  "#7E6148B2","yellow",          
  "darkolivegreen1", "lightskyblue",         
   "darkgreen", "deeppink", "khaki2",          
   "firebrick", "brown1", "darkorange1",         
    "cyan1", "royalblue4", "darksalmon",         
     "darkgoldenrod1", "darkseagreen", "darkorchid") 
color1<-c("#4393c3","#f4a582","#b2182b")
color2<-c("#80cdc1","#35978f","#01665e","#003c30","#8c510a")
#以otu表为起始点
otu <- read.delim('Genus.feature-table.flat.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
bray_dis <- vegdist(otu, method = 'bray') 
#计算bray-curits距离，得到距离矩阵bray_disbray_dis
        
#以距离矩阵为起始点
pcoa <- cmdscale(bray_dis, k = (nrow(otu) - 1), eig = TRUE)
site <- data.frame(pcoa$point)[1:2]
site$name <- rownames(site)
#mapping文件，即分组信息，第一列为样本名，#其他列为分组信息。
#第一列的每一行行名及顺序要和距离矩阵保持一致。
map<-read.table("all.metadata.0522.rst1.DMM.rst1.txt",sep="\t",row.names=1,header=T)  
#将分组文件和数据文件以行名合并
merged=merge(site,map,by="row.names",all.x=TRUE)
#使用 wascores() 被动添加物种得分（坐标），丰度加权平均方法
species <- wascores(pcoa$points[,1:4], otu)
#获取前二主成分
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

#ggplot2 作图
ggplot(data = merged, aes(X1, X2)) +  
    geom_point(aes(color =Type.R_L_Site)) +  
    stat_ellipse(aes(fill = Type.R_L_Site), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +#添加置信椭圆  
    scale_color_manual(values =color[1:length(unique(map$Type.R_L_Site))]) +  
    scale_fill_manual(values =color[1:length(unique(map$Type.R_L_Site))]) +  
    theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +  
    geom_vline(xintercept = 0, color = 'gray', size = 0.5) +  
    geom_hline(yintercept = 0, color = 'gray', size = 0.5) +  
    labs(x = pcoa1, y = pcoa2)
#########################Cumulative distribution of Bray–Curtis dissimilarity#########################
library(ggplot2)
df=read.table("CRC.bray.dist.rst.meta.dedu2.txt",sep="\t",header=T)
ggplot ( df, aes ( V3 , color = V4) )  +  stat_ecdf ( geom="point")+
scale_color_manual(values=c("orange", "black", "brown"))+cuiying_theme
cuiying_theme<- theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',size=1))+
theme(axis.text.x = element_text(size = 12,color="black"),axis.text.y = element_text(size = 12,color="black")) +
theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14))+
theme(plot.margin = margin(t = 3, r = 3,  b = 3, l = 3,unit = "cm"))

#########################SourceTracker & density plot#########################
git clone https://github.com/cozygene/FEAST.git
library(FEAST)
m(list = ls())


#Example data arguments
metadata_file = "../Data_files/metadata_example.txt" #metadata file name
count_matrix = "../Data_files/otu_example.txt" #count_matrix file name
EM_iterations = 1000 # number of EM iterations. 建议使用 1000

# Load sample metadata
metadata <- read.csv("metadata_file.csv" ,header = T, sep = "\t", row.names = 1)

# Load OTU table
otus <- read.table("OTU.csv", header = T, comment = '', check = F, sep = '\t')
otus <- t(as.matrix(otus))

# Extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# Double-check that the mapping file and otu table had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}
envs <- metadata$Env
# Extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='Source')
test.ix <- which(metadata$SourceSink=='Sink')
num_sources <- length(train.ix) #number of sources
COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user
# 根据 COVERAGE 抽平
sources <- as.matrix(rarefy(otus[train.ix,], COVERAGE))
sinks <- as.matrix(rarefy(t(as.matrix(otus[test.ix,])), COVERAGE))
# Estimate source proportions for each sink
FEAST_output<-FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
Proportions_est <- FEAST_output$data_prop[,1]
names(Proportions_est) <- c(as.character(envs[train.ix]), "unknown")




library(ggpubr)
library(ggplot2)
df=read.table("CRC_source.sink.rst.txt",sep="\t",header=T)

pdf("N.density.pdf",height=7,width=8)
ggplot(dfN1,aes(x=ratio,color=Type2,fill=Type2))+geom_density(alpha = 0.5)+scale_fill_manual(
dev.off()
pdf("P.density.pdf",height=7,width=8)
ggplot(dfP1,aes(x=ratio,color=Type2,fill=Type2))+geom_density(alpha = 0.5)+scale_fill_manual(
dev.off()
pdf("T.density.pdf",height=7,width=8)
ggplot(dfT1,aes(x=ratio,color=Type2,fill=Type2))+geom_density(alpha = 0.5)+scale_fill_manual(
dev.off()


#########################Random forest#########################
set.seed(121)
library(randomForest)
library(varSelRF)
library(pROC)
load("df2.Rdata")
df=read.table("inpatient.filter.data1.sever4_08.transform_last.time.no1.last.10.txt",sep="\t",header=T,row.names=24)
df3<-cbind(df2,group=df$severe)

index <- sample(2,nrow(df2),replace = TRUE,prob=c(0.7,0.3))
traindata <- df2[index==1,]
testdata <- df2[index==2,]

rf_ntree<- randomForest(as.factor(group) ~ ., data=traindata,ntree=2000,important=TRUE,proximity=TRUE)
plot(rf_ntree)
rf<- randomForest(as.factor(group) ~ ., data=traindata,ntree=1000,important=TRUE,proximity=TRUE)
pred<-predict(rf, newdata=testdata,type="prob")
f<-cbind(as.factor(testdata$group),pred[,1])
table(pred,testdata$group)
varImpPlot(rf, n.var = min(30, nrow(rf$importance)), main = 'Top 30 - variable importance')


otu_train.cv <- replicate(20, rfcv(otu_train[-ncol(otu_train)], otu_train$group, cv.fold = 10,step = 1.5), simplify = FALSE)




x<-plot.roc(f[,1],f[,2],
            smooth=F,
            lwd=2,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=T,
            main="",
            col="red",print.auc=TRUE)
