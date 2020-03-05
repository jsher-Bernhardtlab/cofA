library("data.table")
Location_of_Coding_folder='C:/cygwin64/home/Bernhardt'
setwd(Location_of_Coding_folder)

#Get a list of every data set
Folders=list.dirs(recursive=FALSE)
#Take the restult file from the first sample
SAMPLE=Folders[1]
SAMPLE_results_fitness=fread(paste(SAMPLE,"/",SAMPLE,"_Fitness_results_1gABC.txt",sep=""))
Results=SAMPLE_results_fitness[,c(1,2,3,16)]
colnames(Results)[dim(Results)[2]]=SAMPLE

#Add the results from all of the other samples
for (i in 2:length(Folders)){
  SAMPLE=Folders[i]
  if (file.exists(paste(SAMPLE,"/",SAMPLE,"_Fitness_results_1gABC.txt",sep=""))){
  SAMPLE_results_fitness=fread(paste(SAMPLE,"/",SAMPLE,"_Fitness_results_1gABC.txt",sep=""))
  Results=cbind(Results,SAMPLE_results_fitness[,16])
  colnames(Results)[dim(Results)[2]]=SAMPLE
}
}

#Clean up data a little, remove essential genes
DATA=Results[,-(1:2)]
DATA=as.matrix(DATA)
rownames(DATA) <- unlist(Results[,1])
DATA = DATA[DATA[,1]!="E",]
DATA = DATA[DATA[,1]!="R",]
DATA = DATA[,-1]
DATA_all=mapply(DATA, FUN=as.numeric)
DATA_all=matrix(DATA_all,ncol=dim(DATA)[2], nrow=dim(DATA)[1])
rownames(DATA_all) <- rownames(DATA)
colnames(DATA_all) <- colnames(DATA)

#Change the order of the samples
DATA_all=DATA_all[,c(1,2,
                     3,4,29,10,11,12,20,21,5,39,40,13,30,18,19,25,26,
                     14,6,7,
                     27,28,
                     31,32,
                     8,9,15,16,35,36,37,38,
                     33,34,17,22,23,24)]

#Make heatmap
colors = c(seq(0,0.8,length=100),seq(0.81,1.20,length=100),seq(1.21,2,length=100))
my_palette <- colorRampPalette(c("red4", "firebrick1","white", "royalblue","blue4"))(n = 299)
library("gplots")
setEPS()
png("Heatmap_1gABC.png", width = 1600, height = 800)
par(oma=c(5,0,0,0))
DATA_all_Heatmap=heatmap.2(DATA_all,
                           Colv=FALSE,
                           na.rm=TRUE, scale="none", breaks=colors, col=my_palette,
                           density.info="none", trace="none",
                           dendrogram="none",
                           symm=F,symkey=F,symbreaks=T,cexCol=2.5,
                           labRow = FALSE, margins=c(7,2),
                           key=FALSE,
                           lwid=c(1,4),lhei=c(1.5,7),
                           cexRow=2)
Heatmap=DATA_all[rev(DATA_all_Heatmap$rowInd),DATA_all_Heatmap$colInd]
dev.off()