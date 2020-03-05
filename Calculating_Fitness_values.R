args = commandArgs(trailingOnly=TRUE)
Location_of_Coding_folder=args[1]
GENELIST_file=args[2]
Blank_results="Blank_Fitness_1gABC_v2.csv"
SAMPLE_name=args[3]
Control_name="1gABC"

library("data.table")
setwd(Location_of_Coding_folder)
SAMPLE=fread(paste(SAMPLE_name,"/",SAMPLE_name,"-trimmed-bowtieMap-TA.txt",sep=""), sep='\t')

#Get the gene names
Names_wdirections=fread(paste("Codes/",GENELIST_file,sep=""), sep='\t',quote="\"")
Names=Names_wdirections[,c(1,2,3,4)]
Names=lapply(Names,as.character)
Names$V5=as.numeric(Names$V4)-as.numeric(Names$V3)
Directions=Names_wdirections[,c(1,2,5)]
Directions=Directions[Directions$V5!="",]
Directions=droplevels(Directions)

#Remove intergenic regions
SAMPLE=SAMPLE[-grep("IG_",SAMPLE$V4),]

#Calculate sample total for proportion
SAMPLE_sum=sum(SAMPLE[,c(5,6)])

#Make the Results file
SAMPLE_results=fread(paste("Codes/", Blank_results,sep=""), sep=',',quote="\"",
                     colClasses=c(rep("character",3), rep("numeric",13)))
for (i in 1:dim(SAMPLE_results)[1]){
  Temp=unique(SAMPLE$V4)[i]
  Data_temp=SAMPLE[SAMPLE$V4==Temp,]
  if (dim(Data_temp)[1]>60){
  Data_temp=Data_temp[-c(1:25),]
  Data_temp=head(Data_temp,-25)
  }
  #1Loci
  #2Name
  #3Essential
  #4Length
  #5TA
  SAMPLE_results[i,6]=colSums(Data_temp[,5]) #Sample F and R reads
  SAMPLE_results[i,7]=colSums(Data_temp[,6]) #Sample F and R reads
  SAMPLE_results[i,8]=sum(Data_temp[,c(5,6)]) #Sample Total reads
  #9Lib F
  #10Lib R
  #11Lib Total
  SAMPLE_results[i,12]=sum(Data_temp$V5>0)+sum(Data_temp$V6>0)   #Unique insertions Sample
  #13Unique Lib
  SAMPLE_results[i,14]=SAMPLE_results[i,8]/SAMPLE_sum #Sample Proportion reads
  #15Lib Proportion
  if(SAMPLE_results[i,14]==0 & SAMPLE_results[i,15]!=0) {
    SAMPLE_results[i,16]=0
  } else if(SAMPLE_results[i,14]==0 & SAMPLE_results[i,15]==0) {
      SAMPLE_results[i,16]=0
  } else {SAMPLE_results[i,16]=log(SAMPLE_results[i,14]*2^10/SAMPLE_results[i,15])/
                        log((1-SAMPLE_results[i,14])*2^10/(1-SAMPLE_results[i,15]))}
}
SAMPLE_results[is.infinite(unlist(SAMPLE_results[,16])),16]=2
SAMPLE_results[unlist(SAMPLE_results[,16])<0,16]=0

colnames(SAMPLE_results) = c('Loci', 'Name','Essential','length',"TA sites",
                             paste(SAMPLE_name, "F reads"),paste(SAMPLE_name, "R reads"),
                             paste(SAMPLE_name, "Total reads"),
                             paste(Control_name,"F reads"),paste(Control_name,"R reads"),
                             paste(Control_name,"Total reads"),
                             paste(SAMPLE_name,"Unique Insertions"),paste(Control_name,"Unique Insertions"),
                             paste(SAMPLE_name,"Proportion"),paste(Control_name,"Proportion"),
                             paste(SAMPLE_name,"Fitness")
)

fwrite(SAMPLE_results,paste(SAMPLE_name,"/",SAMPLE_name,"_Fitness_results_",Control_name,".txt",sep=""),
       sep='\t', row.names=FALSE, col.names = TRUE, quote=FALSE)