## Purpose of this script
##      This script is used to process the raw proteome data supplied by Sergo, generate the proteome data set
## =============================================

###################################
##### FIRST CLEAN WORKSPACE   #####
###################################
rm(list=ls())

###################################
###  Help functions             ###
###################################
## Detecting operating system type
get_os<-function(){
  sysinf<-Sys.info()
  if(!is.null(sysinf)){
    os<-sysinf['sysname']
    if(os=='Darwin')
      os<-"osx"
  }else{ ## mystery machine
    os<-.Platform$OS.type
    if(grep("darwin",R.version$os))
      os<-"osx"
    if(grep("linux-gnu",R.version$os))
      os<-"linux"
  }
  tolower(os)
}

#################################
### Setting up working env   ####
#################################
os_type<-get_os()
if(os_type=="linux"){
  work_dir<-c("/home/jyxia/Work/Proteome paper/DataAnalysis/RScript")
}else{
  work_dir<-c("F:/360syn/Science/Researches/OmicsResearches/ChalmersProject/Proteome paper/DataAnalysis/RScript")
}
setwd(work_dir)

################################
### Useful variables         ###
################################
dateLabel=paste(format(Sys.time(),"%Y%m%d"),"xlsx",sep = ".")

#Load file containing proteome raw data
#读入Sergo提供的原始数据，包括蛋白组数据及IS数据
dataSource<-read.csv("../SergoSuppliedData/ProteinGroupsS288C.csv",header = TRUE,sep = ";",stringsAsFactors = FALSE)
#Load file containing IS concentration data 
dataIS<-read.csv("../SergoSuppliedData/IS-S288C.csv",header = TRUE,sep = ";",stringsAsFactors = FALSE)
dataTotalProteinPercent<-read.csv("../TotalProteinConcMeasurement/TotalProteinMeasured.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)

#Here we only selected the normalized H/L ratio
#我们仅选出Ratio.H.L.normalized.S**那些列，这里共有59个样品所以有59�?
allsamples<-as.character(1:59)
allsamples<-paste("Ratio.H.L.normalized.S",allsamples,sep = "")
columnsWanted<-c("Majority.protein.IDs",allsamples)
dataRelative<-subset(dataSource,select = columnsWanted)

### added to extract information of proteins idendified in individual condition
detectedProteinNum<-apply(dataRelative,2,function(x) length(na.omit(x)))
s1<-which(!is.na(dataRelative$Ratio.H.L.normalized.S1))
s2<-which(!is.na(dataRelative$Ratio.H.L.normalized.S2))
s3<-which(!is.na(dataRelative$Ratio.H.L.normalized.S3))
s4<-which(!is.na(dataRelative$Ratio.H.L.normalized.S4))
s5<-which(!is.na(dataRelative$Ratio.H.L.normalized.S5))
s6<-which(!is.na(dataRelative$Ratio.H.L.normalized.S6))
s7<-which(!is.na(dataRelative$Ratio.H.L.normalized.S7))
s8<-which(!is.na(dataRelative$Ratio.H.L.normalized.S8))
s9<-which(!is.na(dataRelative$Ratio.H.L.normalized.S9))
s10<-which(!is.na(dataRelative$Ratio.H.L.normalized.S10))
s11<-which(!is.na(dataRelative$Ratio.H.L.normalized.S11))
s12<-which(!is.na(dataRelative$Ratio.H.L.normalized.S12))
s13<-which(!is.na(dataRelative$Ratio.H.L.normalized.S13))
s14<-which(!is.na(dataRelative$Ratio.H.L.normalized.S14))
s15<-which(!is.na(dataRelative$Ratio.H.L.normalized.S15))
s16<-which(!is.na(dataRelative$Ratio.H.L.normalized.S16))
s17<-which(!is.na(dataRelative$Ratio.H.L.normalized.S17))
s18<-which(!is.na(dataRelative$Ratio.H.L.normalized.S18))
s19<-which(!is.na(dataRelative$Ratio.H.L.normalized.S19))
s20<-which(!is.na(dataRelative$Ratio.H.L.normalized.S20))
s21<-which(!is.na(dataRelative$Ratio.H.L.normalized.S21))
s22<-which(!is.na(dataRelative$Ratio.H.L.normalized.S22))
s23<-which(!is.na(dataRelative$Ratio.H.L.normalized.S23))
s24<-which(!is.na(dataRelative$Ratio.H.L.normalized.S24))
s25<-which(!is.na(dataRelative$Ratio.H.L.normalized.S25))
s26<-which(!is.na(dataRelative$Ratio.H.L.normalized.S26))
s27<-which(!is.na(dataRelative$Ratio.H.L.normalized.S27))
d025<-intersect(intersect(s1,s2),s3)
d05<-intersect(intersect(s4,s5),s6)
d10<-intersect(intersect(s7,s8),s9)
d15<-intersect(intersect(s10,s11),s12)
d20<-intersect(intersect(s13,s14),s15)
d25<-intersect(intersect(s16,s17),s18)
d30<-intersect(intersect(s19,s20),s21)
d35<-intersect(intersect(s22,s23),s24)
d40<-intersect(intersect(s25,s26),s27)

#added on 01/29/2019 for get information of T1 samples
# especially to see how long it is needed to generate enough ribosome proteins
# when cell growth rate shift from 0.25 to 0.35 h-1
s50<-which(!is.na(dataRelative$Ratio.H.L.normalized.S50))
s51<-which(!is.na(dataRelative$Ratio.H.L.normalized.S51))
s52<-which(!is.na(dataRelative$Ratio.H.L.normalized.S52))
s53<-which(!is.na(dataRelative$Ratio.H.L.normalized.S53))
s54<-which(!is.na(dataRelative$Ratio.H.L.normalized.S54))
s55<-which(!is.na(dataRelative$Ratio.H.L.normalized.S55))
s56<-which(!is.na(dataRelative$Ratio.H.L.normalized.S56))
s57<-which(!is.na(dataRelative$Ratio.H.L.normalized.S57))
s58<-which(!is.na(dataRelative$Ratio.H.L.normalized.S58))
s59<-which(!is.na(dataRelative$Ratio.H.L.normalized.S59))
T11<-intersect(s50,s51)
T12<-intersect(s52,s53)
T13<-intersect(s54,s55)
T14<-intersect(s56,s57)
T15<-intersect(s58,s59)



library(VennDiagram)# for plotting venn diagram
venn.plot<-venn.diagram(x=list(S1=s1,S2=s2,S3=s3),"D025.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.025 [h-1]",
                        col="transparent",
                        fontfamily="arial",
                        rotation.degree=60,
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s4,S2=s5,S3=s6),"D05.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.044 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s7,S2=s8,S3=s9),"D10.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.102 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s10,S2=s11,S3=s12),"D15.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.152 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s13,S2=s14,S3=s15),"D20.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.214 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s16,S2=s17,S3=s18),"D25.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.254 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s19,S2=s20,S3=s21),"D30.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.284 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s22,S2=s23,S3=s24),"D35.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.334 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
venn.plot<-venn.diagram(x=list(S1=s25,S2=s26,S3=s27),"D40.png",
                        imagetype = "png",
                        height = 300,width=300,
                        main="D=0.379 [h-1]",
                        rotation.degree=60,
                        col="transparent",
                        fontfamily="arial",
                        cex=0.2,
                        main.cex = 0.2,
                        cat.cex=0.2,
                        sub.cex=0.1,
                        fill=c("mediumspringgreen","thistle","olivedrab1"))
#Here we select the column iBAQ.H which is the absolute concentration in units of [fmol/6ug total protein]
#从IS表中读取已定量的IS数据，即iBAQ.H列数�?
columnWantedIS<-c("Majority.protein.IDs","iBAQ.H")
dataRef<-subset(dataIS,select=columnWantedIS)

#Get the mutual protein names from the two tables
#从样品表和IS表中取出公共的蛋白数据行
sampleProteins<-dataRelative$Majority.protein.IDs
ISProteins<-dataRef$Majority.protein.IDs
#Get indexes in IS protein id column for each proteins in sample protein list
#首先找出样品蛋白列表中每个蛋白对应的IS蛋白列表中的编号
commonProteinsIdx<-match(ISProteins,sampleProteins)
#Get IS proteins that detected in samples, as their orders in IS table
#然后，将找到的蛋白列表按其在IS蛋白列表中的顺序列出�?
proteinsInBothISAndSample<-ISProteins[!is.na(commonProteinsIdx)]
#Get IS proteins that not deteted in samples, as their orders in IS table
#接着，找出未在样品中检测到的IS蛋白，以其在IS蛋白列表中顺序列�?
proteinsInISButNotInSample<-ISProteins[is.na(commonProteinsIdx)]
#Extract all measured protein 
#然后将所有样品中检测到的IS蛋白的数据保存在extractRowsInSampleData
extractRowsInSampleData<-dataRelative[commonProteinsIdx,]
#Calibrate IS iBAQ.H using bootstrapped regression coefficients
#利用原始文献中iBAQ方法中提到的bootstrapping方法构建IS定量模型模型系数
dataRef$iBAQ.H<-10^(log10(dataRef$iBAQ.H)*1.1252-0.269) #Here 1.1252 and -0.269 were generated using another RMarkDown:BootstrappingRegression.rmd
#根据IS的绝对浓度及样品中的相对值获取计算得到的蛋白浓度
calculatedProteinConcInSampleData<-dataRef$iBAQ.H/extractRowsInSampleData[,2:60]
colnames(calculatedProteinConcInSampleData)<-paste("S",c(1:59),sep="")
#Note:NA or NaN value exist in table
#注意：表格中存在NA和NaN数据
data<-data.frame(proteinIDs=dataRef$Majority.protein.IDs,calculatedProteinConcInSampleData)

#Eliminate ups proteins
#去除ups蛋白�?
upsproteins<-data[grep("*ups",data$proteinIDs),]
data<-data[-grep("*ups",data$proteinIDs),]

#Change unit to [protein molecules/cell], assume yeast cell weight is 13[pg], Avogadro constant is 6.02*10^23, and the change equation is [data]/6*[dataTotalProteinPercent]*13*10^-12*10^6[ug/cell]*10^-15[mol/fmol]*6.02*10^23[molecules/mol]
UnitConvertionConstant=13*1e-6*6.02*1e8/6
# Maybe change to sweep method!
#for (i in 2:60){
#  data[,i]<-data[,i]*dataTotalProteinPercent[,i-1]*UnitConvertionConstant
#}
dataTPN<-sweep(data[,-1],2,colSums(dataTotalProteinPercent*UnitConvertionConstant),"*")
rownames(dataTPN)<-data[,1]

#Output data
#输出数据到excel文件中，分为withoutTPN和withTPN两个表，前者中数据的单位为[fmol/6ugTotalProtein],后者中数据的单位为[protein molecules/cell]
library(openxlsx)
wbData<-createWorkbook()
addWorksheet(wbData,"withoutTPN")
addWorksheet(wbData,"withTPN")
OutputFileName=paste("output/ProteinConc",dateLabel,sep="-")
writeData(wbData,"withoutTPN",data)
writeData(wbData,"withTPN",dataTPN,rowNames = T)
saveWorkbook(wbData,OutputFileName,overwrite = T)
# Substitute by openxlsx functions
#OutputFile=paste("output/ProteinConc",dateLabel,sep="-")
#write.table(data,OutputFile,col.names = TRUE,row.names = FALSE,sep=";")

#Non NA proteins
#将dataframe转换为matrix,并且去除其中的NA数据
dataMatrix<-as.matrix(data[,2:60])
rownames(dataMatrix)<-data[,1]
dataMatrixWithoutWholeRowNAs<-dataMatrix[rowSums(is.na(dataMatrix))!=59,]#Rows that are not totaly NA
dataMatrixNoNARows<-dataMatrix[rowSums(is.na(dataMatrix))==0,]#Rows that are totaly no NA
dataMatrixContainNARows<-dataMatrixWithoutWholeRowNAs[rowSums(is.na(dataMatrixWithoutWholeRowNAs))!=0,]#Rows that contains NA but not total NA
