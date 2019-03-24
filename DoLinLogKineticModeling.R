## PURPOSE ##
# For do linear regression of lin-log kinetic modelling of glycolysis analysis
# Set working directory
PCPrefix="~/Work/Sweden/R"
ScriptLoc=paste(PCPrefix,"RawdataProcess/Proteome/IntegrateMetabolome",sep="/")
setwd(ScriptLoc)

#data<-read.csv("GlycolysisReactions-r1.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
#data<-read.csv("GlycolysisReactions-r2.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)#Add AMP, P, CIT to the regression of PFK
#data<-data[-c(6),]#eliminate the data of D0.3vs0.1
data<-read.csv("GlycolysisReactions-r3.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)#Corrected the calculation error when normalizing with biomass
dataTCA<-read.csv("TCA.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
dataPPP<-read.csv("PPP.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)

#Metabolites
atp<-log(data$ATP)
adp<-log(data$ADP)
amp<-log(data$AMP)
p<-log(data$P)
glc<-log(data$GLC)
g6p<-log(data$G6P)
f6p<-log(data$F6P)
fbp<-log(data$FBP)
dhap<-log(data$DHAP)
g3p<-log(data$G3P)
p13g<-log(data$P13G)
nad<-log(data$NAD)
nadh<-log(data$NADH)
x3pg<-log(data$X3PG)
x2pg<-log(data$X2PG)
pep<-log(data$PEP)
pyr<-log(data$PYR)
cit<-log(data$CIT)

oxo<-log(dataTCA$OXO) #oxoglutarate
fum<-log(dataTCA$FUM) #fumarate
isocit<-log(dataTCA$ISOCIT)
mal<-log(dataTCA$MAL)#malate
suc<-log(dataTCA$SUC)#succinate

x6pg<-log(dataPPP$X6pg)
r5p<-log(dataPPP$r5p)
s7p<-log(dataPPP$s7p)

#EMP pathway
rHXK_specific_1<-data$rHXK/data$GLK1
r<-c(rHXK_specific_1)
fitHXK1<-lm((r-1)~0+atp+glc+adp+g6p)
summary(fitHXK1)

rHXK_specific_2<-data$rHXK/data$HXK1
r<-c(rHXK_specific_2)
fitHXK2<-lm((r-1)~0+atp+glc+adp+g6p)
summary(fitHXK2)

rHXK_specific_3<-data$rHXK/data$HXK2
r<-c(rHXK_specific_3)
fitHXK3<-lm((r-1)~0+atp+glc+adp+g6p)
summary(fitHXK3)

rHXK_specific_sum<-data$rHXK/data$sumR_HXK
r<-c(rHXK_specific_sum)
fitHXKsum<-lm((r-1)~0+atp+glc+adp+g6p)
summary(fitHXKsum)

rPGI_specific<-data$rPGI/data$PGI1
r<-c(rPGI_specific)
fitPGI<-lm((r-1)~0+g6p+f6p)
summary(fitPGI)

rPFK_specific_1<-data$rPFK/data$PFK1
r<-c(rPFK_specific_1)
fitPFK1<-lm((r-1)~0+f6p+fbp+atp+adp+p+amp+cit)
summary(fitPFK1)

#fit3.1<-lm((rOverEnzyme3.1-1)~0+f6p+fbp+adp+atp)
#summary(fit3.1)

rPFK_specific_2<-data$rPFK/data$PFK2
r<-c(rPFK_spcific_2)
fitPFK2<-lm((r-1)~0+f6p+fbp+adp+atp+cit+p+amp)
summary(fitPFK2)

rPFK_specific_sum<-data$rPFK/data$sumR_PFK
r<-c(rPFK_specific_sum)
fitPFKsum<-lm((r-1)~0+f6p+fbp+adp+atp+cit+p+amp)
summary(fitPFKsum)


rFBA_specific<-data$rFBA/data$FBA1
r<-c(rFBA_specific)
fitFBA<-lm((r-1)~0+fbp+dhap)
summary(fitFBA)

rTPI_specific<-data$rTPI/data$TPI1
r<-c(rTPI_specific)
fitTPI<-lm((r-1)~0+dhap)
summary(fitTPI)

rTDH_specific_1<-data$rTDH/data$TDH1
r<-c(rTDH_specific_1)
fitTDH1<-lm((r-1)~0+g3p+nad+p+p13g+nadh)
summary(fitTDH1)

rTDH_specific_2<-data$rTDH/data$TDH2
r<-c(rTDH_specific_2)
fitTDH2<-lm((r-1)~0+g3p+nad+p+p13g+nadh)
summary(fitTDH2)

rTDH_specific_3<-data$rTDH/data$TDH3
r<-c(rTDH_specific_3)
fitTDH3<-lm((r-1)~0+g3p+nad+p+p13g+nadh)
summary(fitTDH3)

rTDH_specific_sum<-data$rTDH/data$sumR_TDH
r<-c(rTDH_specific_sum)
fitTDHsum<-lm((r-1)~0+g3p+nad+p+p13g+nadh)
summary(fitTDHsum)

rPGK_specific<-data$rPGK/data$PGK1
r<-c(rPGK_specific)
fitPGK<-lm((r-1)~0+p13g+x3pg+adp+atp)
summary(fitPGK)

rGPM_specific<-data$rGPM/data$GPM1
r<-c(rGPM_specific)
fitGPM<-lm((r-1)~0+x3pg+x2pg)
summary(fitGPM)

rENO_specific_1<-data$rENO/data$ENO1
r<-c(rENO_specific_1)
fitENO1<-lm((r-1)~0+x2pg+pep)
summary(fitENO1)

rENO_specific_2<-data$rENO/data$ENO2
fitENO2<-lm((r-1)~0+x2pg+pep)
summary(fitENO2)

rENO_specific_sum<-data$rENO/data$sumR_ENO
r<-c(rENO_specific_sum)
fitENOsum<-lm((r-1)~0+x2pg+pep)
summary(fitENOsum)

rPYK<-data$rPYK/data$CDC19
r<-c(rPYK)
x1<-log(data$PEP)
x2<-log(data$PYR)
x3<-log(data$ADP)
x4<-log(data$ATP)
fitPYK<-lm((r-1)~0+pep+pyr+adp+atp)
summary(fitPYK)

#function help to extract p-value of regression model
extractPValue<-function(modelobj){
  if(class(modelobj)!="lm") stop("Not an object of class 'lm' ")
  f<-summary(modelobj)$fstatistic
  p<-pf(f[1],f[2],f[3],lower.tail = F)
  attributes(p)<-NULL
  return(p)
}

adj.r.squared<-c(summary(fitHXK1)$adj.r.squared,summary(fitHXK2)$adj.r.squared,summary(fitHXK3)$adj.r.squared,summary(fitHXKsum)$adj.r.squared,summary(fitPGI)$adj.r.squared,
                 summary(fitPFK1)$adj.r.squared,summary(fitPFK2)$adj.r.squared,summary(fitPFKsum)$adj.r.squared,summary(fitFBA)$adj.r.squared,summary(fitTPI)$adj.r.squared,
                 summary(fitTDH1)$adj.r.squared,summary(fitTDH2)$adj.r.squared,summary(fitTDH3)$adj.r.squared,summary(fitTDHsum)$adj.r.squared,summary(fitPGK)$adj.r.squared,
                 summary(fitGPM)$adj.r.squared,summary(fitENO1)$adj.r.squared,summary(fitENO2)$adj.r.squared,summary(fitENOsum)$adj.r.squared,summary(fitPYK)$adj.r.squared)
p.values<-c(extractPValue(fitHXK1),extractPValue(fitHXK2),extractPValue(fitHXK3),extractPValue(fitHXKsum),extractPValue(fitPGI),
            extractPValue(fitPFK1),extractPValue(fitPFK2),extractPValue(fitPFKsum),extractPValue(fitFBA),extractPValue(fitTPI),
            extractPValue(fitTDH1),extractPValue(fitTDH2),extractPValue(fitTDH3),extractPValue(fitTDHsum),extractPValue(fitPGK),
            extractPValue(fitGPM),extractPValue(fitENO1),extractPValue(fitENO2),extractPValue(fitENOsum),extractPValue(fitPYK))
#model.info<-c(summary(fitHXK1)$coefficients,summary(fitHXK2)$coefficients,summary(fitHXK3)$coefficients,summary(fitHXKsum)$coefficients,summary(fitPGI)$coefficients,
#              summary(fitPFK1)$coefficients,summary(fitPFK2)$coefficients,summary(fitPFKsum)$coefficients,summary(fitFBA)$coefficients,summary(fitTPI)$coefficients,
#              summary(fitTDH1)$coefficients,summary(fitTDH2)$coefficients,summary(fitTDH3)$coefficients,summary(fitTDHsum)$coefficients,summary(fitPGK)$coefficients,
#              summary(fitGPM)$coefficients,summary(fitENO1)$coefficients,summary(fitENO2)$coefficients,summary(fitENOsum)$coefficients,summary(fitPYK)$coefficients)
model_fitness<-data.frame(adj.r.squared,p.values)
rownames(model_fitness)<-c("HXK1","HXK2","HXK3","HXKsum","PGI",
                           "PFK1","PFK2","PFKsum","FBA","TPI",
                           "TDH1","TDH2","TDH3","TDHsum","PGK",
                           "GPM","ENO1","ENO2","ENOsum","PYK")
write.csv(model_fitness,file="EMP-model.fitness.csv")

rownames(model_fitness)[model_fitness$adj.r.squared>0.9]

model_fitness[model_fitness$adj.r.squared>0.95&model_fitness$p.values<0.05,]

result<-list(hxk1=summary(fitHXK1),hxk2=summary(fitHXK2),hxk3=summary(fitHXK3),hxksum=summary(fitHXKsum),pgi=summary(fitPGI),
             pfk1=summary(fitPFK1),pfk2=summary(fitPFK2),pfksum=summary(fitPFKsum),fba=summary(fitFBA),tpi=summary(fitTPI),
             tdh1=summary(fitTDH1),tdh2=summary(fitTDH2),tdh3=summary(fitTDH3),tdhsum=summary(fitTDHsum),pgk=summary(fitPGK),
             gpm=summary(fitGPM),eno1=summary(fitENO1),eno2=summary(fitENO2),enosum=summary(fitENOsum),pyk=summary(fitPYK))

# TCA pathway
rACO_specific<-dataTCA$rACO/dataTCA$pACO1
fitACO<-lm((rACO_specific-1)~0+cit+isocit)
summary(fitACO)

rIDH_specific<-dataTCA$rIDH/dataTCA$PIDH1_2
fitIDH<-lm((rIDH_specific-1)~0+isocit+oxo+nad+nadh)
summary(fitIDH)

rSDH_specific<-dataTCA$rSDH/dataTCA$pSDH1_4
fitSDH<-lm((rSDH_specific-1)~0+suc+fum)
summary(fitSDH)

rFUM_specific<-dataTCA$rFUM/dataTCA$pFUM1
fitFUM<-lm((rFUM_specific-1)~0+mal+fum)
summary(fitFUM)

adj.r.squared<-c(summary(fitACO)$adj.r.squared,summary(fitIDH)$adj.r.squared,summary(fitSDH)$adj.r.squared,summary(fitFUM)$adj.r.squared)
p.values<-c(extractPValue(fitACO),extractPValue(fitIDH),extractPValue(fitSDH),extractPValue(fitFUM))

model_fitness_TCA<-data.frame(adj.r.squared,p.values)
rownames(model_fitness_TCA)<-c("ACO1","IDH1,2","SDH1-4","SUM1")
write.csv(model_fitness_TCA,file="TCA_model.fitness.csv")

# Pentos phosphate pathway
rZWF_specific<-dataPPP$rZWF/dataPPP$pZWF1
fitZWF<-lm((rZWF_specific-1)~0+g6p+x6pg)
summary(fitZWF)

rSOL_specific_2<-dataPPP$rSOL/dataPPP$pSOL2
fitSOL2<-lm((rSOL_specific_2-1)~0+g6p+x6pg)
summary(fitSOL2)

rSOL_specific_1<-dataPPP$rSOL/dataPPP$pSOL1
fitSOL1<-lm((rSOL_specific_1-1)~0+g6p+x6pg)
summary(fitSOL1)

rRPE_specific<-dataPPP$rRPE/dataPPP$pRPE1
fitRPE<-lm((rRPE_specific-1)~0+r5p)
summary(fitRPE)

rTKL_specific<-dataPPP$rTKL/dataPPP$pTKL1
fitTKL<-lm((rTKL_specific-1)~0+r5p+s7p+g3p)
summary(fitTKL)

adj.r.squared<-c(summary(fitZWF)$adj.r.squared,summary(fitSOL2)$adj.r.squared,summary(fitSOL1)$adj.r.squared,summary(fitRPE)$adj.r.squared,summary(fitTKL)$adj.r.squared)
p.values<-c(extractPValue(fitZWF),extractPValue(fitSOL2),extractPValue(fitSOL1),extractPValue(fitRPE),extractPValue(fitTKL))

model_fitness_PPP<-data.frame(adj.r.squared,p.values)
rownames(model_fitness_PPP)<-c("ZWF","SOL2","SOL1","RPE","TKL")
write.csv(model_fitness_PPP,file="PPP_model.fitness.csv")
