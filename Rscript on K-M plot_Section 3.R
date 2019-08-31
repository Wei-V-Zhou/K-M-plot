#=======================================================================================#
# The usages of all packages:                                                           #
#      ggpubr: load a required package:magrittr to offer operators to promote semantics #
#     ggplot2: create an elegant data visualization using the grammar of graphics       #
#    survival: offer the K-M survival curve and operate some survival analysis          #
#   survminer: preform the cox survival analysis  of some parameters                    #
# ggstatsplot: show more statistical details based on ggplot2                           #
#=======================================================================================#

##################
# Load libraries #
##################
# clear objectives and garbage collection
rm(list=ls())
gc()
options(stringsAsFactors=FALSE)
# load packages
{
  library(plyr)
  library(limma)
  library(ggpubr)
  library(stringr)
  library(ggplot2)
  library(Biobase)
  library(survival)
  library(reshape2)
  library(GEOquery)
  library(survminer)
  library(ggfortify)
  library(ggstatsplot)
}

############################################
# Load expression data and gene annotation #
############################################
# load data
if(!file.exists("GSE2603_gSet.Rdata")){
  gSet<-getGEO("GSE2603",destdir=".",GSEMatrix=TRUE,
               AnnotGPL=FALSE,getGPL=FALSE)
  if (length(gSet)>1) idx<-grep("GPL96",attr(gSet,"names")) else idx<-1
  gSet<-gSet[[idx]]
  save(gSet,file="GSE2603_gSet.Rdata")
}
load("GSE2603_gSet.Rdata")
# log2 transform
ex<-exprs(gSet)
qx<-as.numeric(quantile(ex,c(0.,0.25,0.5,0.75,0.99,1.0),na.rm=T))
LogC<-(qx[5]>100)||
  (qx[6]-qx[1]>50 && qx[2]>0)||
  (qx[2]>0 && qx[2]<1 && qx[4]>1 && qx[4]<2)
if(LogC){ex[which(ex<=0)]<-NaN
exprs(gSet)<-log2(ex)}
# load the expression data
exprset<-exprs(gSet)
exprSet<-exprset[,23:ncol(exprset)]
dim(exprSet)
# load the sample name
samples<-sampleNames(gSet)
samples<-samples[23:ncol(exprset)]
# load the annotation
pdata<-pData(gSet)
pdata<-as.data.frame(pdata)
pdata<-pdata[23:ncol(exprset),]

# load the GPL data
if(T){
  gpl <- getGEO("GPL96", destdir = ".")
  probe2gene<-Table(gpl)[,c(1,11)]
}

#############################
# Transform probes to genes #
#############################
# filter the probe with null genes
ids=probe2gene[probe2gene[,2]!='',]
# filter the probe with many genes
a<-strsplit(as.character(ids[,2])," /// ")
tmp<-mapply(cbind,ids[,1],a)
df<-ldply(tmp,data.frame)
probe2gene=df[,2:3]
save(probe2gene,file="probe2gene_K-M.Rdata")

##################
# Data filtering #
##################
# remove the probe without gene annotation
if(T){
  exprSet=exprSet[rownames(exprSet) %in% probe2gene[,1],]
  probe2gene=probe2gene[match(rownames(exprSet),probe2gene[,1]),]
}
dim(exprSet)
dim(probe2gene)
tail(sort(table(probe2gene[,2])),n=12L)
# # select the JAG1 expression
# {
#   x1=exprSet[rownames(exprSet)=="209097_s_at",]
#   x2=exprSet[rownames(exprSet)=="209098_s_at",]
#   x3=exprSet[rownames(exprSet)=="209099_x_at",]
#   x4=exprSet[rownames(exprSet)=="216268_s_at",]
# }
# JAG1=cbind(x1,x2,x3,x4)
# write.csv(JAG1,"JAG1_expression.csv")

# get the maximum expression of the same gene
{
  MAX=by(exprSet,probe2gene[,2],
         function(x) rownames(x)[which.max(rowMeans(x))])
  MAX=as.character(MAX)
  exprSet=exprSet[rownames(exprSet) %in% MAX,]
  rownames(exprSet)=probe2gene[match(rownames(exprSet),probe2gene[,1]),2]
}
dim(exprSet)
exprSet[1:5,1:5]
save(exprSet,file="exprset_K-M.Rdata")

###########################################
# Survival plot analysis of clinical data #
###########################################
# clear objectives and garbage collection
rm(list=ls())
gc()
options(stringsAsFactors=FALSE)
# load the data
if(T){
  load("exprset_K-M.Rdata")
}
exprset<-as.matrix(exprSet[rownames(exprSet)=="SDC4",])
# exprset1<-as.matrix(exprSet[rownames(exprSet)=="HTRA1",])
# exprset2<-as.matrix(exprSet[rownames(exprSet)=="CYR61",])
# EXP<-(exprset+exprset1+exprset2)/3
myclinicaldata<-read.csv("myclinicaldata.csv",header=T)
rownames(myclinicaldata)<-myclinicaldata[,1]
# view the clinical data
{
  clinicaldata_view<-as.matrix(colnames(myclinicaldata))
}

# read the clinical information
{
  choose_columns=c("Path.ER.status","Path.PR.status","Her2.status","van.t.Veer.Signature","Sample.Name",
                   "Sample.Id","MFS..YR.","Met.Event","LM.Event","LMFS..YR.","BM.Event","BMFS..YR.")
  choose_clinicaldata=myclinicaldata[,choose_columns]
}

# choose the specified clinical data
{
  dat1=choose_clinicaldata[!is.na(choose_clinicaldata$BMFS..YR.),]
  dat4=cbind(dat1[,c("BM.Event","BMFS..YR.")],EXP[rownames(dat1),])
  colnames(dat4)[3]<-"EXP_NO_LOG"
}
# write.csv(dat4,"EXP_NO_LOG.csv")
dat4$MEDIAN_LOG2_group=ifelse(dat4$MEDIAN_LOG2>median(dat4$MEDIAN_LOG2),'high','low')

# expressed genes plot
{
  p<-ggboxplot(dat4,x="van.t.Veer.Signature",y="HTRA1",color="van.t.Veer.Signature",palette="jco",add="jitter")
  p+stat_compare_means(method="t.test")
  dat4$JAG1_group=ifelse(dat4$JAG1>median(dat4$JAG1),'high','low')
  ggbetweenstats(data=dat4,x=HTRA1_group,y=HTRA1)
}

# survival plot
attach(dat4)
{
  table(MEDIAN_LOG2_group)
  my.surv<-Surv(BMFS..YR.,BM.Event=='1')
  kmfit<-survfit(my.surv~MEDIAN_LOG2_group,data=dat4)
  plot(kmfit,col=c("red","blue"),xlab="Years",ylab="Bone met-free")
  ggsurvplot(kmfit,palette=c("#E7B800","#2E9FDF"),
             conf.int=TRUE,pval=TRUE,xlab="Years",
             ggtheme=theme_light(),risk.table=TRUE,ncensor.plot=TRUE)
}
detach(dat4)

# cox survival plot
dat4$DFS_STATUS=as.character(dat4$DFS_STATUS)
str(dat4,no.list=T,vec.len=2)
attach(dat4)
{
  my.surv<-Surv(DFS_MONTHS,DFS_STATUS=='Recurred/Progressed')
  m=coxph(my.surv~HTRA1+HTRA1_group,data=dat4)
  ggsurvplot(survfit(m,data=dat4), color="#2E9FDF",
             ggtheme=theme_minimal())
}
detach(dat4)

#============================#
#       Musician: Resonance  #
#           Date: 2019/08/17 #
# Revised author: Resonance  #
#           Time: 2019/08/26 #
#============================#