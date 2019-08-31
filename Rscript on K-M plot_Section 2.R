#=======================================================================================#
# The usages of all packages:                                                           #
#       cgdsr: supply a library to access the data in MSKCC Cancer Genomics Data Server #
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
# load packages
{
  library(cgdsr)
  library(ggpubr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(ggstatsplot)
}

###########################################
# Survival plot analysis of clinical data #
###########################################
# load the data
if(T){
  load("survival_input.Rdata")
}
# view the clinical data
{
  clinicaldata_view<-as.matrix(colnames(myclinicaldata))
}

# read the clinical information
{
  choose_columns=c("DFS_MONTHS","DFS_STATUS","ER_STATUS_BY_IHC","PR_STATUS_BY_IHC","IHC_HER2")
  choose_clinicaldata=myclinicaldata[,choose_columns]
}

# choose the specified clinical data
{
  dat=choose_clinicaldata[choose_clinicaldata$ER_STATUS_BY_IHC=="Negative",]
  dat1=dat[dat$PR_STATUS_BY_IHC=="Negative",]
  dat2=dat1[dat1$IHC_HER2=="Negative",]
  dat3=dat2[dat2$DFS_MONTHS>0,]
  dat3=dat3[!is.na(dat3$DFS_STATUS),]
  dat4=cbind(dat3[,c("DFS_STATUS","DFS_MONTHS")],expr[rownames(dat3),])
  colnames(dat4)[3]<-"SDC4"
}
write.csv(dat4,"SDC4_TCGA_expr.csv")
# expressed genes plot
{
  p<-ggboxplot(dat4,x="DFS_STATUS",y="SDC4",color="DFS_STATUS",palette="jco",add="jitter")
  p+stat_compare_means(method="t.test")
  dat4$SDC4_group=ifelse(dat4$SDC4>median(dat4$SDC4),'high','low')
  ggbetweenstats(data=dat4,x=SDC4_group,y=SDC4)
}

# survival plot
attach(dat4)
{
  table(SDC4_group)
  my.surv<-Surv(DFS_MONTHS,DFS_STATUS=='Recurred/Progressed')
  kmfit3<-survfit(my.surv~SDC4_group,data=dat4)
  plot(kmfit3,col=c("red","blue"))
  ggsurvplot(kmfit3,palette=c("#E7B800","#2E9FDF"),
             conf.int=TRUE,pval=TRUE,xlab="Time / Month",
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