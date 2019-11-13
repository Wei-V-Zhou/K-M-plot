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
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# load packages
{
  library(cgdsr)
  library(ggpubr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(ggstatsplot)
}

########################
# Connect the database #
########################
# create a CGDS connection objext
mycgds <- CGDS("http://www.cbioportal.org/")
# set verbose options to debug and troubleshoot issues
setVerbose(mycgds,TRUE)
# guarantee the unauthorized accessment
mysecurecgds <- CGDS("http://cbioportal.mskcc.org/",
                     token="fd0522cb-7972-40d0-9d83-cb4c14e8a337")

#########################
# Get the required data #
#########################
# get list of cancer studies at server for view to choose
if(F){
  cancerstudy <- getCancerStudies(mycgds)
}
# get available case lists for a given cancer study
{
  mycancerstudy <- getCancerStudies(mycgds)[46, 1]
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[5, 1]
}

# get available genetic profiles
{
  mygeneticprofile <- getGeneticProfiles(mycgds, mycancerstudy)[2, 1]
}

# get data slices for a specified list of genes, genetic profiles and case list
{
  choose_genes <- c("JAG2")
  expr = getProfileData(mycgds, choose_genes, mygeneticprofile, mycaselist)
}

# get clinical data for the case list
{
  myclinicaldata <- getClinicalData(mycgds, mycaselist)
}

# save the data
save(expr, myclinicaldata, file="survival_inputdata.Rdata")

##################
# Load libraries #
##################
# clear objectives and garbage collection
rm(list = ls())
gc()
options(stringsAsFactors = F)
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
  load("survival_inputdata.Rdata")
}
# view the clinical data
{
  clinicaldata_view <- as.matrix(colnames(myclinicaldata))
}

# read the clinical information
{
  choose_columns = c("DFS_MONTHS", "DFS_STATUS", "SUBTYPE")
  choose_clinicaldata = myclinicaldata[ , choose_columns]
}

# choose the specified clinical data
{
  dat2 = choose_clinicaldata[choose_clinicaldata$SUBTYPE == "BRCA_Basal", ]
  dat3 = dat2[dat2$DFS_MONTHS > 0 & dat2$DFS_MONTHS < 100, ]
  dat3 = dat3[!is.na(dat3$DFS_STATUS), ]
  dat4 = cbind(dat3[ , c("DFS_STATUS", "DFS_MONTHS")], expr[rownames(dat3), ])
  colnames(dat4)[3] <- "JAG2"
}
write.csv(dat4, "JAG2_TCGA_expr_luminalB.csv")
# expressed genes plot
{
  p <- ggboxplot(dat4, x="DFS_STATUS", y="JAG2", color="DFS_STATUS", palette="jco", add="jitter")
  p + stat_compare_means(method = "t.test")
  dat4$JAG2_group = ifelse(dat4$JAG2 > median(dat4$JAG2), 'high', 'low')
  # dat4$JAG2_group = ifelse(dat4$JAG2 > quantile(dat4$JAG2)[4], 'high', 'low')
  ggbetweenstats(data = dat4, x = JAG2_group, y = JAG2)
}

# survival plot
attach(dat4)
{
  table(JAG2_group)
  my.surv <- Surv(DFS_MONTHS, DFS_STATUS == 'Recurred/Progressed')
  kmfit3 <- survfit(my.surv~JAG2_group, data = dat4)
  plot(kmfit3, col = c("red", "blue"))
  ggsurvplot(kmfit3, palette=c("#E7B800","#2E9FDF"),
             conf.int = TRUE, pval = TRUE, xlab = "Time / Month",
             ggtheme = theme_light(), risk.table = TRUE, ncensor.plot = TRUE)
}
detach(dat4)

# cox survival plot
dat4$DFS_STATUS=as.character(dat4$DFS_STATUS)
str(dat4,no.list=T,vec.len=2)
attach(dat4)
{
  my.surv<-Surv(DFS_MONTHS, DFS_STATUS=='Recurred/Progressed')
  m=coxph(my.surv~JAG2 + JAG2_group, data = dat4)
  ggsurvplot(survfit(m, data = dat4), color="#2E9FDF",
             ggtheme = theme_minimal())
}
detach(dat4)

#============================#
#       Musician: Resonance  #
#           Date: 2019/09/23 #
# Revised author: Resonance  #
#           Time: 2019/11/11 #
#============================#