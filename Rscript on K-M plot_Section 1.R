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

########################
# Connect the database #
########################
# create a CGDS connection objext
mycgds<-CGDS("http://www.cbioportal.org/")
# set verbose options to debug and troubleshoot issues
setVerbose(mycgds,TRUE)
# guarantee the unauthorized accessment
mysecurecgds<-CGDS("http://cbioportal.mskcc.org/",
                   token="fd0522cb-7972-40d0-9d83-cb4c14e8a337")

#########################
# Get the required data #
#########################
# get list of cancer studies at server for view to choose
if(F){
  cancerstudy<-getCancerStudies(mycgds)
}
# get available case lists for a given cancer study
{
  mycancerstudy<-getCancerStudies(mycgds)[46,1]
  mycaselist<-getCaseLists(mycgds,mycancerstudy)[8,1]
}

# get available genetic profiles
{
  mygeneticprofile<-getGeneticProfiles(mycgds,mycancerstudy)[8,1]
}

# get data slices for a specified list of genes, genetic profiles and case list
{
  choose_genes<-c("SDC4")
  expr=getProfileData(mycgds,choose_genes,mygeneticprofile,mycaselist)
}

# get clinical data for the case list
{
  myclinicaldata<-getClinicalData(mycgds,mycaselist)
}

# save the data
save(expr,myclinicaldata,file="survival_input.Rdata")

#============================#
#       Musician: Resonance  #
#           Date: 2019/08/17 #
# Revised author: Resonance  #
#           Time: 2019/08/26 #
#============================#