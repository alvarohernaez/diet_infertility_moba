rm(list=ls())

library(chron)
library(colorspace)
library(mime)
library(munsell)
library(labeling)
library(rlang)
library(stringi)
library(evaluate)
library(highr)
library(markdown)
library(yaml)
library(backports)
library(jsonlite)
library(digest)
library(plyr)
library(reshape2)
library(scales)
library(tibble)
library(lazyeval)
library(RColorBrewer)
library(stringr)
library(knitr)
library(magrittr)
library(checkmate)
library(htmlwidgets)
library(viridisLite)
library(Rcpp)
library(Formula)
library(ggplot2)
library(latticeExtra)
library(acepack)
library(gtable)
library(data.table)
library(htmlTable)
library(viridis)
library(htmltools)
library(base64enc)
library(minqa)
library(RcppEigen)
library(lme4)
library(SparseM)
library(MatrixModels)
library(pbkrtest)
library(quantreg)
library(car)
library(Hmisc)
library(survival)
library(foreign)
library(bitops)
library(caTools)
library(gplots)
library(ROCR)
library(mice)
library(writexl)
library(HardyWeinberg)
library(officer)
library(uuid)
library(compareGroups)
library(nlme)
library(vcd)
library(boot)
library(tibble)
library(haven)
library(icenReg)
library(MASS)
library(sandwich)   
library(lmtest)
library(gam)
library(smoothHR)
library(metafor)
library(DBI)
library(mitools)
library(RcppArmadillo)
library(miceadds)
library(dplyr)
library(estimatr)
library(lubridate)
library(snakecase)
library(janitor)
library(fmsb)

# Not possible to install "RadialMR": copy installed folder from local library outside TSD to library in TSD
# C:\Users\p1775-alvarohc\Documents\R\win-library\4.1
library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)
library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR) 
library(ieugwasr)
library(TwoSampleMR)

RutinesLocals<- "N:/durable/Project/ALHE_Infertility/R_packages/Routines"
source(file.path(RutinesLocals,"carrega.llibreria.r"))
source(file.path(RutinesLocals,"merge2.r"))
source(file.path(RutinesLocals,"fix2.r"))
source(file.path(RutinesLocals,"table2.r"))
source(file.path(RutinesLocals,"subset2.r"))
source(file.path(RutinesLocals,"format2.r"))
source(file.path(RutinesLocals,"order2.r"))
source(file.path(RutinesLocals,"intervals.r"))


### GUAPAS ###
##############

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.00001,signif(x,1),
                   ifelse(abs(x)<0.0001,signif(x,1),
                          ifelse(abs(x)<0.001,signif(x,1),
                                 ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                        ifelse(abs(x)<1,sprintf("%.2f",round(x,2)),
                                               ifelse(abs(x)<10,sprintf("%.2f",round(x,2)),
                                                      ifelse(abs(x)<100,sprintf("%.1f",round(x,1)),
                                                             ifelse(abs(x)>=100,round(x,0),round(x,0)))))))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

ic_guapa2<-function(x,y,z)
{
  ic<-paste(x," (",y," to ",z,")",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",
                      ifelse(abs(x)<0.01,sprintf("%.3f",round(x,3)),
                             ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                    ifelse(abs(x)<1,sprintf("%.3f",round(x,3)),guapa(x))))))
  return(pval)
}

pval_guapa2<-function(x)
{
  pval<-ifelse(x<0.00001," < 0.00001",
               ifelse(x<0.001," < 0.001",
                      ifelse(abs(x)<0.01,sprintf("%.3f",round(x,3)),
                             ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                    ifelse(abs(x)<1,sprintf("%.3f",round(x,3)),guapa(x))))))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

beta_se_ic_guapa2 <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa2(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa3 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),3)
  ic95a<-round(exp(x-(z*y)),3)
  ic95b<-round(exp(x+(z*y)),3)
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

# Non-linear MR function
# source("N:/data/durable/Projects/Magnus_MR_BMI/R/Old/nlme_summ_aes.R")

dir.create("N:/durable/Project/ALHE_diet_infertility/Outputs")
dir.create("N:/durable/Project/ALHE_diet_infertility/Outputs/Descriptive")
dir.create("N:/durable/Project/ALHE_diet_infertility/Outputs/Results")

memory.limit(35000)
setwd("N:/durable/Project/ALHE_diet_infertility")


### HORN'S PARALLELL ANALYSIS OF PRINCIPAL COMPONENTS ###
#########################################################

library(paran)
load("./Data/MoBa_diet.RData")

datx<-dat[!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$diet_mom==1,]
dat2<-datx[,c("fruit_mom","veg_mom","nutleg_mom","wholeg_mom","refg_mom","dairy_tot_mom","cof_tea_mom",
              "sbev_mom","snack_proxy_mom","eggs_mom","fish_sf_mom","rpmeat_mom","wmeat_mom")]
pc<-princomp(na.omit(dat2),cor=TRUE)
summary(pc)
paran<-paran(na.omit(dat2))
paran$Retained

# 4 components retained for mother's diet #

datx<-dat[!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$diet_dad==1,]
dat2<-datx[,c("fruit_dad","veg_dad","nutleg_dad","wholeg_dad","refg_dad","dairy_tot_dad","cof_tea_dad",
              "sbev_dad","snack_proxy_dad","eggs_dad","fish_dad","rpmeat_dad","wmeat_dad")]
pc<-princomp(na.omit(dat2),cor=TRUE)
summary(pc)
paran<-paran(na.omit(dat2))
paran$Retained

# 4 components retained for father's diet #


### DESCRIPTIVE ###
###################

### POPULATION DESCRIPTION ###
##############################

load("./Data/MoBa_diet.RData")

datx<-dat[!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$diet_mom==1 & dat$preg_plan==1,]
xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Mothers-All")
write.table(tab1,file="./Outputs/Descriptive/descr_diet_mothers.csv",sep=";",col.names=NA)


datx<-dat[!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$diet_dad==1 & dat$preg_plan==1,]
xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Fathers-All")
write.table(tab1,file="./Outputs/Descriptive/descr_diet_fathers.csv",sep=";",col.names=NA)


datx<-dat[!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$metabolomics==1 & dat$preg_plan==1,]
xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Mothers-All")
write.table(tab1,file="./Outputs/Descriptive/descr_metab_mothers.csv",sep=";",col.names=NA)


datx<-dat[!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$metabolomics==1 & dat$preg_plan==1,]
xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Fathers-All")
write.table(tab1,file="./Outputs/Descriptive/descr_metab_fathers.csv",sep=";",col.names=NA)


### SELECTION BIAS - DIET vs NO DIET ###

datx<-dat[!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom),]
datx$diet_mom<-with(datx,ifelse(preg_plan==1,diet_mom,
                                ifelse(preg_plan==0,0,diet_mom)))
datx$metabolomics<-with(datx,ifelse(preg_plan==1,metabolomics,
                                    ifelse(preg_plan==0,0,metabolomics)))
xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom","diet_mom","metabolomics")]

all<-NULL
all<-createTable(compareGroups(diet_mom~. -metabolomics,
                               xxx, method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Non-included","Included","P-value","N")
write.table(tab1,file="./Outputs/Descriptive/selectionbias_diet_mothers.csv",sep=";",col.names=NA)


all<-NULL
all<-createTable(compareGroups(metabolomics~. -diet_mom,
                               xxx, method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Non-included","Included","P-value","N")
write.table(tab1,file="./Outputs/Descriptive/selectionbias_metab_mothers.csv",sep=";",col.names=NA)


datx<-dat[!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad),]
datx$diet_dad<-with(datx,ifelse(preg_plan==1,diet_dad,
                                ifelse(preg_plan==0,0,diet_dad)))
datx$metabolomics<-with(datx,ifelse(preg_plan==1,metabolomics,
                                    ifelse(preg_plan==0,0,metabolomics)))
xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad","diet_dad","metabolomics")]

all<-NULL
all<-createTable(compareGroups(diet_dad~. -metabolomics,
                               xxx, method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Non-included","Included","P-value","N")
write.table(tab1,file="./Outputs/Descriptive/selectionbias_diet_fathers.csv",sep=";",col.names=NA)


all<-NULL
all<-createTable(compareGroups(metabolomics~. -diet_dad,
                               xxx, method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Non-included","Included","P-value","N")
write.table(tab1,file="./Outputs/Descriptive/selectionbias_metab_fathers.csv",sep=";",col.names=NA)


### LOGISTIC REGRESSION / MENDELIAN RANDOMIZATION: LINEAR ASSOCIATIONS ###
##########################################################################

load("./Data/MoBa_diet.RData")

vars00<-c("dash_mom","dash_dad",
          "fruit_mom","fruit_dad","veg_mom","veg_dad","fv_lowpest_mom","fv_highpest_mom","nutleg_mom","nutleg_dad",
          "wholeg_mom","wholeg_dad","refg_mom","refg_dad","cof_tea_mom","cof_tea_dad","sbev_mom","sbev_dad",
          "dairy_tot_mom","dairy_tot_dad","dairy_lf_mom","dairy_lf_dad","fish_mom","fish_dad","rpmeat_mom","rpmeat_dad","wmeat_mom","wmeat_dad",
          "snack_proxy_mom","snack_proxy_dad","eggs_mom","eggs_dad")
vars01<-c("dash_mom_z","dash_dad_z",
          "fruit_mom_z","fruit_dad_z","veg_mom_z","veg_dad_z","fv_lowpest_mom_z","fv_highpest_mom_z","nutleg_mom_z","nutleg_dad_z",
          "wholeg_mom_z","wholeg_dad_z","refg_mom_z","refg_dad_z","cof_tea_mom_z","cof_tea_dad_z","sbev_mom_z","sbev_dad_z",
          "dairy_tot_mom_z","dairy_tot_dad_z","dairy_lf_mom_z","dairy_lf_dad_z","fish_mom_z","fish_dad_z","rpmeat_mom_z","rpmeat_dad_z","wmeat_mom_z","wmeat_dad_z",
          "snack_proxy_mom_z","snack_proxy_dad_z","eggs_mom_z","eggs_dad_z")
vars02<-c("agedelivery_mom","agedelivery_dad",
          "agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_mom","agedelivery_mom","agedelivery_dad",
          "agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad",
          "agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad",
          "agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad")
vars03<-c("eduyears_mom","eduyears_dad",
          "eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_mom","eduyears_mom","eduyears_dad",
          "eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad",
          "eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad",
          "eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad")
vars04<-c("bmi_mom","bmi_dad",
          "bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_mom","bmi_mom","bmi_dad",
          "bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_dad",
          "bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_dad",
          "bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_mom")
vars05<-c("smkinit_mom","smkinit_dad",
          "smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_mom","smkinit_mom","smkinit_dad",
          "smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad",
          "smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad",
          "smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad")
vars06<-c("m_id_2535","f_id_2535",
          "m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","m_id_2535","m_id_2535","f_id_2535",
          "m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","f_id_2535",
          "m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","f_id_2535",
          "m_id_2535","f_id_2535","m_id_2535","f_id_2535")
vars07<-c("diet_mom","diet_dad",
          "diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_mom","diet_mom","diet_dad",
          "diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_dad",
          "diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_dad",
          "diet_mom","diet_dad","diet_mom","diet_dad")
vars08<-c("DASH diet score, women","DASH diet score, male partners",
          "Fruit intake (g/d), women","Fruit intake (g/d), male partners",
          "Vegetable intake (g/d), women","Vegetable intake (g/d), male partners",
          "Low-pesticide fruits and vegetables (g/d), women","High-pesticide fruits and vegetables (g/d), women",
          "Nut and legume intake (g/d), women","Nut and legume intake (g/d), male partners",
          "Whole grain intake (g/d), women","Whole grain intake (g/d), male partners",
          "Refined grain intake (g/d), women","Refined grain intake (g/d), male partners",
          "Coffee and tea intake (mL/d), women","Coffee and tea intake (mL/d), male partners",
          "Sweetened beverage intake (mL/d), women","Sweetened beverage intake (mL/d), male partners",
          "Total dairy intake (g/d), women","Total dairy intake (g/d), male partners",
          "Low-fat dairy intake (g/d), women","Low-fat dairy intake (g/d), male partners",
          "Fish intake (g/d), women","Fish intake (g/d), male partners",
          "Red and processed meat intake (g/d), women","Red and processed meat intake (g/d), male partners",
          "Poultry intake (g/d), women","Poultry intake (g/d), male partners",
          "Proxy of snack intake (g/d), women","Proxy of snack intake (g/d), male partners",
          "Egg intake (g/d), women","Egg intake (g/d), male partners")


tab<-NULL
tab2<-NULL
for(i in 1:length(vars01))
  
{
  datx<-dat[dat[,vars07[i]]==1 & dat$preg_plan==1 & !is.na(dat[,vars06[i]]),]
  datx[,vars01[i]]<-as.numeric(with(datx,scale(datx[,vars00[i]])))
  sample<-length(which(!is.na(datx[,vars01[i]])))
  
  mod01<-miceadds::glm.cluster(formula=as.factor(subf)~datx[,vars01[i]],
                               data=datx, cluster=datx[,vars06[i]], family="binomial")
  coef01<-risk_se_ic_guapa2(as.numeric(summary(mod01)[2,1]),as.numeric(summary(mod01)[2,2]))
  pval01<-pval_guapa(as.numeric(summary(mod01)[2,4]))
  or01<-exp(as.numeric(summary(mod01)[2,1]))
  orlo01<-exp(as.numeric(summary(mod01)[2,1])-(z*as.numeric(summary(mod01)[2,2])))
  orhi01<-exp(as.numeric(summary(mod01)[2,1])+(z*as.numeric(summary(mod01)[2,2])))
  
  mod02<-miceadds::glm.cluster(formula=as.factor(subf)~datx[,vars01[i]]
                               +datx[,vars02[i]]+datx[,vars03[i]]+datx[,vars04[i]]+datx[,vars05[i]]+parity,
                               data=datx, cluster=datx[,vars06[i]], family="binomial")
  coef02<-risk_se_ic_guapa2(as.numeric(summary(mod02)[2,1]),as.numeric(summary(mod02)[2,2]))
  pval02<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  or02<-exp(as.numeric(summary(mod02)[2,1]))
  orlo02<-exp(as.numeric(summary(mod02)[2,1])-(z*as.numeric(summary(mod02)[2,2])))
  orhi02<-exp(as.numeric(summary(mod02)[2,1])+(z*as.numeric(summary(mod02)[2,2])))
  
  aaa<-datx[,vars00[i]]
  mod_lin<-glm(formula=as.factor(subf)~aaa
               +datx[,vars02[i]]+datx[,vars03[i]]+datx[,vars04[i]]+datx[,vars05[i]]+parity,
               data=datx, family="binomial")
  mod_nlin<-gam(formula=as.factor(subf)~bs(aaa,df=4)
                +datx[,vars02[i]]+datx[,vars03[i]]+datx[,vars04[i]]+datx[,vars05[i]]+parity,
                data=datx, family="binomial")
  p_nonlin_lrtest<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  tab<-rbind(tab,cbind(coef01,pval01,or01,orlo01,orhi01,coef02,pval02,or02,orlo02,orhi02,p_nonlin_lrtest))
  
  ptemp<-termplot(mod_nlin,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  #value<-closest(temp$x,mean(aaa,na.rm=TRUE))[1] # Centered around the mean
  value<-closest(temp$x,min(aaa,na.rm=TRUE))[1] # Centered around the minimim value
  center<-with(temp, y[x==value])
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./Outputs/Results/",vars00[i],".jpg",sep="")
  labely<-c("Infertility risk (odds ratio, 95% CI)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  plot.data<-subset2(plot.data,"plot.data$lci>=0.1 & plot.data$uci<=3")
  
  p_lin2<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  p_nonlin2<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  p_lin2<-ifelse(p_lin2=="<0.001"," < 0.001",
                 ifelse(p_lin2=="<0.00001"," < 0.00001",paste(" = ",p_lin2,sep="")))
  p_nonlin2<-ifelse(p_nonlin2=="<0.001"," < 0.001",
                    ifelse(p_nonlin2=="<0.00001"," < 0.00001",paste(" = ",p_nonlin2,sep="")))
  leg<-paste("p-value for linearity",p_lin2,
             "\np-value for non-linearity",p_nonlin2,sep="")
  
  figure<-ggplot(data=plot.data, aes_string(x=plot.data$x, y=plot.data$y)) + 
    geom_ribbon(aes_string(ymin=plot.data$lci, ymax=plot.data$uci), alpha=0.25, fill="Black") +
    geom_line(aes_string(x=plot.data$x, y=plot.data$y), color='Black') +
    geom_hline(yintercept=1, linetype=2) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0)) +
    labs(x=vars08[i],y=labely) +
    annotate("text", x=max(plot.data$x,na.rm=TRUE)*0.98, y=max(plot.data$uci,na.rm=TRUE), label=leg, vjust=1, hjust=1, size=6) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())  
  
  ggsave(filename=name, dpi=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
}
      
rownames(tab)<-vars00
write.table(tab,file="./Outputs/Results/results.csv",sep=";",col.names=NA)

