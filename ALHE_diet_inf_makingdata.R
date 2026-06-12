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

dir.create("N:/durable/Project/ALHE_diet_infertility")
dir.create("N:/durable/Project/ALHE_diet_infertility/Data")

memory.limit(35000)
setwd("N:/durable/Project/ALHE_diet_infertility")


### DATABASE CONSTRUCTION ###
#############################

### MoBa VARIABLES ###
######################

dat<-spss.get("N:/durable/RAW/Pregnancy_questionnaires/PDB2535_Skjema1_v12.sav",
              use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(dat)<-tolower(names(dat))

# BMI #

dat$bmi_mom<-dat$aa85/((dat$aa87/100)^2)
dat$bmi_dad<-dat$aa89/((dat$aa88/100)^2)
dat$bmi_mom<-with(dat,ifelse(bmi_mom<13 | bmi_mom>60,NA,bmi_mom))
dat$bmi_dad<-with(dat,ifelse(bmi_dad<13 | bmi_dad>60,NA,bmi_dad))


# HIGHEST EDUCATIONAL LEVEL, COMPLETED OR ONGOING ###

# Mothers #

dat$aa1125<-with(dat,ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125>aa1124),aa1125,
                            ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125<aa1124),NA,aa1125)))
dat$aa1128x<-with(dat,ifelse(aa1128==1 | aa1129==1,1,0))
dat$aa1128x<-with(dat,ifelse(is.na(aa1128x),0,aa1128x))

dat$edu_mom<-with(dat,ifelse(aa1124==1,1,
                             ifelse(aa1124==2,2,
                                    ifelse(aa1124==3,3,
                                           ifelse(aa1124==4,4,
                                                  ifelse(aa1124==5,5,
                                                         ifelse(aa1124==6,6,NA)))))))
dat$edu_mom<-with(dat,ifelse(is.na(aa1125),edu_mom,
                             ifelse(aa1125==1,1,
                                    ifelse(aa1125==2,2,
                                           ifelse(aa1125==3,3,
                                                  ifelse(aa1125==4,4,
                                                         ifelse(aa1125==5,5,
                                                                ifelse(aa1125==6,6,NA))))))))
dat$edu_mom<-with(dat,ifelse(is.na(edu_mom),0,edu_mom))

dat$edu_mom<-with(dat,ifelse(edu_mom==0 & aa1128x==1,3,
                             ifelse((edu_mom>0 & edu_mom<3) & aa1128x==1,3,
                                    ifelse(edu_mom>=3 & aa1128x==1,edu_mom,
                                           ifelse(edu_mom==0 & aa1128x==0,NA,
                                                  ifelse((edu_mom>0 & edu_mom<3) & aa1128x==0,edu_mom,
                                                         ifelse(edu_mom>=3 & aa1128x==0,edu_mom,NA)))))))

dat$eduyears_mom<-with(dat,ifelse(edu_mom==1,10,
                                  ifelse(edu_mom==2,10,
                                         ifelse(edu_mom==3,13,
                                                ifelse(edu_mom==4,13,
                                                       ifelse(edu_mom==5,19,
                                                              ifelse(edu_mom==6,20,
                                                                     ifelse(is.na(edu_mom),NA,NA))))))))

dat$edu_mom<-with(dat,ifelse(edu_mom==1,1,
                             ifelse(edu_mom==2,1,
                                    ifelse(edu_mom==3,2,
                                           ifelse(edu_mom==4,2,
                                                  ifelse(edu_mom==5,3,
                                                         ifelse(edu_mom==6,4,
                                                                ifelse(is.na(edu_mom),NA,NA))))))))

# Fathers #

dat$aa1127<-with(dat,ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127>aa1126),aa1127,
                            ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127<aa1126),NA,aa1127)))
dat$aa1130x<-with(dat,ifelse(aa1130==1 | aa1131==1,1,0))
dat$aa1130x<-with(dat,ifelse(is.na(aa1130x),0,aa1130x))

dat$edu_dad<-with(dat,ifelse(aa1126==1,1,
                             ifelse(aa1126==2,2,
                                    ifelse(aa1126==3,3,
                                           ifelse(aa1126==4,4,
                                                  ifelse(aa1126==5,5,
                                                         ifelse(aa1126==6,6,NA)))))))
dat$edu_dad<-with(dat,ifelse(is.na(aa1127),edu_dad,
                             ifelse(aa1127==1,1,
                                    ifelse(aa1127==2,2,
                                           ifelse(aa1127==3,3,
                                                  ifelse(aa1127==4,4,
                                                         ifelse(aa1127==5,5,
                                                                ifelse(aa1127==6,6,NA))))))))
dat$edu_dad<-with(dat,ifelse(is.na(edu_dad),0,edu_dad))

dat$edu_dad<-with(dat,ifelse(edu_dad==0 & aa1130x==1,3,
                             ifelse((edu_dad>0 & edu_dad<3) & aa1130x==1,3,
                                    ifelse(edu_dad>=3 & aa1130x==1,edu_dad,
                                           ifelse(edu_dad==0 & aa1130x==0,NA,
                                                  ifelse((edu_dad>0 & edu_dad<3) & aa1130x==0,edu_dad,
                                                         ifelse(edu_dad>=3 & aa1130x==0,edu_dad,NA)))))))

dat$eduyears_dad<-with(dat,ifelse(edu_dad==1,10,
                                  ifelse(edu_dad==2,10,
                                         ifelse(edu_dad==3,13,
                                                ifelse(edu_dad==4,13,
                                                       ifelse(edu_dad==5,19,
                                                              ifelse(edu_dad==6,20,
                                                                     ifelse(is.na(edu_dad),NA,NA))))))))
dat$edu_dad<-with(dat,ifelse(edu_dad==1,1,
                             ifelse(edu_dad==2,1,
                                    ifelse(edu_dad==3,2,
                                           ifelse(edu_dad==4,2,
                                                  ifelse(edu_dad==5,3,
                                                         ifelse(edu_dad==6,4,
                                                                ifelse(is.na(edu_dad),NA,NA))))))))


# PARITY (number of previous deliveries) #

vars01<-c("aa95","aa101","aa107","aa113","aa119","aa125","aa131","aa137","aa143","aa149")

for(i in 1:length(vars01))
{
  dat[,vars01[i]]<-with(dat,ifelse(is.na(dat[,vars01[i]]),0,dat[,vars01[i]]))
  dat[,vars01[i]]<-with(dat,ifelse(dat[,vars01[i]]>0,1,dat[,vars01[i]]))
}

dat$parity<-with(dat,aa95+aa101+aa107+aa113+aa119+aa125+aa131+aa137+aa143+aa149)


# SMOKING #

dat$smkinit_mom<-with(dat,ifelse(aa1355==2 | aa1356>1 | aa1357!=0 | aa1358!=0 | aa1359>1 | aa1360!=0 | aa1361!=0 | 
                                   aa1362!=0 | aa1363!=0 | aa1364!=0 | aa1979>1,1,0))
dat$smkinit_mom<-with(dat,ifelse(is.na(smkinit_mom),0,smkinit_mom))
dat$smkinit_mom<-with(dat,ifelse(is.na(aa1355) & is.na(aa1356) & is.na(aa1357) & is.na(aa1358) & is.na(aa1359) & is.na(aa1360) & is.na(aa1361) & 
                                   is.na(aa1362) & is.na(aa1363) & is.na(aa1364) & is.na(aa1979),NA,smkinit_mom))

dad<-spss.get("N:/durable/RAW/Pregnancy_questionnaires/PDB2535_SkjemaFar_v12.sav",
              use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(dad)<-tolower(names(dad))
dad<-dad[,c("preg_id_2535","ff214","ff215","ff216","ff217","ff218","ff219","ff220")]
dat<-merge2(dat,dad,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)

dat$smkinit_dad<-with(dat,ifelse(aa1353==2 | aa1354==2 | ff214==2 | ff215>1 | ff216!=0 | ff217!=0 | ff218>1 | ff219!=0 | ff220!=0,1,0))
dat$smkinit_dad<-with(dat,ifelse(is.na(smkinit_dad),0,smkinit_dad))
dat$smkinit_dad<-with(dat,ifelse(is.na(aa1353) & is.na(aa1354) & is.na(ff214) & is.na(ff215) & is.na(ff216) & is.na(ff217) & is.na(ff218) & is.na(ff219) & is.na(ff220),NA,smkinit_dad))
dad<-NULL


# MERGE WITH BIRTH REGISTRY #
#############################

mbrnx<-spss.get("N:/durable/RAW/Birth_registry/PDB2535_MFR_541_v12_extra_vars.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(mbrnx)<-tolower(names(mbrnx))
mbrnx<-mbrnx[,c("preg_id_2535","mors_alder_desimalt","fars_alder","dodkat","dodfodte","spabort_12","spabort_23")]
mbrnx$mors_alder<-floor(mbrnx$mors_alder_desimalt)

mbrn<-spss.get("N:/durable/RAW/Birth_registry/PDB2535_MFR_541_v12.sav",
                use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(mbrn)<-tolower(names(mbrn))
mbrn<-merge2(mbrn,mbrnx,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
mbrn$flerfodsel<-with(mbrn,ifelse(is.na(flerfodsel),0,flerfodsel))
mbrn<-subset2(mbrn,"mbrn$flerfodsel==0")


# PREECLAMPSIA / ECLAMPSIA #

vars<-c("preekl","preekltidl","eklampsi","hellp","hypertensjon_alene","hypertensjon_kronisk","svlen","fstart","diabetes_mellitus",
        "spabort_12_5","spabort_23_5","dodfodte","dodkat")
for(i in 1:length(vars))
{
  mbrn[,vars[i]]<-with(mbrn,ifelse(is.na(mbrn[,vars[i]]),0,mbrn[,vars[i]]))
}

mbrn$eclampsia<-with(mbrn,ifelse(preekl>0,1,
                                 ifelse(preekltidl>0,1,
                                        ifelse(eklampsi>0,1,
                                               ifelse(hellp>0,1,0)))))
mbrn$htapreg<-with(mbrn,ifelse(preekl>0,1,
                               ifelse(preekltidl>0,1,
                                      ifelse(eklampsi>0,1,
                                             ifelse(hellp>0,1,
                                                    ifelse(hypertensjon_alene>0,1,0))))))
mbrn$hta_chronic<-mbrn$hypertensjon_kronisk


# PRETERM BIRTH #                                               

mbrn$preterm<-with(mbrn,ifelse(svlen<37,1,0))
mbrn$preterm_sp<-with(mbrn,ifelse(preterm==1 & fstart==1,1,
                                  ifelse(preterm==1 & fstart==2,NA,
                                         ifelse(preterm==1 & fstart==3,NA,
                                                ifelse(preterm==0,0,NA)))))

# STILLBIRTH #

mbrn$stillbirth_index<-with(mbrn,ifelse(dodkat==7 & svlen>=23,1,
                                        ifelse(dodkat==8 & svlen>=23,1,
                                               ifelse(dodkat==9 & svlen>=23,1,0))))
mbrn$stillbirth_index<-with(mbrn,ifelse(is.na(stillbirth_index),9,stillbirth_index))
mbrn$stillbirth_history<-with(mbrn,ifelse(dodfodte>0,1,0))
mbrn$stillbirth_history<-with(mbrn,ifelse(is.na(stillbirth_history),9,stillbirth_history))
mbrn$stillbirth<-with(mbrn,ifelse(stillbirth_index==1 & stillbirth_history==1,1,
                                  ifelse(stillbirth_index==1 & stillbirth_history==0,1,
                                         ifelse(stillbirth_index==1 & stillbirth_history==9,1,
                                                ifelse(stillbirth_index==0 & stillbirth_history==1,1,
                                                       ifelse(stillbirth_index==0 & stillbirth_history==0,0,
                                                              ifelse(stillbirth_index==0 & stillbirth_history==9,0,
                                                                     ifelse(stillbirth_index==9 & stillbirth_history==1,1,
                                                                            ifelse(stillbirth_index==9 & stillbirth_history==0,0,
                                                                                   ifelse(stillbirth_index==9 & stillbirth_history==9,9,NA))))))))))
mbrn$stillbirth<-with(mbrn,ifelse(stillbirth==9,NA,stillbirth))


# MISCARRIAGE #

mbrn$miscarriage_index<-with(mbrn,ifelse(dodkat==7 & svlen<23,1,
                                        ifelse(dodkat==8 & svlen<23,1,
                                               ifelse(dodkat==9 & svlen<23,1,0))))
mbrn$miscarriage_index<-with(mbrn,ifelse(is.na(miscarriage_index),9,miscarriage_index))
mbrn$spabort_12_5<-with(mbrn,ifelse(spabort_12_5>0,1,spabort_12_5))
mbrn$spabort_12_5<-with(mbrn,ifelse(is.na(spabort_12_5),9,spabort_12_5))
mbrn$spabort_23_5<-with(mbrn,ifelse(spabort_23_5>0,1,spabort_23_5))
mbrn$spabort_23_5<-with(mbrn,ifelse(is.na(spabort_23_5),9,spabort_23_5))
mbrn$miscarriage_history<-with(mbrn,ifelse(spabort_12_5==1 & spabort_23_5==1,1,
                                   ifelse(spabort_12_5==1 & spabort_23_5==0,1,
                                          ifelse(spabort_12_5==1 & spabort_23_5==9,1,
                                                 ifelse(spabort_12_5==0 & spabort_23_5==1,1,
                                                        ifelse(spabort_12_5==0 & spabort_23_5==0,0,
                                                               ifelse(spabort_12_5==0 & spabort_23_5==9,0,
                                                                      ifelse(spabort_12_5==9 & spabort_23_5==1,1,
                                                                             ifelse(spabort_12_5==9 & spabort_23_5==0,0,
                                                                                    ifelse(spabort_12_5==9 & spabort_23_5==9,9,NA))))))))))
mbrn$miscarriage_history<-with(mbrn,ifelse(is.na(miscarriage_history),9,miscarriage_history))
mbrn$miscarriage<-with(mbrn,ifelse(miscarriage_index==1 & miscarriage_history==1,1,
                                  ifelse(miscarriage_index==1 & miscarriage_history==0,1,
                                         ifelse(miscarriage_index==1 & miscarriage_history==9,1,
                                                ifelse(miscarriage_index==0 & miscarriage_history==1,1,
                                                       ifelse(miscarriage_index==0 & miscarriage_history==0,0,
                                                              ifelse(miscarriage_index==0 & miscarriage_history==9,0,
                                                                     ifelse(miscarriage_index==9 & miscarriage_history==1,1,
                                                                            ifelse(miscarriage_index==9 & miscarriage_history==0,0,
                                                                                   ifelse(miscarriage_index==9 & miscarriage_history==9,9,NA))))))))))
mbrn$miscarriage<-with(mbrn,ifelse(miscarriage==9,NA,miscarriage))


# GESTATIONAL DIABETES #                                               

mbrn$gdm<-with(mbrn,ifelse(diabetes_mellitus==4,1,0))


# SMALL FOR GESTATIONAL AGE (SGA) # 

# Erase children with sex not specified (0), uncertain (3) or missing (9) #

mbrn$kjonn<-with(mbrn,ifelse(kjonn==1,1,
                             ifelse(kjonn==2,2,NA)))
mbrn<-subset2(mbrn,"!is.na(mbrn$kjonn)")

mbrn$svlen2<-with(mbrn,ifelse(svlen<25,24,
                              ifelse(svlen>42,43,svlen)))
mbrn$vekt<-with(mbrn,ifelse(vekt<500,NA,
                            ifelse(vekt>6500,NA,vekt)))

mbrn$barn_id<-paste(mbrn$preg_id_2535,mbrn$barn_nr,sep="_")

mbrn_boys<-subset2(mbrn,"mbrn$kjonn==1")
mbrn_boys<-mbrn_boys[,c("barn_id","svlen2","vekt")]
mbrn_girls<-subset2(mbrn,"mbrn$kjonn==2")
mbrn_girls<-mbrn_girls[,c("barn_id","svlen2","vekt")]

# SGA in boys #

mbrn_boys$s24<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==24,mbrn_boys$vekt,NA))
mbrn_boys$s25<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==25,mbrn_boys$vekt,NA))
mbrn_boys$s26<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==26,mbrn_boys$vekt,NA))
mbrn_boys$s27<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==27,mbrn_boys$vekt,NA))
mbrn_boys$s28<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==28,mbrn_boys$vekt,NA))
mbrn_boys$s29<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==29,mbrn_boys$vekt,NA))
mbrn_boys$s30<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==30,mbrn_boys$vekt,NA))
mbrn_boys$s31<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==31,mbrn_boys$vekt,NA))
mbrn_boys$s32<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==32,mbrn_boys$vekt,NA))
mbrn_boys$s33<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==33,mbrn_boys$vekt,NA))
mbrn_boys$s34<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==34,mbrn_boys$vekt,NA))
mbrn_boys$s35<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==35,mbrn_boys$vekt,NA))
mbrn_boys$s36<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==36,mbrn_boys$vekt,NA))
mbrn_boys$s37<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==37,mbrn_boys$vekt,NA))
mbrn_boys$s38<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==38,mbrn_boys$vekt,NA))
mbrn_boys$s39<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==39,mbrn_boys$vekt,NA))
mbrn_boys$s40<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==40,mbrn_boys$vekt,NA))
mbrn_boys$s41<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==41,mbrn_boys$vekt,NA))
mbrn_boys$s42<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==42,mbrn_boys$vekt,NA))
mbrn_boys$s43<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==43,mbrn_boys$vekt,NA))

mbrn_boys$sga24<-as.numeric(ntile(mbrn_boys$s24,10))
mbrn_boys$sga25<-as.numeric(ntile(mbrn_boys$s25,10))
mbrn_boys$sga26<-as.numeric(ntile(mbrn_boys$s26,10))
mbrn_boys$sga27<-as.numeric(ntile(mbrn_boys$s27,10))
mbrn_boys$sga28<-as.numeric(ntile(mbrn_boys$s28,10))
mbrn_boys$sga29<-as.numeric(ntile(mbrn_boys$s29,10))
mbrn_boys$sga30<-as.numeric(ntile(mbrn_boys$s30,10))
mbrn_boys$sga31<-as.numeric(ntile(mbrn_boys$s31,10))
mbrn_boys$sga32<-as.numeric(ntile(mbrn_boys$s32,10))
mbrn_boys$sga33<-as.numeric(ntile(mbrn_boys$s33,10))
mbrn_boys$sga34<-as.numeric(ntile(mbrn_boys$s34,10))
mbrn_boys$sga35<-as.numeric(ntile(mbrn_boys$s35,10))
mbrn_boys$sga36<-as.numeric(ntile(mbrn_boys$s36,10))
mbrn_boys$sga37<-as.numeric(ntile(mbrn_boys$s37,10))
mbrn_boys$sga38<-as.numeric(ntile(mbrn_boys$s38,10))
mbrn_boys$sga39<-as.numeric(ntile(mbrn_boys$s39,10))
mbrn_boys$sga40<-as.numeric(ntile(mbrn_boys$s40,10))
mbrn_boys$sga41<-as.numeric(ntile(mbrn_boys$s41,10))
mbrn_boys$sga42<-as.numeric(ntile(mbrn_boys$s42,10))
mbrn_boys$sga43<-as.numeric(ntile(mbrn_boys$s43,10))

vars<-c("sga24","sga25","sga26","sga27","sga28","sga29","sga30","sga31","sga32","sga33",
        "sga34","sga35","sga36","sga37","sga38","sga39","sga40","sga41","sga42","sga43")

for(i in 1:length(vars))
  
{
  mbrn_boys[,vars[i]]<-with(mbrn_boys,ifelse(is.na(mbrn_boys[,vars[i]]),0,mbrn_boys[,vars[i]]))
}


mbrn_boys$sga_boys<-with(mbrn_boys,ifelse(sga24==1,1,
                                          ifelse(sga25==1,1,
                                                 ifelse(sga26==1,1,
                                                        ifelse(sga27==1,1,
                                                               ifelse(sga28==1,1,
                                                                      ifelse(sga29==1,1,
                                                                             ifelse(sga30==1,1,
                                                                                    ifelse(sga31==1,1,
                                                                                           ifelse(sga32==1,1,
                                                                                                  ifelse(sga33==1,1,
                                                                                                         ifelse(sga34==1,1,
                                                                                                                ifelse(sga35==1,1,
                                                                                                                       ifelse(sga36==1,1,
                                                                                                                              ifelse(sga37==1,1,
                                                                                                                                     ifelse(sga38==1,1,
                                                                                                                                            ifelse(sga39==1,1,
                                                                                                                                                   ifelse(sga40==1,1,
                                                                                                                                                          ifelse(sga41==1,1,
                                                                                                                                                                 ifelse(sga42==1,1,
                                                                                                                                                                        ifelse(sga43==1,1,0)))))))))))))))))))))
mbrn_boys<-mbrn_boys[,c("barn_id","sga_boys")]

# SGA in girls #

mbrn_girls$s24<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==24,mbrn_girls$vekt,NA))
mbrn_girls$s25<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==25,mbrn_girls$vekt,NA))
mbrn_girls$s26<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==26,mbrn_girls$vekt,NA))
mbrn_girls$s27<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==27,mbrn_girls$vekt,NA))
mbrn_girls$s28<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==28,mbrn_girls$vekt,NA))
mbrn_girls$s29<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==29,mbrn_girls$vekt,NA))
mbrn_girls$s30<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==30,mbrn_girls$vekt,NA))
mbrn_girls$s31<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==31,mbrn_girls$vekt,NA))
mbrn_girls$s32<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==32,mbrn_girls$vekt,NA))
mbrn_girls$s33<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==33,mbrn_girls$vekt,NA))
mbrn_girls$s34<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==34,mbrn_girls$vekt,NA))
mbrn_girls$s35<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==35,mbrn_girls$vekt,NA))
mbrn_girls$s36<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==36,mbrn_girls$vekt,NA))
mbrn_girls$s37<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==37,mbrn_girls$vekt,NA))
mbrn_girls$s38<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==38,mbrn_girls$vekt,NA))
mbrn_girls$s39<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==39,mbrn_girls$vekt,NA))
mbrn_girls$s40<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==40,mbrn_girls$vekt,NA))
mbrn_girls$s41<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==41,mbrn_girls$vekt,NA))
mbrn_girls$s42<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==42,mbrn_girls$vekt,NA))
mbrn_girls$s43<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==43,mbrn_girls$vekt,NA))

mbrn_girls$sga24<-as.numeric(ntile(mbrn_girls$s24,10))
mbrn_girls$sga25<-as.numeric(ntile(mbrn_girls$s25,10))
mbrn_girls$sga26<-as.numeric(ntile(mbrn_girls$s26,10))
mbrn_girls$sga27<-as.numeric(ntile(mbrn_girls$s27,10))
mbrn_girls$sga28<-as.numeric(ntile(mbrn_girls$s28,10))
mbrn_girls$sga29<-as.numeric(ntile(mbrn_girls$s29,10))
mbrn_girls$sga30<-as.numeric(ntile(mbrn_girls$s30,10))
mbrn_girls$sga31<-as.numeric(ntile(mbrn_girls$s31,10))
mbrn_girls$sga32<-as.numeric(ntile(mbrn_girls$s32,10))
mbrn_girls$sga33<-as.numeric(ntile(mbrn_girls$s33,10))
mbrn_girls$sga34<-as.numeric(ntile(mbrn_girls$s34,10))
mbrn_girls$sga35<-as.numeric(ntile(mbrn_girls$s35,10))
mbrn_girls$sga36<-as.numeric(ntile(mbrn_girls$s36,10))
mbrn_girls$sga37<-as.numeric(ntile(mbrn_girls$s37,10))
mbrn_girls$sga38<-as.numeric(ntile(mbrn_girls$s38,10))
mbrn_girls$sga39<-as.numeric(ntile(mbrn_girls$s39,10))
mbrn_girls$sga40<-as.numeric(ntile(mbrn_girls$s40,10))
mbrn_girls$sga41<-as.numeric(ntile(mbrn_girls$s41,10))
mbrn_girls$sga42<-as.numeric(ntile(mbrn_girls$s42,10))
mbrn_girls$sga43<-as.numeric(ntile(mbrn_girls$s43,10))

vars<-c("sga24","sga25","sga26","sga27","sga28","sga29","sga30","sga31","sga32","sga33",
        "sga34","sga35","sga36","sga37","sga38","sga39","sga40","sga41","sga42","sga43")

for(i in 1:length(vars))
  
{
  mbrn_girls[,vars[i]]<-with(mbrn_girls,ifelse(is.na(mbrn_girls[,vars[i]]),0,mbrn_girls[,vars[i]]))
}


mbrn_girls$sga_girls<-with(mbrn_girls,ifelse(sga24==1,1,
                                             ifelse(sga25==1,1,
                                                    ifelse(sga26==1,1,
                                                           ifelse(sga27==1,1,
                                                                  ifelse(sga28==1,1,
                                                                         ifelse(sga29==1,1,
                                                                                ifelse(sga30==1,1,
                                                                                       ifelse(sga31==1,1,
                                                                                              ifelse(sga32==1,1,
                                                                                                     ifelse(sga33==1,1,
                                                                                                            ifelse(sga34==1,1,
                                                                                                                   ifelse(sga35==1,1,
                                                                                                                          ifelse(sga36==1,1,
                                                                                                                                 ifelse(sga37==1,1,
                                                                                                                                        ifelse(sga38==1,1,
                                                                                                                                               ifelse(sga39==1,1,
                                                                                                                                                      ifelse(sga40==1,1,
                                                                                                                                                             ifelse(sga41==1,1,
                                                                                                                                                                    ifelse(sga42==1,1,
                                                                                                                                                                           ifelse(sga43==1,1,0)))))))))))))))))))))

mbrn_girls<-mbrn_girls[,c("barn_id","sga_girls")]

mbrn<-merge2(mbrn,mbrn_boys,by.id=c("barn_id"),all.x=TRUE,sort=FALSE)
mbrn<-merge2(mbrn,mbrn_girls,by.id=c("barn_id"),all.x=TRUE,sort=FALSE)

vars<-c("sga_boys","sga_girls")
for(i in 1:length(vars))
{
  mbrn[,vars[i]]<-with(mbrn,ifelse(is.na(mbrn[,vars[i]]),0,mbrn[,vars[i]]))
}

mbrn$sga<-with(mbrn,ifelse(sga_boys==1,1,
                           ifelse(sga_girls==1,1,0)))
mbrn_boys<-NULL
mbrn_girls<-NULL
mbrn<-rename.vars(mbrn,from=c("mors_alder","fars_alder"),to=c("agedelivery_mom","agedelivery_dad"))
mbrn$agedelivery_dad<-with(mbrn,ifelse(agedelivery_dad>100,NA,agedelivery_dad))

mbrn<-mbrn[,c("preg_id_2535","agedelivery_mom","agedelivery_dad","art","flerfodsel",
              "eclampsia","htapreg","hta_chronic","preterm","preterm_sp","gdm","sga","miscarriage","stillbirth")]
dat<-merge2(dat,mbrn,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
attributes(dat$agedelivery_mom)$label<-c("Age, women")
attributes(dat$agedelivery_dad)$label<-c("Age, men")
dat$flerfodsel<-with(dat,ifelse(is.na(flerfodsel),0,flerfodsel))
dat<-dat[dat$flerfodsel==0,]
dad<-NULL
mbrn<-NULL
mbrnx<-NULL


# INFERTILITY #

dat$preg_plan<-with(dat,ifelse(is.na(aa46),NA,
                               ifelse(aa46==0,NA,
                                      ifelse(aa46==1,0,
                                             ifelse(aa46==2,1,NA)))))
attributes(dat$preg_plan)$value.label<-c("0=No","1=Yes")
attributes(dat$preg_plan)$label<-c("Planned pregnancy")

dat$ttp<-with(dat,ifelse(is.na(aa48) & !is.na(preg_plan),0,
                         ifelse(is.na(aa48) & is.na(preg_plan),NA,
                                ifelse(!is.na(aa48) & !is.na(preg_plan),aa48,
                                       ifelse(!is.na(aa48) & is.na(preg_plan),NA,NA)))))

dat$art<-with(dat,ifelse(is.na(art),0,
                         ifelse(art>0,1,NA)))
dat$subf<-with(dat,ifelse(ttp<12,0,
                          ifelse(ttp>=12,1,NA)))
dat$subf<-with(dat,ifelse(art==0,subf,
                          ifelse(art==1,1,NA)))

dat<-dat[,c("preg_id_2535","agedelivery_mom","agedelivery_dad","bmi_mom","bmi_dad","eduyears_mom","eduyears_dad","smkinit_mom","smkinit_dad",
            "parity","preg_plan","subf","art",
            "eclampsia","htapreg","hta_chronic","preterm","preterm_sp","gdm","sga","miscarriage","stillbirth")]


# MERGE WITH KEY FILE (m_id AND f_id) #

keys<-spss.get("N:/durable/RAW/Valid_consent/PDB2535_SV_INFO_V12_20230417.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(keys)<-tolower(names(keys))
dat<-merge2(keys,dat,by.id=c("preg_id_2535"),all.y=TRUE,sort=FALSE)
dat$m_id_2535<-gsub(" ", "", dat$m_id_2535)
dat$f_id_2535<-gsub(" ", "", dat$f_id_2535)


# Pregnancy outcomes (unique variable for mothers and fathers, all pregnancies) #

mom_agemax<-dat[,c("m_id_2535","agedelivery_mom")]
mom_agemax<-mom_agemax[order(mom_agemax$m_id_2535,-abs(mom_agemax$agedelivery_mom)),]
mom_agemax<-mom_agemax[!duplicated(mom_agemax$m_id_2535),]

mom_bmimax<-dat[,c("m_id_2535","bmi_mom")]
mom_bmimax<-mom_bmimax[order(mom_bmimax$m_id_2535,-abs(mom_bmimax$bmi_mom)),]
mom_bmimax<-mom_bmimax[!duplicated(mom_bmimax$m_id_2535),]

mom_eduyearsmax<-dat[,c("m_id_2535","eduyears_mom")]
mom_eduyearsmax<-mom_eduyearsmax[order(mom_eduyearsmax$m_id_2535,-abs(mom_eduyearsmax$eduyears_mom)),]
mom_eduyearsmax<-mom_eduyearsmax[!duplicated(mom_eduyearsmax$m_id_2535),]

mom_smokingmax<-dat[,c("m_id_2535","smkinit_mom")]
mom_smokingmax<-mom_smokingmax[order(mom_smokingmax$m_id_2535,-abs(mom_smokingmax$smkinit_mom)),]
mom_smokingmax<-mom_smokingmax[!duplicated(mom_smokingmax$m_id_2535),]

mom_paritymax<-dat[,c("m_id_2535","parity")]
mom_paritymax<-mom_paritymax[order(mom_paritymax$m_id_2535,-abs(mom_paritymax$parity)),]
mom_paritymax<-mom_paritymax[!duplicated(mom_paritymax$m_id_2535),]

mom_eclampsia<-dat[,c("m_id_2535","eclampsia")]
mom_eclampsia<-mom_eclampsia[order(mom_eclampsia$m_id_2535,-abs(mom_eclampsia$eclampsia)),]
mom_eclampsia<-mom_eclampsia[!duplicated(mom_eclampsia$m_id_2535),]

mom_htapreg<-dat[,c("m_id_2535","htapreg")]
mom_htapreg<-mom_htapreg[order(mom_htapreg$m_id_2535,-abs(mom_htapreg$htapreg)),]
mom_htapreg<-mom_htapreg[!duplicated(mom_htapreg$m_id_2535),]

mom_hta_chronic<-dat[,c("m_id_2535","hta_chronic")]
mom_hta_chronic<-mom_hta_chronic[order(mom_hta_chronic$m_id_2535,-abs(mom_hta_chronic$hta_chronic)),]
mom_hta_chronic<-mom_hta_chronic[!duplicated(mom_hta_chronic$m_id_2535),]

mom_preterm<-dat[,c("m_id_2535","preterm")]
mom_preterm<-mom_preterm[order(mom_preterm$m_id_2535,-abs(mom_preterm$preterm)),]
mom_preterm<-mom_preterm[!duplicated(mom_preterm$m_id_2535),]

mom_preterm_sp<-dat[,c("m_id_2535","preterm_sp")]
mom_preterm_sp<-mom_preterm_sp[order(mom_preterm_sp$m_id_2535,-abs(mom_preterm_sp$preterm_sp)),]
mom_preterm_sp<-mom_preterm_sp[!duplicated(mom_preterm_sp$m_id_2535),]

mom_gdm<-dat[,c("m_id_2535","gdm")]
mom_gdm<-mom_gdm[order(mom_gdm$m_id_2535,-abs(mom_gdm$gdm)),]
mom_gdm<-mom_gdm[!duplicated(mom_gdm$m_id_2535),]

mom_sga<-dat[,c("m_id_2535","sga")]
mom_sga<-mom_sga[order(mom_sga$m_id_2535,-abs(mom_sga$sga)),]
mom_sga<-mom_sga[!duplicated(mom_sga$m_id_2535),]

mom_miscarriage<-dat[,c("m_id_2535","miscarriage")]
mom_miscarriage<-mom_miscarriage[order(mom_miscarriage$m_id_2535,-abs(mom_miscarriage$miscarriage)),]
mom_miscarriage<-mom_miscarriage[!duplicated(mom_miscarriage$m_id_2535),]

mom_stillbirth<-dat[,c("m_id_2535","stillbirth")]
mom_stillbirth<-mom_stillbirth[order(mom_stillbirth$m_id_2535,-abs(mom_stillbirth$stillbirth)),]
mom_stillbirth<-mom_stillbirth[!duplicated(mom_stillbirth$m_id_2535),]

mom_subf<-dat[,c("m_id_2535","subf")]
mom_subf<-mom_subf[order(mom_subf$m_id_2535,-abs(mom_subf$subf)),]
mom_subf<-mom_subf[!duplicated(mom_subf$m_id_2535),]

mom_art<-dat[,c("m_id_2535","art")]
mom_art<-mom_art[order(mom_art$m_id_2535,-abs(mom_art$art)),]
mom_art<-mom_art[!duplicated(mom_art$m_id_2535),]

dat_mom<-dat[,c("preg_id_2535","m_id_2535")]
dat_mom<-merge2(dat_mom,mom_agemax,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_bmimax,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_eduyearsmax,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_smokingmax,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_paritymax,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_htapreg,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_eclampsia,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_hta_chronic,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_preterm,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_preterm_sp,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_gdm,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_sga,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_miscarriage,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_stillbirth,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_subf,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-merge2(dat_mom,mom_art,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
names(dat_mom)<-c("m_id_2535","preg_id_2535","agemax_mom","bmimax_mom","eduyearsmax_mom","smkinitmax_mom","paritymax_mom",
                  "htapreg_mom","eclampsia_mom","hta_chronic_mom","preterm_mom","preterm_sp_mom",
                  "gdm_mom","sga_mom","miscarriage_mom","stillbirth_mom","subf_mom","art_mom")


dad_agemax<-dat[,c("f_id_2535","agedelivery_dad")]
dad_agemax<-dad_agemax[order(dad_agemax$f_id_2535,-abs(dad_agemax$agedelivery_dad)),]
dad_agemax<-dad_agemax[!duplicated(dad_agemax$f_id_2535),]

dad_bmimax<-dat[,c("f_id_2535","bmi_dad")]
dad_bmimax<-dad_bmimax[order(dad_bmimax$f_id_2535,-abs(dad_bmimax$bmi_dad)),]
dad_bmimax<-dad_bmimax[!duplicated(dad_bmimax$f_id_2535),]

dad_eduyearsmax<-dat[,c("f_id_2535","eduyears_dad")]
dad_eduyearsmax<-dad_eduyearsmax[order(dad_eduyearsmax$f_id_2535,-abs(dad_eduyearsmax$eduyears_dad)),]
dad_eduyearsmax<-dad_eduyearsmax[!duplicated(dad_eduyearsmax$f_id_2535),]

dad_smokingmax<-dat[,c("f_id_2535","smkinit_dad")]
dad_smokingmax<-dad_smokingmax[order(dad_smokingmax$f_id_2535,-abs(dad_smokingmax$smkinit_dad)),]
dad_smokingmax<-dad_smokingmax[!duplicated(dad_smokingmax$f_id_2535),]

dad_paritymax<-dat[,c("f_id_2535","parity")]
dad_paritymax<-dad_paritymax[order(dad_paritymax$f_id_2535,-abs(dad_paritymax$parity)),]
dad_paritymax<-dad_paritymax[!duplicated(dad_paritymax$f_id_2535),]

dad_eclampsia<-dat[,c("f_id_2535","eclampsia")]
dad_eclampsia<-dad_eclampsia[order(dad_eclampsia$f_id_2535,-abs(dad_eclampsia$eclampsia)),]
dad_eclampsia<-dad_eclampsia[!duplicated(dad_eclampsia$f_id_2535),]

dad_htapreg<-dat[,c("f_id_2535","htapreg")]
dad_htapreg<-dad_htapreg[order(dad_htapreg$f_id_2535,-abs(dad_htapreg$htapreg)),]
dad_htapreg<-dad_htapreg[!duplicated(dad_htapreg$f_id_2535),]

dad_hta_chronic<-dat[,c("f_id_2535","hta_chronic")]
dad_hta_chronic<-dad_hta_chronic[order(dad_hta_chronic$f_id_2535,-abs(dad_hta_chronic$hta_chronic)),]
dad_hta_chronic<-dad_hta_chronic[!duplicated(dad_hta_chronic$f_id_2535),]

dad_preterm<-dat[,c("f_id_2535","preterm")]
dad_preterm<-dad_preterm[order(dad_preterm$f_id_2535,-abs(dad_preterm$preterm)),]
dad_preterm<-dad_preterm[!duplicated(dad_preterm$f_id_2535),]

dad_preterm_sp<-dat[,c("f_id_2535","preterm_sp")]
dad_preterm_sp<-dad_preterm_sp[order(dad_preterm_sp$f_id_2535,-abs(dad_preterm_sp$preterm_sp)),]
dad_preterm_sp<-dad_preterm_sp[!duplicated(dad_preterm_sp$f_id_2535),]

dad_gdm<-dat[,c("f_id_2535","gdm")]
dad_gdm<-dad_gdm[order(dad_gdm$f_id_2535,-abs(dad_gdm$gdm)),]
dad_gdm<-dad_gdm[!duplicated(dad_gdm$f_id_2535),]

dad_sga<-dat[,c("f_id_2535","sga")]
dad_sga<-dad_sga[order(dad_sga$f_id_2535,-abs(dad_sga$sga)),]
dad_sga<-dad_sga[!duplicated(dad_sga$f_id_2535),]

dad_miscarriage<-dat[,c("f_id_2535","miscarriage")]
dad_miscarriage<-dad_miscarriage[order(dad_miscarriage$f_id_2535,-abs(dad_miscarriage$miscarriage)),]
dad_miscarriage<-dad_miscarriage[!duplicated(dad_miscarriage$f_id_2535),]

dad_stillbirth<-dat[,c("f_id_2535","stillbirth")]
dad_stillbirth<-dad_stillbirth[order(dad_stillbirth$f_id_2535,-abs(dad_stillbirth$stillbirth)),]
dad_stillbirth<-dad_stillbirth[!duplicated(dad_stillbirth$f_id_2535),]

dad_subf<-dat[,c("f_id_2535","subf")]
dad_subf<-dad_subf[order(dad_subf$f_id_2535,-abs(dad_subf$subf)),]
dad_subf<-dad_subf[!duplicated(dad_subf$f_id_2535),]

dad_art<-dat[,c("f_id_2535","art")]
dad_art<-dad_art[order(dad_art$f_id_2535,-abs(dad_art$art)),]
dad_art<-dad_art[!duplicated(dad_art$f_id_2535),]

dat_dad<-dat[,c("preg_id_2535","f_id_2535")]
dat_dad<-merge2(dat_dad,dad_agemax,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_bmimax,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_eduyearsmax,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_smokingmax,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_paritymax,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_htapreg,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_eclampsia,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_hta_chronic,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_preterm,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_preterm_sp,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_gdm,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_sga,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_miscarriage,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_stillbirth,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_subf,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
dat_dad<-merge2(dat_dad,dad_art,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
names(dat_dad)<-c("f_id_2535","preg_id_2535","agemax_dad","bmimax_dad","eduyearsmax_dad","smkinitmax_dad","paritymax_dad",
                  "htapreg_dad","eclampsia_dad","hta_chronic_dad","preterm_dad","preterm_sp_dad",
                  "gdm_dad","sga_dad","miscarriage_dad","stillbirth_dad","subf_dad","art_dad")

dat_mom$m_id_2535<-NULL
dat_dad$f_id_2535<-NULL

dat<-merge2(dat,dat_mom,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,dat_dad,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
dat_mom<-NULL
dat_dad<-NULL


### SENTRIXIDs ###
##################

mom_gen<-spss.get("N:/durable/Raw/Linkage_files_genotype_data/2022_04_25_MoBaGeneticsTot_Mother_PDB2535.sav",
                  use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(mom_gen)<-tolower(names(mom_gen))
mom_gen$m_id_2535<-gsub(" ", "", mom_gen$m_id_2535)
mom_gen$sentrixid_mom<-gsub(" ", "", mom_gen$sentrix_id)
mom_gen<-mom_gen[,c("m_id_2535","sentrixid_mom")]
mom_gen<-mom_gen[!duplicated(mom_gen$m_id_2535),]

dad_gen<-spss.get("N:/durable/Raw/Linkage_files_genotype_data/2022_04_25_MoBaGeneticsTot_Father_PDB2535.sav",
                  use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(dad_gen)<-tolower(names(dad_gen))
dad_gen$f_id_2535<-gsub(" ", "", dad_gen$f_id_2535)
dad_gen$sentrixid_dad<-gsub(" ", "", dad_gen$sentrix_id)
dad_gen<-dad_gen[,c("f_id_2535","sentrixid_dad")]
dad_gen<-dad_gen[!duplicated(dad_gen$f_id_2535),]

dat<-merge2(dat,mom_gen,by.id=c("m_id_2535"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,dad_gen,by.id=c("f_id_2535"),all.x=TRUE,sort=FALSE)
mom_gen<-NULL
dad_gen<-NULL


### FOOD GROUPS ###
###################

# MOTHERS - FFQ #

food<-spss.get("N:/durable/Raw/Pregnancy_questionnaires/PDB2535_Skjema2_beregning_CDW_foodgroups255_v12.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(food)<-tolower(names(food))

ffq<-spss.get("N:/durable/Raw/Pregnancy_questionnaires/PDB2535_Skjema2CDW_v12.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(ffq)<-tolower(names(ffq))
ffq<-ffq[,c("preg_id_2535","bb825","bb826","bb827","bb828","bb829","bb830","bb831","bb832","bb833")]

vars<-c("bb825","bb828","bb831")
for(i in 1:length(vars))
{
  ffq[,vars[i]]<-with(ffq,ifelse(is.na(ffq[,vars[i]]),0,ffq[,vars[i]]))
  ffq[,vars[i]]<-with(ffq,ifelse(ffq[,vars[i]]==0,0,
                                 ifelse(ffq[,vars[i]]==1,4,
                                        ifelse(ffq[,vars[i]]==2,2.5,
                                               ifelse(ffq[,vars[i]]==3,1,NA)))))
}

vars<-c("bb826","bb829","bb832")
for(i in 1:length(vars))
{
  ffq[,vars[i]]<-with(ffq,ifelse(is.na(ffq[,vars[i]]),0,ffq[,vars[i]]))
  ffq[,vars[i]]<-with(ffq,ifelse(ffq[,vars[i]]==0,0,
                                 ifelse(ffq[,vars[i]]==1,5.5/7,
                                        ifelse(ffq[,vars[i]]==2,3.5/7,
                                               ifelse(ffq[,vars[i]]==3,1.5/7,NA)))))
}

vars<-c("bb827","bb830","bb833")
for(i in 1:length(vars))
{
  ffq[,vars[i]]<-with(ffq,ifelse(is.na(ffq[,vars[i]]),0,ffq[,vars[i]]))
  ffq[,vars[i]]<-with(ffq,ifelse(ffq[,vars[i]]==0,0,
                                 ifelse(ffq[,vars[i]]==1,2.5/30,
                                        ifelse(ffq[,vars[i]]==2,1/30,
                                               ifelse(ffq[,vars[i]]==3,0,NA)))))
}

ffq$snack_proxy_mom<-with(ffq,bb825+bb826+bb827+bb828+bb829+bb830+bb831+bb832+bb833)
ffq<-ffq[,c("preg_id_2535","snack_proxy_mom")]
food<-merge2(food,ffq,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)


food$juice_mom<-with(food,spis171+spis172+spis173)
food$juice_fruit_mom<-with(food,ifelse(juice_mom>=250,250,
                                       ifelse(juice_mom<250,juice_mom,NA)))
food$juice_mom<-with(food,ifelse(juice_mom>=250,juice_mom-250,
                                       ifelse(juice_mom<250,0,NA)))
food$fruit_mom<-with(food,spis148+spis149+spis150+spis153+spis154+spis155+spis156+spis157+spis158+spis159+spis160+spis161+spis162
                 +spis163+spis174+spis175+juice_fruit_mom)
food$veg_mom<-with(food,spis120+spis121+spis122+spis123+spis124+spis126+spis127+spis128+spis129+spis130+spis131+spis132+spis133+spis134
               +spis137+spis138+spis139+spis141+spis143+spis145+spis152+spis176+spis181+spis182+spis183+spis184+spis185+spis186
               +spis187+(spis178/3))
food$nuts_mom<-with(food,spis166+spis168+spis169)
food$legume_mom<-with(food,spis125+(spis177/2))
food$nutleg_mom<-with(food,nuts_mom+legume_mom)
food$wholeg_mom<-with(food,spis98+spis101+spis109+spis113+spis118)
food$refg_mom<-with(food,spis97+spis99+spis100+spis102+spis108+spis111+spis114+spis115+spis116
                +(spis59/3)+(spis96/3)+(spis178/3))
food$potato_mom<-with(food,spis144+spis189+spis190)
food$dairy_lf_mom<-with(food,spis02+spis03+spis11+spis13+spis14+spis21+spis22)
food$dairy_tot_mom<-with(food,spis01+spis02+spis03+spis04+spis05+spis06+spis07+spis08+spis09+spis10+spis11+spis12+spis13
                     +spis14+spis15+spis16+spis17+spis18+spis19+spis20+spis21+spis22+spis23+spis24+spis25)
food$dairy_fat_mom<-with(food,dairy_tot_mom-dairy_lf_mom)
food$oils_mom<-with(food,spis194+spis195+spis198+spis201+spis212)
food$cof_tea_mom<-with(food,spis238+spis239+spis245+spis246+spis248+spis253)
food$juice_mom<-with(food,spis171+spis172+spis173)
food$ssb_mom<-food$spis235
food$asb_mom<-food$spis236
food$sbev_mom<-with(food,juice_mom+ssb_mom+asb_mom)
food$sweets_mom<-with(food,spis216+spis217+spis218+spis219+spis220+spis221+spis222+spis223+spis224+spis225+spis226)
food$ssb_sweets_mom<-with(food,ssb_mom+sweets_mom)
food$dessert_mom<-with(food,spis103+spis104+spis105+spis106+spis107+spis110+spis112+spis117)
food$animalfat_mom<-with(food,spis193+(spis214*0.7))
food$eggs_mom<-with(food,spis28+spis29)
food$fish_mom<-with(food,spis71+spis74+spis76+spis79+spis80+spis81+spis82+spis83+spis84+spis85
                +spis86+spis87+spis88+spis89+spis90+spis91+spis92+spis93+spis94+(spis96/9))
food$seafood_mom<-with(food,spis72+spis73+spis75+spis77+spis78+spis95+(spis96*2/9))
food$fish_sf_mom<-with(food,fish_mom+seafood_mom)
food$rpmeat_mom<-with(food,spis30+spis31+spis32+spis33+spis37+spis38+spis39+spis40+spis41+spis42+spis43+spis44+spis45+spis46
                  +spis56+spis57+spis58+(spis59/3)+spis60+spis61+spis62+spis63+spis64+spis65+spis66+spis67+spis68+spis69)
food$wmeat_mom<-with(food,spis48+spis49+spis50+spis51+spis52+spis53+spis54+spis55)
food$tmeat_mom<-with(food,rpmeat_mom+wmeat_mom)
food$misc_animal_mom<-with(food,spis34+spis35+spis36+spis47+spis70)
food$snacks_mom<-with(food,spis119+spis233+spis234)
food$fv_lowpest_mom<-with(food,spis125+spis151+spis165+spis131+spis186+spis152+spis188+spis127
                          +spis184+spis123+spis183+spis159+spis175+spis153+spis122+spis150
                          +spis124+spis187+spis126+spis182+spis129+spis133+spis134)
food$fv_highpest_mom<-with(food,spis139+spis140+spis148+spis143+spis155+spis163+spis145+spis154
                           +spis164+spis144+spis137+spis156+spis162+spis149+spis138+spis132+spis181)

#quantile(food$rpmeat_mom,probs=seq(0,1,0.025),na.rm=TRUE)
#hist(food$cof_tea_mom)
#sort(food$cof_tea_mom,decreasing=TRUE)
food$fruit_mom<-with(food,ifelse(fruit_mom>1000,NA,fruit_mom))
food$veg_mom<-with(food,ifelse(veg_mom>500,NA,veg_mom))
food$wholeg_mom<-with(food,ifelse(wholeg_mom>600,NA,wholeg_mom))
food$refg_mom<-with(food,ifelse(refg_mom>600,NA,refg_mom))
food$dairy_tot_mom<-with(food,ifelse(dairy_tot_mom>1500,NA,dairy_tot_mom))
food$dairy_lf_mom<-with(food,ifelse(dairy_lf_mom>1500,NA,dairy_lf_mom))
food$sbev_mom<-with(food,ifelse(sbev_mom>2000,NA,sbev_mom))
food$rpmeat_mom<-with(food,ifelse(rpmeat_mom>250,NA,rpmeat_mom))
food$fish_mom<-with(food,ifelse(fish_mom>250,NA,fish_mom))
food$cof_tea_mom<-with(food,ifelse(cof_tea_mom>2001,NA,cof_tea_mom))

food$dash_mom_ok<-with(food,ifelse(!is.na(fruit_mom) & !is.na(veg_mom) & !is.na(nutleg_mom) & !is.na(wholeg_mom)
                                   & !is.na(dairy_lf_mom) & !is.na(snack_proxy_mom) & !is.na(sbev_mom) & !is.na(rpmeat_mom),1,0))
food<-food[food$dash_mom_ok==1,]

food$dashq_fruit_mom<-cut(food$fruit_mom,quantile(food$fruit_mom,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
food$dashq_veg_mom<-cut(food$veg_mom,quantile(food$veg_mom,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
food$dashq_nutleg_mom<-cut(food$nutleg_mom,quantile(food$nutleg_mom,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
food$dashq_dairy_lf_mom<-cut(food$dairy_lf_mom,quantile(food$dairy_lf_mom,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
food$dashq_wholeg_mom<-cut(food$wholeg_mom,quantile(food$wholeg_mom,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
food$dashq_snacks_mom<-cut(food$snack_proxy_mom,unique(quantile(food$snack_proxy_mom,probs=seq(0,1,0.2),na.rm=TRUE)),include.lowest=TRUE,labels=FALSE)
food$dashq_snacks_mom<-with(food,ifelse(dashq_snacks_mom==4,1,
                                        ifelse(dashq_snacks_mom==3,2,
                                               ifelse(dashq_snacks_mom==2,3,
                                                      ifelse(dashq_snacks_mom==1,4.5,NA)))))
food$dashq_sbev_mom<-cut(food$sbev_mom,quantile(food$sbev_mom,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
food$dashq_rpmeat_mom<-cut(food$rpmeat_mom,quantile(food$rpmeat_mom,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)


vars<-c("dashq_sbev_mom","dashq_rpmeat_mom")

for(i in 1:length(vars))
{
  food[,vars[i]]<-with(food,ifelse(food[,vars[i]]==1,5,
                                   ifelse(food[,vars[i]]==2,4,
                                          ifelse(food[,vars[i]]==3,3,
                                                 ifelse(food[,vars[i]]==4,2,
                                                        ifelse(food[,vars[i]]==5,1,NA))))))
}

food$dash_mom<-with(food,dashq_fruit_mom+dashq_veg_mom+dashq_nutleg_mom+dashq_dairy_lf_mom
                    +dashq_wholeg_mom+dashq_snacks_mom+dashq_sbev_mom+dashq_rpmeat_mom)

food<-food[,c("preg_id_2535","fruit_mom","veg_mom","nuts_mom","legume_mom","nutleg_mom","wholeg_mom","refg_mom","potato_mom","oils_mom",
              "dairy_lf_mom","dairy_tot_mom","dairy_fat_mom","cof_tea_mom","juice_mom","ssb_mom","asb_mom","sbev_mom","sweets_mom",
              "dessert_mom","snacks_mom","snack_proxy_mom","animalfat_mom","eggs_mom","fish_mom","seafood_mom","fish_sf_mom",
              "rpmeat_mom","wmeat_mom","tmeat_mom","misc_animal_mom","fv_lowpest_mom","fv_highpest_mom","dash_mom")]
food$diet_mom<-1


# FATHERS - QUESTIONS #

dad<-spss.get("N:/durable/RAW/Pregnancy_questionnaires/PDB2535_SkjemaFar_v12.sav",
              use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(dad)<-tolower(names(dad))

dad$diet_dad<-with(dad,ifelse(is.na(dad$ff403) & is.na(dad$ff404) & is.na(dad$ff405) & is.na(dad$ff406)
                             & is.na(dad$ff407) & is.na(dad$ff408) & is.na(dad$ff409)
                             & is.na(dad$ff410) & is.na(dad$ff411) & is.na(dad$ff412) & is.na(dad$ff413) 
                             & is.na(dad$ff414) & is.na(dad$ff415) & is.na(dad$ff416)
                             & is.na(dad$ff417) & is.na(dad$ff418) & is.na(dad$ff419)
                             & is.na(dad$ff420) & is.na(dad$ff421) & is.na(dad$ff422) & is.na(dad$ff423) 
                             & is.na(dad$ff424) & is.na(dad$ff425) & is.na(dad$ff426)
                             & is.na(dad$ff427) & is.na(dad$ff428) & is.na(dad$ff429)
                             & is.na(dad$ff430) & is.na(dad$ff431) & is.na(dad$ff432) & is.na(dad$ff433) 
                             & is.na(dad$ff434) & is.na(dad$ff435) & is.na(dad$ff436)
                             & is.na(dad$ff437) & is.na(dad$ff438) & is.na(dad$ff439)
                             & is.na(dad$ff440) & is.na(dad$ff441) & is.na(dad$ff442) & is.na(dad$ff443) 
                             & is.na(dad$ff444) & is.na(dad$ff445) & is.na(dad$ff446)
                             & is.na(dad$ff447) & is.na(dad$ff448) & is.na(dad$ff449) & is.na(dad$ff450) & is.na(dad$ff451),0,1))
dad<-dad[dad$diet_dad==1,]

vars<-c("ff403","ff404","ff405","ff406")

for(i in 1:length(vars))
{
  dad[,vars[i]]<-with(dad,ifelse(is.na(dad[,vars[i]]),0,dad[,vars[i]]))
}

vars<-c("ff408","ff409","ff410","ff411","ff412","ff413","ff414","ff415","ff416")

for(i in 1:length(vars))
{
  dad[,vars[i]]<-with(dad,ifelse(is.na(dad[,vars[i]]),1,dad[,vars[i]]))
  dad[,vars[i]]<-with(dad,ifelse(dad[,vars[i]]==0,0,
                                 ifelse(dad[,vars[i]]==1,0,
                                        ifelse(dad[,vars[i]]==2,1.5/7,
                                               ifelse(dad[,vars[i]]==3,3.5/7,
                                                      ifelse(dad[,vars[i]]==4,6/7,
                                                             ifelse(dad[,vars[i]]==5,1.5,NA)))))))
}

vars<-c("ff417","ff418","ff419","ff420","ff421","ff422","ff423","ff424","ff425","ff426","ff427")

for(i in 1:length(vars))
{
  dad[,vars[i]]<-with(dad,ifelse(is.na(dad[,vars[i]]),1,dad[,vars[i]]))
  dad[,vars[i]]<-with(dad,ifelse(dad[,vars[i]]==0,0,
                                 ifelse(dad[,vars[i]]==1,0,
                                        ifelse(dad[,vars[i]]==2,0.5,
                                               ifelse(dad[,vars[i]]==3,1,
                                                      ifelse(dad[,vars[i]]==4,2.5,
                                                             ifelse(dad[,vars[i]]==5,4,NA)))))))
}

vars<-c("ff428","ff429","ff430","ff431","ff432","ff433","ff434","ff435","ff436","ff437","ff438")

for(i in 1:length(vars))
{
  dad[,vars[i]]<-with(dad,ifelse(is.na(dad[,vars[i]]),1,dad[,vars[i]]))
  dad[,vars[i]]<-with(dad,ifelse(dad[,vars[i]]==0,0,
                                 ifelse(dad[,vars[i]]==1,0,
                                        ifelse(dad[,vars[i]]==2,0.375/7,
                                               ifelse(dad[,vars[i]]==3,0.875/7,
                                                      ifelse(dad[,vars[i]]==4,2.5/7,
                                                             ifelse(dad[,vars[i]]==5,4/7,NA)))))))
}

vars<-c("ff439","ff440","ff441","ff442")

for(i in 1:length(vars))
{
  dad[,vars[i]]<-with(dad,ifelse(is.na(dad[,vars[i]]),1,dad[,vars[i]]))
  dad[,vars[i]]<-with(dad,ifelse(dad[,vars[i]]==0,0,
                                 ifelse(dad[,vars[i]]==1,0,
                                        ifelse(dad[,vars[i]]==2,0.5/7,
                                               ifelse(dad[,vars[i]]==3,1.5/7,
                                                      ifelse(dad[,vars[i]]==4,3.5/7,
                                                             ifelse(dad[,vars[i]]==5,1,NA)))))))
}

vars<-c("ff445","ff446","ff447")

for(i in 1:length(vars))
{
  dad[,vars[i]]<-with(dad,ifelse(is.na(dad[,vars[i]]),1,dad[,vars[i]]))
  dad[,vars[i]]<-with(dad,ifelse(dad[,vars[i]]==0,0,
                                 ifelse(dad[,vars[i]]==1,0,
                                        ifelse(dad[,vars[i]]==2,0.5/7,
                                               ifelse(dad[,vars[i]]==3,2.5/7,
                                                      ifelse(dad[,vars[i]]==4,6/7,
                                                             ifelse(dad[,vars[i]]==5,1,NA)))))))
}

# Serving sizes in Norway: https://www.matportalen.no/verktoy/the_norwegian_food_composition_table/weights_measures_and_portion_sizes_for_foods #

dad$juice_fruit<-with(dad,ifelse(ff419<=1,ff419,
                                 ifelse(ff419>1,1,NA)))
dad$juice_juice<-with(dad,ff419-juice_fruit)
dad$fruit_dad<-with(dad,(ff442*150)+(juice_fruit*250))
dad$veg_dad<-with(dad,(ff439*150)+(ff440*75)+(ff441*175))
dad$nutleg_dad<-with(dad,ff437*150)
dad$wholeg_dad<-with(dad,(ff404*40)+(ff405*40))
dad$refg_dad<-with(dad,(ff403*40)+(ff406*15)+(ff430*150))
dad$snack_proxy_dad<-with(dad,ff445+ff446+ff447)
dad$sbev_dad<-with(dad,(juice_juice*250)+(ff420*250)+(ff421*250)+(ff422*250)+(ff423*250))
dad$cof_tea_dad<-with(dad,(ff424*250)+(ff425*250)+(ff426*100)+(ff427*250))
dad$dairy_lf_dad<-with(dad,(ff418*250)+(ff408*15))
dad$dairy_fat_dad<-with(dad,(ff417*250)+(ff409*15))
dad$dairy_tot_dad<-with(dad,dairy_lf_dad+dairy_fat_dad)
dad$rpmeat_dad<-with(dad,(ff411*30)+(ff412*30)+(ff413*20)+(ff428*180)+(ff429*75)+(ff430*50)+(ff431*75)+(ff432*150))
dad$wmeat_dad<-with(dad,ff433*150)
dad$tmeat_dad<-with(dad,rpmeat_dad+wmeat_dad)
dad$eggs_dad<-with(dad,ff416*50)
dad$fish_dad<-with(dad,(ff414*25)+(ff434*200)+(ff435*150)+(ff436*150))

#quantile(dad$refg_dad,probs=seq(0,1,0.025),na.rm=TRUE)
#hist(dad$cof_tea_dad)
#sort(dad$cof_tea_dad,decreasing=TRUE)
dad$wholeg_dad<-with(dad,ifelse(wholeg_dad>600,NA,wholeg_dad))
dad$refg_dad<-with(dad,ifelse(refg_dad>600,NA,refg_dad))
dad$dairy_tot_dad<-with(dad,ifelse(dairy_tot_dad>1500,NA,dairy_tot_dad))
dad$dairy_lf_dad<-with(dad,ifelse(dairy_lf_dad>1500,NA,dairy_lf_dad))
dad$sbev_dad<-with(dad,ifelse(sbev_dad>2000,NA,sbev_dad))
dad$rpmeat_dad<-with(dad,ifelse(rpmeat_dad>250,NA,rpmeat_dad))
dad$fish_dad<-with(dad,ifelse(fish_dad>250,NA,fish_dad))
dad$cof_tea_dad<-with(dad,ifelse(cof_tea_dad>2000,NA,cof_tea_dad))

dad<-dad[!is.na(dad$fruit_dad) & !is.na(dad$veg_dad) & !is.na(dad$nutleg_dad) & !is.na(dad$dairy_lf_dad) & 
           !is.na(dad$wholeg_dad) & !is.na(dad$snack_proxy_dad) & !is.na(dad$sbev_dad) & !is.na(dad$rpmeat_dad),]

dad$dashq_fruit_dad<-cut(dad$fruit_dad,quantile(dad$fruit_dad,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
dad$dashq_veg_dad<-cut(dad$veg_dad,quantile(dad$veg_dad,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
dad$dashq_nutleg_dad<-cut(dad$nutleg_dad,unique(quantile(dad$nutleg_dad,probs=seq(0,1,0.2),na.rm=TRUE)),include.lowest=TRUE,labels=FALSE)
dad$dashq_nutleg_dad<-with(dad,ifelse(dashq_nutleg_dad==1,2.5,
                                      ifelse(dashq_nutleg_dad==2,5,NA)))
dad$dashq_dairy_lf_dad<-cut(dad$dairy_lf_dad,quantile(dad$dairy_lf_dad,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
dad$dashq_wholeg_dad<-cut(dad$wholeg_dad,quantile(dad$wholeg_dad,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
dad$dashq_snacks_dad<-cut(dad$snack_proxy_dad,unique(quantile(dad$snack_proxy_dad,probs=seq(0,1,0.2),na.rm=TRUE)),include.lowest=TRUE,labels=FALSE)
dad$dashq_snacks_dad<-with(dad,ifelse(dashq_snacks_dad==1,4.5,
                                      ifelse(dashq_snacks_dad==2,3,
                                             ifelse(dashq_snacks_dad==3,2,
                                                    ifelse(dashq_snacks_dad==4,1,NA)))))
dad$dashq_sbev_dad<-cut(dad$sbev_dad,unique(quantile(dad$sbev_dad,probs=seq(0,1,0.2),na.rm=TRUE)),include.lowest=TRUE,labels=FALSE)
dad$dashq_sbev_dad<-with(dad,ifelse(dashq_sbev_dad==1,4.5,
                                    ifelse(dashq_sbev_dad==2,3,
                                           ifelse(dashq_sbev_dad==3,2,
                                                  ifelse(dashq_sbev_dad==4,1,NA)))))
dad$dashq_rpmeat_dad<-cut(dad$rpmeat_dad,quantile(dad$rpmeat_dad,probs=seq(0,1,0.2),na.rm=TRUE),include.lowest=TRUE,labels=FALSE)
dad$dashq_rpmeat_dad<-with(dad,ifelse(dashq_rpmeat_dad==1,5,
                                      ifelse(dashq_rpmeat_dad==2,4,
                                             ifelse(dashq_rpmeat_dad==3,3,
                                                    ifelse(dashq_rpmeat_dad==4,2,
                                                           ifelse(dashq_rpmeat_dad==5,1,NA))))))
dad$dash_dad<-with(dad,dashq_fruit_dad+dashq_veg_dad+dashq_nutleg_dad+dashq_dairy_lf_dad
                   +dashq_wholeg_dad+dashq_snacks_dad+dashq_sbev_dad+dashq_rpmeat_dad)

dad<-dad[,c("preg_id_2535","fruit_dad","veg_dad","nutleg_dad","wholeg_dad","refg_dad",
            "dairy_lf_dad","dairy_tot_dad","dairy_fat_dad","cof_tea_dad","sbev_dad",
            "eggs_dad","fish_dad","rpmeat_dad","wmeat_dad","tmeat_dad","snack_proxy_dad","dash_dad","diet_dad")]

dat<-merge2(dat,food,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,dad,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
dat$diet_mom<-with(dat,ifelse(is.na(diet_mom),0,diet_mom))
dat$diet_dad<-with(dat,ifelse(is.na(diet_dad),0,diet_dad))
food<-NULL
dad<-NULL


### METABOLOMICS DATA ###
#########################

metab<-spss.get("N:/durable/RAW/Metabolomics_data/PDB2535_Analyseresultater.sav",
                use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(metab)<-tolower(names(metab))

metab$rolle<-with(metab,ifelse(rolle=="SU2PT_MOTHER        ",1,
                               ifelse(rolle=="SU2PT_FATHER        ",2,NA)))
for(i in 1:length(metab))
{
  metab[,i]<-as.numeric(metab[,i])
}

metab_mom<-subset2(metab,"metab$rolle==1")
names(metab_mom)<-paste(names(metab_mom),"_mom",sep="")
metab_mom<-rename.vars(metab_mom,from=c("preg_id_2535_mom"),to=c("preg_id_2535"))
metab_mom$rolle_mom<-NULL

metab_dad<-subset2(metab,"metab$rolle==2")
names(metab_dad)<-paste(names(metab_dad),"_dad",sep="")
metab_dad<-rename.vars(metab_dad,from=c("preg_id_2535_dad"),to=c("preg_id_2535"))
metab_dad$rolle_dad<-NULL

metab<-merge2(metab_mom,metab_dad,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
metab_mom<-NULL
metab_dad<-NULL


# Selection bias #

metab$metabolomics<-1
dat_metab<-metab[,c("preg_id_2535","metabolomics")]
dat<-merge2(dat,dat_metab,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
dat$metabolomics<-with(dat,ifelse(is.na(metabolomics),0,metabolomics))
dat_metab<-NULL
save(dat,file="./Data/MoBa_diet.RData")

metab$metabolomics<-NULL
dat<-merge2(metab,dat,by.id=c("preg_id_2535"),all.x=TRUE,sort=FALSE)
save(dat,file="./Data/MoBa_diet_metab.RData")


### FLOW CHART ###
##################

load("./Data/MoBa_diet.RData")

# Total number of pregnancies with consent ok (n = 112052) #
length(which(!is.na(keys$preg_id_2535) & !duplicated(keys$preg_id_2535)))

# Total unique mothers with consent ok and infertility information (n = 68432) #
length(which(!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$preg_plan==1))

# Total unique mothers with consent ok, infertility and diet information (n = 56869) #
length(which(!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$diet_mom==1 & dat$preg_plan==1))

# Total unique mothers with consent ok, infertility and metabolomics information (n = 4066) #
length(which(!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$metabolomics==1 & dat$preg_plan==1))

# Total unique mothers with consent ok, infertility, diet, and metabolome information (n = 3880) #
length(which(!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$diet_mom==1 & dat$metabolomics==1 & dat$preg_plan==1))

# Total unique mothers with consent ok, infertility, diet, and genotype information (n = 53591) #
length(which(!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$diet_mom==1 & !is.na(dat$sentrixid_mom) & dat$preg_plan==1))

# Total unique mothers with consent ok, infertility, diet, genotype and metabolome information (n = 3868) #
length(which(!duplicated(dat$m_id_2535) & !is.na(dat$subf_mom) & dat$diet_mom==1 & !is.na(dat$sentrixid_mom) & dat$metabolomics==1 & dat$preg_plan==1))

# Total unique fathers with consent ok and infertility information (n = 56643) #
length(which(!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$preg_plan==1))

# Total unique fathers with consent ok, infertility and diet information (n = 26448) #
length(which(!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$diet_dad==1 & dat$preg_plan==1))

# Total unique fathers with consent ok, infertility and metabolomics information (n = 4066) #
length(which(!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$metabolomics==1 & dat$preg_plan==1))

# Total unique fathers with consent ok, infertility, diet, and metabolome information (n = 2011) #
length(which(!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$diet_dad==1 & dat$metabolomics==1 & dat$preg_plan==1))

# Total unique fathers with consent ok, infertility, diet, and genotype information (n = 21247) #
length(which(!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$diet_dad==1 & !is.na(dat$sentrixid_dad) & dat$preg_plan==1))

# Total unique fathers with consent ok, infertility, diet, genotype and metabolome information (n = 2005) #
length(which(!duplicated(dat$f_id_2535) & !is.na(dat$subf_dad) & dat$diet_dad==1 & !is.na(dat$sentrixid_dad) & dat$metabolomics==1 & dat$preg_plan==1))

