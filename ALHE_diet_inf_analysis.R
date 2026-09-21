rm(list=ls())

library(compareGroups)
library(dplyr)
library(gam)
library(ggplot2)
library(Hmisc)
library(lmtest)
library(miceadds)
library(paran)
library(splines)

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
                      ifelse(abs(x)<0.01,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),
                             ifelse(abs(x)<0.1,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),
                                    ifelse(abs(x)<1,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),guapa(x))))))
  return(pval)
}

mean_values <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  n <- sum(!is.na(x))
  se <- sqrt(var(x, na.rm=TRUE)/n)
  z <- qnorm(1-0.05/2)
  media <- mean(x, na.rm=TRUE)
  ic95a <- media - (z*se)
  ic95b <- media + (z*se)
  ic_ok <- c(media, ic95a, ic95b)
  return(ic_ok)
}

prev_values_complete <- function(x, category = 1) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) {
    return(c(
      cases = NA_real_,
      n = 0,
      prev = NA_real_,
      lo = NA_real_,
      hi = NA_real_))}
  k <- sum(x == category)
  ci <- prop.test(k, n, correct = FALSE)$conf.int
  c(cases = k,
    n = n,
    prev = k / n,
    lo = ci[1],
    hi = ci[2])}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x, na.rm=TRUE)
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

risk_se_ic_guapa4 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-exp(x)
  ic95a<-exp(x-(z*y))
  ic95b<-exp(x+(z*y))
  ic_ok<-c(hr,ic95a,ic95b)
  return(ic_ok)
}


header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)
se <- function(x) sqrt(var(x) / length(x))

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

options(scipen=999)

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

# 5 components retained for mother's diet #

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


# MOTHERS #

# 1. Total unique mothers #

datx<-dat[!is.na(dat$m_id_2535),]
datx<-datx[!duplicated(datx$m_id_2535),]

xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_mom1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_mom1)<-c("Total unique mothers")


# 2. Unique mothers - planned pregnancy #

datx<-dat[!is.na(dat$m_id_2535) & dat$preg_plan==1,]
datx<-datx[!duplicated(datx$m_id_2535),]

xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_mom2<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_mom2)<-c("Planned pregnancy")


# 3. Unique mothers - planned pregnancy + infertility information #

datx<-dat[!is.na(dat$m_id_2535) & dat$preg_plan==1 & !is.na(dat$subf_mom),]
datx<-datx[!duplicated(datx$m_id_2535),]

xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_mom3<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_mom3)<-c("Planned + infertility")


# 4. Unique mothers - planned pregnancy + infertility + own diet #

datx<-dat[!is.na(dat$m_id_2535) & dat$preg_plan==1 & !is.na(dat$subf_mom) & dat$diet_mom==1,]
datx<-datx[!duplicated(datx$m_id_2535),]

xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_mom4<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_mom4)<-c("Planned + infertility + own diet")


# 5. Unique mothers - planned pregnancy + infertility + own diet + diet in both parents #

datx<-dat[!is.na(dat$m_id_2535) & dat$preg_plan==1 & !is.na(dat$subf_mom) & dat$diet_mom==1 & dat$diet_dad==1,]
datx<-datx[!duplicated(datx$m_id_2535),]

xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_mom5<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_mom5)<-c("Planned + infertility + own diet + both diets")


# Combine mothers #

tab_mothers<-cbind(tab_mom1,tab_mom2,tab_mom3,tab_mom4,tab_mom5)

write.table(tab_mothers,file="./Outputs/Descriptive/descr_subgroups_mothers.csv",sep=";",col.names=NA)



# FATHERS #

# 1. Total unique fathers #

datx<-dat[!is.na(dat$f_id_2535),]
datx<-datx[!duplicated(datx$f_id_2535),]

xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_dad1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_dad1)<-c("Total unique fathers")


# 2. Unique fathers - planned pregnancy #

datx<-dat[!is.na(dat$f_id_2535) & dat$preg_plan==1,]
datx<-datx[!duplicated(datx$f_id_2535),]

xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_dad2<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_dad2)<-c("Planned pregnancy")


# 3. Unique fathers - planned pregnancy + infertility information #

datx<-dat[!is.na(dat$f_id_2535) & dat$preg_plan==1 & !is.na(dat$subf_dad),]
datx<-datx[!duplicated(datx$f_id_2535),]

xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_dad3<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_dad3)<-c("Planned + infertility")


# 4. Unique fathers - planned pregnancy + infertility + own diet #

datx<-dat[!is.na(dat$f_id_2535) & dat$preg_plan==1 & !is.na(dat$subf_dad) & dat$diet_dad==1,]
datx<-datx[!duplicated(datx$f_id_2535),]

xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_dad4<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_dad4)<-c("Planned + infertility + own diet")


# 5. Unique fathers - planned pregnancy + infertility + own diet + diet in both parents #

datx<-dat[!is.na(dat$f_id_2535) & dat$preg_plan==1 & !is.na(dat$subf_dad) & dat$diet_dad==1 & dat$diet_mom==1,]
datx<-datx[!duplicated(datx$f_id_2535),]

xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx,method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE,show.p.overall=FALSE,show.p.trend=FALSE,hide.no=0)

tab_dad5<-as.data.frame(cbind(all$descr[,1]))
colnames(tab_dad5)<-c("Planned + infertility + own diet + both diets")


# Combine fathers #

tab_fathers<-cbind(tab_dad1,tab_dad2,tab_dad3,tab_dad4,tab_dad5)

write.table(tab_fathers,file="./Outputs/Descriptive/descr_subgroups_fathers.csv",sep=";",col.names=NA)


### SELECTION BIAS - DIET vs NO DIET ###
########################################


# MOTHERS #

datx<-dat[!is.na(dat$subf_mom),]

datx$diet_mom<-with(datx,ifelse(preg_plan==1,diet_mom,
                                ifelse(preg_plan==0,0,diet_mom)))

datx<-datx[!duplicated(datx$m_id_2535),]

xxx<-datx[,c("agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom","subf_mom","art_mom","diet_mom")]

all<-NULL
all<-createTable(compareGroups(diet_mom~.,
                               xxx,method=c("bmimax_mom"=2,"smkinitmax_mom"=3,"subf_mom"=3,"art_mom"=3)),
                 show.n=TRUE,show.p.overall=TRUE,show.p.trend=FALSE,hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Non-included","Included","P-value","N")
write.table(tab1,file="./Outputs/Descriptive/selectionbias_diet_mothers.csv",sep=";",col.names=NA)


# FATHERS #

datx<-dat[!is.na(dat$subf_dad),]

datx$diet_dad<-with(datx,ifelse(preg_plan==1,diet_dad,
                                ifelse(preg_plan==0,0,diet_dad)))

datx<-datx[!duplicated(datx$f_id_2535),]

xxx<-datx[,c("agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad","subf_dad","art_dad","diet_dad")]

all<-NULL
all<-createTable(compareGroups(diet_dad~.,
                               xxx,method=c("bmimax_dad"=2,"smkinitmax_dad"=3,"subf_dad"=3,"art_dad"=3)),
                 show.n=TRUE,show.p.overall=TRUE,show.p.trend=FALSE,hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Non-included","Included","P-value","N")
write.table(tab1,file="./Outputs/Descriptive/selectionbias_diet_fathers.csv",sep=";",col.names=NA)


### LOGISTIC REGRESSION ###
###########################

load("./Data/MoBa_diet.RData")

vars00<-c("dash_mom","dash_dad",
          "fruit_mom","fruit_dad","veg_mom","veg_dad","nutleg_mom","nutleg_dad",
          "wholeg_mom","wholeg_dad","sbev_mom","sbev_dad","dairy_lf_mom","dairy_lf_dad",
          "fish_mom","fish_dad","rpmeat_mom","rpmeat_dad","snack_proxy_mom","snack_proxy_dad")
vars01<-c("dash_mom_z","dash_dad_z",
          "fruit_mom_z","fruit_dad_z","veg_mom_z","veg_dad_z","nutleg_mom_z","nutleg_dad_z",
          "wholeg_mom_z","wholeg_dad_z","sbev_mom_z","sbev_dad_z","dairy_lf_mom_z","dairy_lf_dad_z",
          "fish_mom_z","fish_dad_z","rpmeat_mom_z","rpmeat_dad_z","snack_proxy_mom_z","snack_proxy_dad_z")
vars02<-c("agedelivery_mom","agedelivery_dad",
          "agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_mom",
          "agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad",
          "agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad","agedelivery_mom","agedelivery_dad")
vars03<-c("eduyears_mom","eduyears_dad",
          "eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_mom",
          "eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad",
          "eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad","eduyears_mom","eduyears_dad")
vars04<-c("bmi_mom","bmi_dad",
          "bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_mom",
          "bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_dad",
          "bmi_mom","bmi_dad","bmi_mom","bmi_dad","bmi_mom","bmi_dad")
vars05<-c("smkinit_mom","smkinit_dad",
          "smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_mom",
          "smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad",
          "smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad","smkinit_mom","smkinit_dad")
vars06<-c("m_id_2535","f_id_2535",
          "m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","m_id_2535",
          "m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","f_id_2535",
          "m_id_2535","f_id_2535","m_id_2535","f_id_2535","m_id_2535","f_id_2535")
vars07<-c("diet_mom","diet_dad",
          "diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_mom",
          "diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_dad",
          "diet_mom","diet_dad","diet_mom","diet_dad","diet_mom","diet_dad")
vars08<-c("DASH diet score, women","DASH diet score, male partners",
          "Fruit intake (g/d), women","Fruit intake (g/d), male partners",
          "Vegetable intake (g/d), women","Vegetable intake (g/d), male partners",
          "Nut and legume intake (g/d), women","Nut and legume intake (g/d), male partners",
          "Whole grain intake (g/d), women","Whole grain intake (g/d), male partners",
          "Sweetened beverage intake (mL/d), women","Sweetened beverage intake (mL/d), male partners",
          "Low-fat dairy intake (g/d), women","Low-fat dairy intake (g/d), male partners",
          "Fish intake (g/d), women","Fish intake (g/d), male partners",
          "Red and processed meat intake (g/d), women","Red and processed meat intake (g/d), male partners",
          "Proxy of snack intake, women","Proxy of snack intake, male partners")

# Spline display settings
spline_ylim<-c(0.5,2)
serving_divisor<-c(fruit=100,veg=100,nutleg=30,wholeg=40,dairy_lf=200,fish=100,sbev=250,rpmeat=100)
xmax_serv<-c(fruit=5,veg=4,nutleg=3,wholeg=10,dairy_lf=5,fish=1.5,sbev=6,rpmeat=2)
spline_groups<-names(serving_divisor)

tab<-NULL

for(i in 1:length(vars01))
{
  datx<-dat[dat[,vars07[i]]==1 & dat$preg_plan==1 & !is.na(dat[,vars06[i]]),]
  
  # Standardize exposure in the complete sex-specific dietary sample
  exp_mean<-mean(datx[,vars00[i]],na.rm=TRUE)
  exp_sd<-sd(datx[,vars00[i]],na.rm=TRUE)
  datx[,vars01[i]]<-(datx[,vars00[i]]-exp_mean)/exp_sd
  sample<-length(which(!is.na(datx[,vars01[i]])))
  
  # MODEL 01: main sample, unadjusted
  mod01<-miceadds::glm.cluster(formula=as.factor(subf)~datx[,vars01[i]],
                               data=datx,cluster=datx[,vars06[i]],family="binomial")
  coef01<-risk_se_ic_guapa2(as.numeric(summary(mod01)[2,1]),as.numeric(summary(mod01)[2,2]))
  pval01<-pval_guapa(as.numeric(summary(mod01)[2,4]))
  or01<-exp(as.numeric(summary(mod01)[2,1]))
  orlo01<-exp(as.numeric(summary(mod01)[2,1])-(z*as.numeric(summary(mod01)[2,2])))
  orhi01<-exp(as.numeric(summary(mod01)[2,1])+(z*as.numeric(summary(mod01)[2,2])))
  
  # MODEL 02: main sample, adjusted
  mod02<-miceadds::glm.cluster(formula=as.factor(subf)~datx[,vars01[i]]
                               +datx[,vars02[i]]+datx[,vars03[i]]+datx[,vars04[i]]+datx[,vars05[i]]+parity,
                               data=datx,cluster=datx[,vars06[i]],family="binomial")
  coef02<-risk_se_ic_guapa2(as.numeric(summary(mod02)[2,1]),as.numeric(summary(mod02)[2,2]))
  pval02<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  or02<-exp(as.numeric(summary(mod02)[2,1]))
  orlo02<-exp(as.numeric(summary(mod02)[2,1])-(z*as.numeric(summary(mod02)[2,2])))
  orhi02<-exp(as.numeric(summary(mod02)[2,1])+(z*as.numeric(summary(mod02)[2,2])))
  
  # Identify partner diet variables
  is_mom<-grepl("_mom$",vars00[i])
  partner_diet<-if(is_mom) "diet_dad" else "diet_mom"
  partner_dash<-if(is_mom) "dash_dad" else "dash_mom"
  
  # Subset with dietary data in the partner
  # Exposure z-score is inherited from datx: same SD as models 01-02
  datp<-datx[datx[,partner_diet]==1 & !is.na(datx[,partner_dash]),]
  sample_partner<-length(which(!is.na(datp[,vars01[i]])))
  
  # MODEL 03: partner-diet subset, unadjusted
  mod03<-miceadds::glm.cluster(formula=as.factor(subf)~datp[,vars01[i]],
                               data=datp,cluster=datp[,vars06[i]],family="binomial")
  coef03<-risk_se_ic_guapa2(as.numeric(summary(mod03)[2,1]),as.numeric(summary(mod03)[2,2]))
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  or03<-exp(as.numeric(summary(mod03)[2,1]))
  orlo03<-exp(as.numeric(summary(mod03)[2,1])-(z*as.numeric(summary(mod03)[2,2])))
  orhi03<-exp(as.numeric(summary(mod03)[2,1])+(z*as.numeric(summary(mod03)[2,2])))
  
  # MODEL 04: partner-diet subset, adjusted as model 02
  mod04<-miceadds::glm.cluster(formula=as.factor(subf)~datp[,vars01[i]]
                               +datp[,vars02[i]]+datp[,vars03[i]]+datp[,vars04[i]]+datp[,vars05[i]]+parity,
                               data=datp,cluster=datp[,vars06[i]],family="binomial")
  coef04<-risk_se_ic_guapa2(as.numeric(summary(mod04)[2,1]),as.numeric(summary(mod04)[2,2]))
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  or04<-exp(as.numeric(summary(mod04)[2,1]))
  orlo04<-exp(as.numeric(summary(mod04)[2,1])-(z*as.numeric(summary(mod04)[2,2])))
  orhi04<-exp(as.numeric(summary(mod04)[2,1])+(z*as.numeric(summary(mod04)[2,2])))
  
  # MODEL 05: partner-diet subset, adjusted as model 02 + partner DASH
  mod05<-miceadds::glm.cluster(formula=as.factor(subf)~datp[,vars01[i]]
                               +datp[,vars02[i]]+datp[,vars03[i]]+datp[,vars04[i]]+datp[,vars05[i]]+parity
                               +datp[,partner_dash],
                               data=datp,cluster=datp[,vars06[i]],family="binomial")
  coef05<-risk_se_ic_guapa2(as.numeric(summary(mod05)[2,1]),as.numeric(summary(mod05)[2,2]))
  pval05<-pval_guapa(as.numeric(summary(mod05)[2,4]))
  or05<-exp(as.numeric(summary(mod05)[2,1]))
  orlo05<-exp(as.numeric(summary(mod05)[2,1])-(z*as.numeric(summary(mod05)[2,2])))
  orhi05<-exp(as.numeric(summary(mod05)[2,1])+(z*as.numeric(summary(mod05)[2,2])))
  
  # MAIN SPLINE
  aaa<-datx[,vars00[i]]
  basevar<-sub("_(mom|dad)$","",vars00[i])
  is_snack<-basevar=="snack_proxy"
  is_serv<-basevar %in% spline_groups
  aaa_spline<-if(is_snack) datx[,vars01[i]] else if(is_serv) aaa/serving_divisor[basevar] else aaa
  
  mod_lin<-glm(formula=as.factor(subf)~aaa_spline
               +datx[,vars02[i]]+datx[,vars03[i]]+datx[,vars04[i]]+datx[,vars05[i]]+parity,
               data=datx,family="binomial")
  mod_nlin<-gam(formula=as.factor(subf)~bs(aaa_spline,df=4)
                +datx[,vars02[i]]+datx[,vars03[i]]+datx[,vars04[i]]+datx[,vars05[i]]+parity,
                data=datx,family="binomial")
  
  p_nonlin_lrtest<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  tab<-rbind(tab,cbind(sample,
                       coef01,pval01,or01,orlo01,orhi01,
                       coef02,pval02,or02,orlo02,orhi02,
                       sample_partner,
                       coef03,pval03,or03,orlo03,orhi03,
                       coef04,pval04,or04,orlo04,orhi04,
                       coef05,pval05,or05,orlo05,orhi05,
                       p_nonlin_lrtest))
  
  ptemp<-termplot(mod_nlin,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp[[1]]
  value<-if(is_snack) closest(temp$x,0)[1] else closest(temp$x,min(aaa_spline,na.rm=TRUE))[1]
  center<-with(temp,y[x==value])
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  ci<-exp(ytemp-center)
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  
  if(is_serv) plot.data<-plot.data[plot.data$x>=0 & plot.data$x<=xmax_serv[basevar],]
  plot.data<-plot.data[plot.data$lci>=0.1 & plot.data$uci<=3,]
  
  name<-paste("./Outputs/Results/",vars00[i],".jpg",sep="")
  labely<-"Infertility risk (odds ratio, 95% CI)"
  labelx<-vars08[i]
  
  if(is_serv) labelx<-paste0(c(fruit="Fruit",veg="Vegetable",nutleg="Nut and legume",wholeg="Whole grain",
                               dairy_lf="Low-fat dairy",fish="Fish",sbev="Sweetened beverage",
                               rpmeat="Red and processed meat")[basevar],
                             " intake (servings/day), ",ifelse(is_mom,"women","male partners"))
  
  if(is_snack) labelx<-paste0("Fast-food/snack proxy (z-score), ",ifelse(is_mom,"women","male partners"))
  
  p_lin2<-pval_guapa(as.numeric(summary(mod_lin)$coefficients[2,4]))
  p_nonlin2<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  p_lin2<-ifelse(p_lin2=="<0.001"," < 0.001",
                 ifelse(p_lin2=="<0.00001"," < 0.00001",paste(" = ",p_lin2,sep="")))
  p_nonlin2<-ifelse(p_nonlin2=="<0.001"," < 0.001",
                    ifelse(p_nonlin2=="<0.00001"," < 0.00001",paste(" = ",p_nonlin2,sep="")))
  leg<-paste("p-value for linearity",p_lin2,"\np-value for non-linearity",p_nonlin2,sep="")
  
  figure<-ggplot(plot.data,aes(x=x,y=yest))+
    geom_ribbon(aes(ymin=lci,ymax=uci),alpha=0.25,fill="Black")+
    geom_line(color="Black")+
    geom_hline(yintercept=1,linetype=2)+
    theme_bw()+
    scale_x_continuous(expand=c(0,0))+
    labs(x=labelx,y=labely)+
    annotate("text",x=max(plot.data$x,na.rm=TRUE)*0.98,
             y=if(is_serv || is_snack) spline_ylim[2] else max(plot.data$uci,na.rm=TRUE),
             label=leg,vjust=1,hjust=1,size=6)+
    theme(axis.title.x=element_text(vjust=0.5,size=18,face="bold"),
          axis.title.y=element_text(vjust=0.5,size=18,face="bold"),
          axis.text.x=element_text(size=16,colour="black"),
          axis.text.y=element_text(size=16,colour="black"),
          axis.ticks.x=element_line(colour="black"),
          axis.ticks.y=element_line(colour="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y.right=element_blank(),
          axis.ticks.y.right=element_blank())
  
  if(is_serv || is_snack) figure<-figure+coord_cartesian(ylim=spline_ylim)
  if(is_snack) figure<-figure+geom_vline(xintercept=0,linetype=3)
  
  ggsave(filename=name,plot=figure,dpi=1200)
  figure
  
  
  # PARTNER-DIET SPLINE, ADJUSTED AS MODEL 02 + PARTNER DASH
  aaa_partner<-datp[,vars00[i]]
  
  # For snack, retain the z-score defined in the full sex-specific dietary sample
  aaa_spline_partner<-if(is_snack) datp[,vars01[i]] else if(is_serv) aaa_partner/serving_divisor[basevar] else aaa_partner
  
  mod_lin_partner<-glm(formula=as.factor(subf)~aaa_spline_partner
                       +datp[,vars02[i]]+datp[,vars03[i]]+datp[,vars04[i]]+datp[,vars05[i]]+parity
                       +datp[,partner_dash],
                       data=datp,family="binomial")
  
  mod_nlin_partner<-gam(formula=as.factor(subf)~bs(aaa_spline_partner,df=4)
                        +datp[,vars02[i]]+datp[,vars03[i]]+datp[,vars04[i]]+datp[,vars05[i]]+parity
                        +datp[,partner_dash],
                        data=datp,family="binomial")
  
  p_nonlin_partner<-pval_guapa(lrtest(mod_lin_partner,mod_nlin_partner)[2,5])
  
  ptemp<-termplot(mod_nlin_partner,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp[[1]]
  value<-if(is_snack) closest(temp$x,0)[1] else closest(temp$x,min(aaa_spline_partner,na.rm=TRUE))[1]
  center<-with(temp,y[x==value])
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  ci<-exp(ytemp-center)
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  
  if(is_serv) plot.data<-plot.data[plot.data$x>=0 & plot.data$x<=xmax_serv[basevar],]
  plot.data<-plot.data[plot.data$lci>=0.1 & plot.data$uci<=3,]
  
  name<-paste("./Outputs/Results/",vars00[i],"_partnerDASH.jpg",sep="")
  
  p_lin_partner<-pval_guapa(as.numeric(summary(mod_lin_partner)$coefficients[2,4]))
  p_nonlin_partner2<-pval_guapa(lrtest(mod_lin_partner,mod_nlin_partner)[2,5])
  
  p_lin_partner<-ifelse(p_lin_partner=="<0.001"," < 0.001",
                        ifelse(p_lin_partner=="<0.00001"," < 0.00001",paste(" = ",p_lin_partner,sep="")))
  p_nonlin_partner2<-ifelse(p_nonlin_partner2=="<0.001"," < 0.001",
                            ifelse(p_nonlin_partner2=="<0.00001"," < 0.00001",paste(" = ",p_nonlin_partner2,sep="")))
  leg<-paste("p-value for linearity",p_lin_partner,
             "\np-value for non-linearity",p_nonlin_partner2,sep="")
  
  figure_partner<-ggplot(plot.data,aes(x=x,y=yest))+
    geom_ribbon(aes(ymin=lci,ymax=uci),alpha=0.25,fill="Black")+
    geom_line(color="Black")+
    geom_hline(yintercept=1,linetype=2)+
    theme_bw()+
    scale_x_continuous(expand=c(0,0))+
    labs(x=labelx,y=labely)+
    annotate("text",x=max(plot.data$x,na.rm=TRUE)*0.98,
             y=if(is_serv || is_snack) spline_ylim[2] else max(plot.data$uci,na.rm=TRUE),
             label=leg,vjust=1,hjust=1,size=6)+
    theme(axis.title.x=element_text(vjust=0.5,size=18,face="bold"),
          axis.title.y=element_text(vjust=0.5,size=18,face="bold"),
          axis.text.x=element_text(size=16,colour="black"),
          axis.text.y=element_text(size=16,colour="black"),
          axis.ticks.x=element_line(colour="black"),
          axis.ticks.y=element_line(colour="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y.right=element_blank(),
          axis.ticks.y.right=element_blank())
  
  if(is_serv || is_snack) figure_partner<-figure_partner+coord_cartesian(ylim=spline_ylim)
  if(is_snack) figure_partner<-figure_partner+geom_vline(xintercept=0,linetype=3)
  
  ggsave(filename=name,plot=figure_partner,dpi=1200)
  figure_partner
}

rownames(tab)<-vars00
colnames(tab)<-c("N_main",
                 "coef01","pval01","or01","orlo01","orhi01",
                 "coef02","pval02","or02","orlo02","orhi02",
                 "N_partner",
                 "coef03","pval03","or03","orlo03","orhi03",
                 "coef04","pval04","or04","orlo04","orhi04",
                 "coef05","pval05","or05","orlo05","orhi05",
                 "p_nonlin_main")
write.table(tab,file="./Outputs/Results/results.csv",sep=";",col.names=NA)
