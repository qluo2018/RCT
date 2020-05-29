# RCT: 3-month oral Bumetanide for ASD 
# R code for Data Analysis 
# Dr. Yuan Dai, Dr. Lingli Zhang @ Xinhua Hospital 
# Dr. Qiang Luo @ Fudan University
# Date of release: 20th Feb 2020
# Date of modification: 29th May 2020 to address reviewers' comments
# Version: v2.0.0.1
# Email: qluo@fudan.edu.cn

library(dplyr)
library(tidyr)
library(ppcor)
library(psych)
library(psycho)
library(car)
library(lmerTest)
library(naniar)
library(xlsx)
library(ggpubr)
library(tibble)
library(gdata)
library(predictmeans)
library(DHARMa)
library(sjstats)
library(ggplot2)
library(Rmisc)
library(lmPerm)

#################################################
#################################################
## 1. modified intended to treat analysis      ##
#################################################
#################################################



#######################
#######################
## Step 1: read data ##
#######################
#######################

# ## 1.1 read clinical data: in the long format

bumedata <- read.csv("bumedata_behav_mrs_20200216.csv", header = TRUE)

##################################
# MANUALLY CHECKING IS NECESSARY #
##################################

bumedata <- bumedata[-which(bumedata$id==1104),]   #this subject did not take any drug

# 3) check whether the variables are numeric or factor. If it is factor, we need to know which levels it has; 
# we only allow subject names and sex to be factor in this case.
# BE CAREFUL, when converting factor to numeric. Replacing the NaN with NA in the level set before converting 
for (i in c(1:length(bumedata))){
  if(is.factor(bumedata[,i])){
    print(colnames(bumedata)[i])
    print(levels(bumedata[,i]))
  }
}

bumedata$group <- ifelse(bumedata$group=='1',1,0)
bumedata$stage<- ifelse(bumedata$stage=='3',1,0)
bumedata$sex <- ifelse(bumedata$sex_1male2female=='1',1,0)
bumedata$sex <- as.factor(bumedata$sex)
bumedata$cutoffDQ75IQ70 <- ifelse(bumedata$cutoffDQ75IQ70=='1',1,0)
bumedata$cutoffDQ75IQ70 <- as.factor(bumedata$cutoffDQ75IQ70)

# write the data after cleanning 
write.csv(bumedata, file = "bumedata_behav_mrs_0216_clean.csv")



####################################
####################################
## Step 2: Description Statistics ##
####################################
####################################


######################################################
# table 1. characteristics and background ############
######################################################
# check the distribution of the data
# check if any outlier 
# CAUTION: it is not easily to define an outlier by statistic only. We need a really good reason to exclude
# some one from the following analysis. For example, this data point has failed in some quality check. Otherwise,
# we might be asked to report the results using the full data and do sensitivity analysis (maybe a bootstrap)
bumedata$group4plot <- as.factor(ifelse(bumedata$group==1,'bumetanide','control'))
bumedata$stage4plot <- as.factor(ifelse(bumedata$stage==0,'baseline','3 months'))
bumedata$stage4plot = with(bumedata, relevel(stage4plot, "baseline"))
bumedata$group4plot = with(bumedata, relevel(group4plot, "control"))
levels(bumedata$group4plot)
levels(bumedata$stage4plot)

bumedata.baseline <- bumedata[which(bumedata$stage == 0),]
bumedata.3mont <- bumedata[which(bumedata$stage == 1),]


# 2.1 baseline comparison
table(bumedata.baseline$group, bumedata.baseline$sex)
prop.test(table(bumedata.baseline$group, bumedata.baseline$sex), correct = FALSE)

table(bumedata.baseline$group, bumedata.baseline$cutoffDQ75IQ70)
prop.test(table(bumedata.baseline$group, bumedata.baseline$cutoffDQ75IQ70), correct = FALSE)


print(colnames(bumedata)[17:43])
desmat <- matrix(rep(0,16*27), ncol = 16, nrow = 27)
rownames(desmat) <- colnames(bumedata.baseline)[c(17:43)]
colnames(desmat) <- c("n","min", "max", "mean", "sd",
                      "n","min", "max", "mean", "sd",
                      "t.df", "t.t", "t.p",
                      "kw.df", "kw.chi-squared", "kw.p")
for (i in 17:43){
  # descriptive table:: desmat
  description <- describeBy(bumedata.baseline[,i], group = bumedata.baseline$group, mat = TRUE, na.rm = TRUE )
  desmat[i-16,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")],
                                         description[2, c("n","min", "max", "mean", "sd")]))
  # baseline comparison: testmat
  ttest <- t.test(bumedata.baseline[,i]~group, data = bumedata.baseline)
  desmat[i-16,11] <- ttest$parameter # degree-of-freedom
  desmat[i-16,12] <- ttest$statistic # t-statistic
  desmat[i-16,13] <- ttest$p.value #p value
  kwtest <- kruskal.test(bumedata.baseline[,i]~group, data = bumedata.baseline)
  desmat[i-16,14] <- kwtest$parameter # degree-of-freedom
  desmat[i-16,15] <- kwtest$statistic # t-statistic
  desmat[i-16,16] <- kwtest$p.value #p value
}

print(desmat)
write.xlsx(desmat, file="output_rct_behav119.xls", sheetName = "baselinebehaviour")


# 2.2 follow-up
desmat <- matrix(rep(0,16*17), ncol = 16, nrow = 17)
rownames(desmat) <- colnames(bumedata.3mont)[c(27:43)]
colnames(desmat) <- c("n","min", "max", "mean", "sd",
                      "n","min", "max", "mean", "sd",
                      "t.df", "t.t", "t.p",
                      "kw.df", "kw.chi-squared", "kw.p")
for (i in 27:43){
  # descriptive table:: desmat
  description <- describeBy(bumedata.3mont[,i], group = bumedata.3mont$group, mat = TRUE, na.rm = TRUE )
  desmat[i-26,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")],
                                         description[2, c("n","min", "max", "mean", "sd")]))
  # follow-up comparison: testmat
  ttest <- t.test(bumedata.3mont[,i]~group, data = bumedata.3mont)
  desmat[i-26,11] <- ttest$parameter # degree-of-freedom
  desmat[i-26,12] <- ttest$statistic # t-statistic
  desmat[i-26,13] <- ttest$p.value #p value
  kwtest <- kruskal.test(bumedata.3mont[,i]~group, data = bumedata.3mont)
  desmat[i-26,14] <- kwtest$parameter # degree-of-freedom
  desmat[i-26,15] <- kwtest$statistic # t-statistic
  desmat[i-26,16] <- kwtest$p.value #p value
}
print(desmat)
write.xlsx(desmat, file="output_rct_behav119.xls", sheetName = "month3behaviour", append = T)


# ######################################################
# ######################################################
# ## Step 3: test the drug effect ##
# ######################################################
# ######################################################
# 
# # for those with the shaptio or the levent test being significant, we need to do permutation 
# # LMM: liner mixed effect model with nuisance variables
stats3 <- data.frame(t = rep(0,17), t.p = rep(1,17), f=rep(0,17), f.p=rep(1,17), f.df=rep(1,17), eta.sq.partial=rep(1,17))
rownames(stats3) <- colnames(bumedata[27:43])

print(colnames(bumedata)[27:43])
for (i in 27:43){
  score <- bumedata[,i]
  fit <- lmer( score ~  stage + group + group*stage  + (1|id), data=bumedata)
  stats3$shapitestp[i-26] <- shapiro.test(residuals(fit))$p.value
  stats3$leventestp[i-26] <- leveneTest(score ~ as.factor(group*stage), data=bumedata)[1,3]
  stats3$t[i-26] <- summary(fit)$coefficients["stage:group","t value"]
  stats3$t.p[i-26] <- summary(fit)$coefficients["stage:group","Pr(>|t|)"]
  stats3$t.df[i-26] <- summary(fit)$coefficients["stage:group","df"]
  ftests <- anova(fit, type=1, ddf="Kenward-Roger")
  stats3$f[i-26] <- ftests["stage:group","F value"]
  stats3$f.p[i-26] <- ftests["stage:group","Pr(>F)"]
  stats3$f.df[i-26] <- ftests["stage:group","DenDF"]
  stats3$eta.sq.partial[i-26] <- eta_sq(fit, partial = TRUE)[3,2] 
}

print(stats3)
write.xlsx(stats3, file="output_rct_behav119.xls", sheetName = "drugeffectonbehaviourModel3", append = T)


# # permutation test for subscales
set.seed(300)
currentdata <- bumedata
perm.p.cars <- matrix(rep(0,15), ncol=1, nrow=15)
# suggested by our open-label RCT findings, we focused on the following 6 
rownames(perm.p.cars) <- colnames(currentdata)[c(29,31,32,33,35,41)]

for (i in rownames(perm.p.cars)){
  score <- (currentdata[[i]])
  # linear model comparison
  fit1 <- lmer(score ~ stage * group + stage + group + (1|id), data = currentdata)
  fit0 <- lmer(score ~ stage + group + (1|id), data = currentdata)
  perm.p.cars[[i,1]] <- permlmer(fit0, fit1, perms = 3000, ncore = 8, plot = FALSE)$"Perm-p"[2]
}

print(perm.p.cars)

perm.p.cars.fdr  <- p.adjust(perm.p.cars, method = "fdr") 
print(perm.p.cars.fdr)
write.xlsx(cbind(perm.p.cars, perm.p.cars.fdr), file="output_rct_behav119.xls",
           sheetName = "drugbehModel3PermFDR", append = T)


# 3.4 compare CGI-I

# nonparametric test
stats.cgi <- matrix(rep(0, 3), ncol = 3, nrow = 1)
rownames(stats.cgi) <- colnames(bumedata.3mont)[44]
colnames(stats.cgi) <- c("degree of freedom", "Kruskal-Wallis chi-squared", "p-value")

for (i in c(44)){
  kstest <- kruskal.test(bumedata.3mont[,i]~group, data=bumedata.3mont)  # significant
  stats.cgi[i-43,c(1:3)] <- c(kstest$parameter, kstest$statistic, kstest$p.value)
}
print(stats.cgi)
write.xlsx(stats.cgi, file="output_rct_behav119.xls", sheetName = "month3CGI2", append = T)


##########################
##########################
## Step 4: MRS analysis ##
##########################
##########################



## 4.0 quality control
ROI.names = "ic"
NAAplusSD = paste0(ROI.names, '_NAA_NAAG_SD')
GABAplusSD = paste0(ROI.names, '_GABA_SD')
fwhm = paste0(ROI.names, '_fwhm')
snr = paste0(ROI.names, '_snr')

# some info of data quality
des.QC <- matrix(rep(0,8), ncol = 2*4, nrow = 1)
rownames(des.QC) <- ROI.names
colnames(des.QC) <- c("NAA.SD.mean", "NAA.SD.sd","GABA.SD.mean", "GABA.SD.sd","FW.mean", "FW.sd","SNR.mean", "SNR.sd")
des.QC[,c(1,2)] <- as.matrix(describe(bumedata[[NAAplusSD]]))[c(3,4)]
des.QC[,c(3,4)] <- as.matrix(describe(bumedata[[GABAplusSD]]))[c(3,4)]
des.QC[,c(5,6)] <- as.matrix(describe(bumedata[[fwhm]]))[c(3,4)]
des.QC[,c(7,8)] <- as.matrix(describe(bumedata[[snr]]))[c(3,4)]

##############################################################################
## qc threshold,NAA+NAAG %SD < 20 & GABA %SD < 20 &  FWHM < 0.05 & SNR > 15 ##
##############################################################################
bumedata$QC <- matrix(rep(0,dim(bumedata)[1]), nrow = dim(bumedata)[1])
bumedata$QC <- (bumedata[[NAAplusSD]] < 20 & bumedata[[GABAplusSD]] < 20 & bumedata[[fwhm]] < 0.05 & bumedata[[snr]] > 15) 
bumedata$QC <- ifelse(bumedata$QC=='TRUE',1,0)

# subjects ID 1071 and 1045, T1 images were missing at month 3
bumedata$QC[which(bumedata$id==1071 & which(bumedata$stage == 1))] = 0
bumedata$QC[which(bumedata$id==1045 & which(bumedata$stage == 1))] = 0
# subject ID 1012 IC voxel,wrong position at both time points
bumedata$QC[which(bumedata$id==1012)] = 0


# count Num, surviving qc
counts <- matrix(rep(0,3*1), nrow = 3, ncol = 1, dimnames = list(c("baseline", "3 month", 'paired'), ROI.names))
temp1 <- bumedata$QC[which(bumedata$stage==0)]
counts[1,] <- length(which(temp1 == 1))
temp2 <- bumedata$QC[which(bumedata$stage==1)]
counts[2,] <- length(which(temp2 == 1))

mrsqc <- data.frame(cbind(bumedata$id, bumedata$QC, bumedata$stage))
colnames(mrsqc) <- c('id', 'ic.QC','stage')
mrsqc_paired <- merge(mrsqc[which(mrsqc$stage==0),], mrsqc[which(mrsqc$stage==1),], by = "id")
mrsqc_paired$ic.QC <- ((mrsqc_paired$ic.QC.x == mrsqc_paired$ic.QC.y) & mrsqc_paired$ic.QC.x == 1)
bumedata$QC.paired <- matrix(rep(0,dim(bumedata)[1]), nrow = dim(bumedata)[1])
for (i in c(1:dim(bumedata)[1])){
  bumedata$QC.paired[i] <- mrsqc_paired$ic.QC[which(mrsqc_paired$id == bumedata$id[i])]
}

counts[3,] <- colSums(bumedata$QC.paired, na.rm = T)/2
print(counts)
write.csv(bumedata, file = "bumedata_behav_mrs_0509_qc.csv")
write.xlsx(counts, file="output_rct_mrs.xls", sheetName = "MRSafterQC", append = T)


# 4.1 data description 
currentdata = bumedata[bumedata$stage==0,]
ROI.names = c("ic")
MRS.trans = c("_naa_corr", "_gaba_corr", "_glx_corr","_gaba_glx_corr")
desmat <- matrix(rep(0,4*18), ncol = 18, nrow = 4)  
colnames(desmat) <- c("n","min", "max", "mean", "sd", "n2","min2", "max2", 
                      "mean2", "sd2", 
                      "shapi.p", "levene.p", "t", "t.p", "F", "F.p", 
                      "F.Df", "ks.p")
rownames(desmat) <- c(1:4)
tableCount <- matrix(rep(0,length(MRS.trans)*length(ROI.names)), ncol=length(MRS.trans), nrow=length(ROI.names))
dimnames(tableCount) <- list(ROI.names, MRS.trans)
testmat <- list()
bp <- list()
s <- 0
z <- 0
index = matrix(c(1,1), ncol = 1, nrow = 1)
rownames(index) <- ROI.names
t <- tableCount

for (j in MRS.trans){ 
  
  s <- s + 1
  
  NOI <- paste0(ROI.names,j)
  validity <- (currentdata$QC == 1) 
  tableCount[ROI.names,j] = sum(validity==TRUE, na.rm=T)
  score <- (currentdata[[NOI]][validity==1])
  # descriptive table:: desma
  description <- describeBy(score, group = currentdata$group[validity==1], mat = TRUE, na.rm = TRUE )
  rownames(desmat)[s] <- NOI
  desmat[s,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")], 
                                       description[2, c("n","min", "max", "mean", "sd")]))
  
  # baseline comparison: testmat
  
  # linear model comparison
  fit <- lm(score~group + age + sex + cutoffDQ75IQ70, data = currentdata[validity==1,])
  desmat[s,11] <- shapiro.test(residuals(fit))$p.value
  desmat[s,12]<- leveneTest(score ~ as.factor(group), data=currentdata[validity==1,])[1,3]
  desmat[s,13] <- summary(fit)$coefficients["group","t value"]
  t[ROI.names,j] <- desmat[s,13]
  desmat[s,14] <- summary(fit)$coefficients["group","Pr(>|t|)"]
  ftests <- anova(fit)
  desmat[s,15]<- ftests["group","F value"]
  desmat[s,16] <- ftests["group","Pr(>F)"]
  desmat[s,17] <- ftests["Residuals", "Df"]
  # non-parametric comparison
  testmat[[NOI]] <- kruskal.test(score~currentdata$group[validity==1])
  desmat[s,18] <- testmat[[NOI]]$p.value
  
  bp[[s]] <- ggplot(currentdata[validity==1,], aes_string(x = NOI)) + 
    geom_histogram() + labs(x=sub(".tissue.corrected", "", NOI))
  
}


ggarrange(plotlist = bp, 
          ncol = length(MRS.trans), nrow = length(ROI.names))

print(desmat)
write.xlsx(desmat, file="output_rct_mrs.xls", sheetName = "MRSbaseline", append = T)

######################################################################################### 
## AJP revision: add-in analysis about NAA################################################# 
## studying whether baseline NAA differed!!! bl-NAA, norm with (X-mean(X))./std(X), within each grp
######################################################################################### 
j = "_naa_corr"
NOI <- paste0(ROI.names,j)

currentdata = bumedata[bumedata$stage==0,]
m <- is.na(currentdata$QC)
currentdata <- currentdata[!m,]
validity <- (currentdata$QC == 1) 
data = currentdata[validity==1,]
NAAscore <- (data[[NOI]])

# baseline NAA norm
NAAnorm = (NAAscore-min(NAAscore))/(max(NAAscore)-min(NAAscore))
desmat <- matrix(rep(0,1*18), ncol = 18, nrow = 1)
colnames(desmat) <- c("n","min", "max", "mean", "sd", "n2","min2", "max2", 
                      "mean2", "sd2", 
                      "shapi.p", "levene.p", "t", "t.p", "F", "F.p", 
                      "F.Df", "ks.p")

# descriptive table:: desma
description <- describeBy(NAAnorm, group =currentdata$group[validity==1], mat = TRUE, na.rm = TRUE )
# rownames(desmat)[s] <- NOI
desmat[,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")], 
                                    description[2, c("n","min", "max", "mean", "sd")]))
# baseline comparison: testmat

# linear model comparison
fit <- lm(NAAnorm~group + age + sex + cutoffDQ75IQ70, data)
desmat[,11] <- shapiro.test(residuals(fit))$p.value
desmat[,12]<- leveneTest(NAAnorm ~ as.factor(group), data)[1,3]
desmat[,13] <- summary(fit)$coefficients["group","t value"]
t[ROI.names,j] <- desmat[,13]
desmat[,14] <- summary(fit)$coefficients["group","Pr(>|t|)"]
ftests <- anova(fit)
desmat[,15]<- ftests["group","F value"]
desmat[,16] <- ftests["group","Pr(>F)"]
desmat[,17] <- ftests["Residuals", "Df"]
# non-parametric comparison
testmat[[NOI]] <- kruskal.test(NAAnorm~data$group)
desmat[,18] <- testmat[[NOI]]$p.value
desmat
# figure
ggboxplot(data,x = 'group', y = "NAAnorm", file = 'group', bxp.errorbar = T,bxp.errorbar.width = 0.2, palette = 'npg',
                    add = 'point')


######################################################################################### 
## AJP revision: add-in analysis about NAA################################################# 
## studying whether bum influenced NAA!!! NAA norm, within each timepoint
######################################################################################### 
j = "_naa_corr"
NOI <- paste0(ROI.names,j)
bume <- bumedata
m <- is.na(bume$QC.paired)
bume <- bume[!m,]
validity <- (bume$QC.paired == 1) 
data = bume[validity==1,]
NAAscore <- (data[[NOI]])
# NAA norm, within each timepoint
NAAnorm = NAAscore
NAAnorm[data$stage==0] <- (NAAscore[data$stage==0]-min(NAAscore[data$stage==0]))/(max(NAAscore[data$stage==0])-min(NAAscore[data$stage==0]))
NAAnorm[data$stage==1] <- (NAAscore[data$stage==1]-min(NAAscore[data$stage==1]))/(max(NAAscore[data$stage==1])-min(NAAscore[data$stage==1]))


# LMM
desmat <- matrix(rep(0,1*9), ncol = 9, nrow = 1)
colnames(desmat) <- c("t","t.p","f","f.p", "f.df", "shapitestp", "leventestp", "eta.sq.partial", "perm.p")

fit <- lmer( NAAnorm ~  stage + group + group*stage + age + sex + cutoffDQ75IQ70 + (1|id), data)  
desmat[,1]<- summary(fit)$coefficients["stage:group","t value"]
desmat[,2]<- summary(fit)$coefficients["stage:group","Pr(>|t|)"]
# t.df[ROI.names,j] <- summary(fit)$coefficients["stage:group","df"]
ftests <- anova(fit, type=1, ddf="Kenward-Roger")
desmat[,3] <- ftests["stage:group","F value"]
desmat[,4] <- ftests["stage:group","Pr(>F)"]
desmat[,5] <- ftests["stage:group","DenDF"]
desmat[,6] <- shapiro.test(residuals(fit))$p.value
desmat[,7]<- leveneTest(NAAnorm ~ as.factor(group), data=bume[validity==1,])[1,3]
desmat[,8] <- eta_sq(fit, partial = TRUE)[6,2] 

# permutation test
set.seed(300)
fit1 <- lmer(NAAnorm  ~ stage + group + stage * group + age + sex + cutoffDQ75IQ70 + (1|id), data = bume[validity==1,])
fit0 <- lmer(NAAnorm  ~ stage + group + age + sex + cutoffDQ75IQ70 + (1|id), data = bume[validity==1,])
desmat[,9] <- permlmer(fit0, fit1, perms = 3000, ncore = 8, plot = FALSE)$"Perm-p"[2]
desmat

#############AJP revision#################################
#Diff of GABA change btw 2 grps
validity <- (bume$QC.paired == 1) 
data = bume[validity==1,]

desmat <- matrix(rep(0,17*1), ncol = 17, nrow = 1)
colnames(desmat) <- c("n","min", "max", "mean", "sd", "n2","min2", "max2", 
                      "mean2", "sd2", 
                      "shapi.p", "levene.p", "t", "t.p", "F", "F.p", 
                      "F.Df")

description <- describeBy((data$ic_gaba_corr[data$stage==0]-data$ic_gaba_corr[data$stage==1]), group = (data$group[data$stage==0]), mat = TRUE, na.rm = TRUE )
desmat[1,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")],
                                     description[2, c("n","min", "max", "mean", "sd")]))
# change comparison: testmat
ttest <- t.test((data$ic_gaba_corr[data$stage==1]-data$ic_gaba_corr[data$stage==0])~data$group[data$stage==0], data)
desmat[1,11] <- ttest$parameter # degree-of-freedom
desmat[1,12] <- ttest$statistic # t-statistic
desmat[1,13] <- ttest$p.value #p value

data2 = data[data$stage==0,]
fit <- lm((data$ic_gaba_corr[data$stage==1]-data$ic_gaba_corr[data$stage==0]) ~ group+ age+ sex + cutoffDQ75IQ70, data2)
desmat[1,11] <- shapiro.test(residuals(fit))$p.value
desmat[1,12]<- leveneTest((data$ic_gaba_corr[data$stage==1]-data$ic_gaba_corr[data$stage==0]) ~ as.factor(group),data2)[1,3]
desmat[1,13] <- summary(fit)$coefficients["group","t value"]
desmat[1,14] <- summary(fit)$coefficients["group","Pr(>|t|)"]
ftests <- anova(fit)
desmat[1,15]<- ftests["group","F value"]
desmat[1,16] <- ftests["group","Pr(>F)"]
desmat[1,17] <- ftests["Residuals", "Df"]
print(desmat)
#############Above, AJP revision#######################

# 3) LMM
MRS.trans = c("_naa_corr", "_gaba_corr", "_glx_corr","_gaba_glx_corr")
MRS.trans <- MRS.trans[-c(1,3,4)] # our previous open-label RCT suggests the changes in GABA, so we focused on GABA in the current analysis

#MRS.trans = c("_naa_corr", "_gaba_corr", "_glx_corr","_gaba_glx_corr")
#MRS.trans <- MRS.trans[-c(2:4)]  ##keep NAA

#MRS.trans = c("_naa_corr", "_gaba_corr", "_glx_corr","_gaba_glx_corr")
#MRS.trans <- MRS.trans[-c(1:2,4)]  ##keep glx

t <- matrix(rep(0, length(ROI.names) * length(MRS.trans)), nrow = length(ROI.names), ncol = length(MRS.trans))
dimnames(t) <- list(ROI.names, MRS.trans)
shapitestp <- t
leventestp <- t
eta.sq.partial <- t
f <- t
t.p <- t + 1
t.df <- t
f.p <- t.p
f.df <- t

for (j in MRS.trans){
  
  NOI <- paste0(ROI.names,j)
  print(NOI)
  
  validity <- (bumedata$QC.paired == 1 & !is.na(bumedata$QC.paired)) 
  
  
  dat <- bumedata[validity==TRUE,]
  score <- (dat[[NOI]])
  #baseline.score <- score
  #baseline.score[dat$stage==1] <- baseline.score[dat$stage==0]
  fit <- lmer( score ~  stage + group + group*stage + age + sex + cutoffDQ75IQ70 + (1|id), data=dat)  
  shapitestp[ROI.names,j] <- shapiro.test(residuals(fit))$p.value
  leventestp[ROI.names,j] <- leveneTest(score ~ as.factor(group), data=bumedata[validity==1,])[1,3]
  t[ROI.names,j] <- summary(fit)$coefficients["stage:group","t value"]
  t.p[ROI.names,j] <- summary(fit)$coefficients["stage:group","Pr(>|t|)"]
  t.df[ROI.names,j] <- summary(fit)$coefficients["stage:group","df"]
  ftests <- anova(fit, type=1, ddf="Kenward-Roger")
  f[ROI.names,j] <- ftests["stage:group","F value"]
  f.p[ROI.names,j] <- ftests["stage:group","Pr(>F)"]
  f.df[ROI.names,j] <- ftests["stage:group","DenDF"]
  eta.sq.partial[ROI.names,j] <- eta_sq(fit, partial = TRUE)[6,2] 
}


# permutation test
set.seed(300)
perm.p <- matrix(rep(0, length(ROI.names) * length(MRS.trans)), nrow = length(ROI.names), ncol = length(MRS.trans))
dimnames(perm.p) <- list(ROI.names, MRS.trans)
currentdata = bumedata

for (j in MRS.trans){
  NOI <- paste0(ROI.names,j)
  print(NOI)
  validity <- (currentdata$QC.paired == 1 & !is.na(currentdata$QC.paired)) #(abs(scale(currentdata[[NOI]])) < 3) 
  score <- (currentdata[[NOI]][validity==1])
  #baseline.score <- score
  #baseline.score[currentdata$stage[validity==1]==1] <- baseline.score[currentdata$stage[validity==1]==0]
  # linear model comparison
  fit1 <- lmer(score ~ stage + group + stage * group + age + sex + cutoffDQ75IQ70 + (1|id), data = currentdata[validity==1,])
  fit0 <- lmer(score ~ stage + group + age + sex + cutoffDQ75IQ70 + (1|id), data = currentdata[validity==1,])
  perm.p[ROI.names,j] <- permlmer(fit0, fit1, perms = 3000, ncore = 8, plot = FALSE)$"Perm-p"[2]
}

stats3.mrs <- data.frame(t = t, t.p = t.p, f=f, f.p=f.p, f.df = f.df,
                         shapitestp = shapitestp, leventestp = leventestp, eta.sq.partial = eta.sq.partial,
                         perm.p = perm.p)
colnames(stats3.mrs) <- c("t","t.p","f","f.p", "f.df", "shapitestp", "leventestp", "eta.sq.partial", "perm.p")
print(stats3.mrs)
write.xlsx(stats3.mrs, file="output_rct_mrs.xls", sheetName = "DEMRSM3baselineCov2", append = T)


#################
# Figure 2B
#################
ICGABA <- bumedata[which(bumedata$QC.paired==1), c("stage4plot","group4plot","id","ic_gaba_corr")]
colnames(ICGABA)[c(1,2)] <- c("stage","group")
#ICGABA <- summarySE(ICGABA, measurevar="ic_gaba_corr",groupvars=c("group","stage"))
rm(p)
p <- ggbarplot(ICGABA, x="group", y="ic_gaba_corr", fill = "stage", 
               add = "mean_se", add.params = list(color = "black", group="stage"),
               position=position_dodge(0.8)) + 
  stat_compare_means(aes(group = stage), paired = TRUE, label = "p.signif", label.y = 0.5)
p+labs(x="Group",y="Change of GABA/NAA ratio from baseline to month 3")
ggsave("Figure2B.pdf")
dev.off()


# 4.3 correlations between behaviour and MRS association

# 1) cross-sectional correlations at both baseline and follow-up 
stage = levels(bumedata$stage4plot)[c(1,2)]
group = levels(bumedata$group4plot)
beh = colnames(bumedata)[c(27,28,29,31,32,33,35,41)]
dat <- data.frame()
num.records <- length(stage) * length(group) * length(ROI.names) * length(MRS.trans)
lookupindex <- array(rep(0,num.records*length(beh)), 
                     dim = c(length(stage), length(group), length(ROI.names), length(MRS.trans), length(beh)), 
                     dimnames=list(stage, group, ROI.names, MRS.trans, beh))
stats.corr <- data.frame(r = rep(0,num.records*length(beh)), 
                         r.p = rep(1,num.records*length(beh)), 
                         n = rep(0,num.records*length(beh)))
significant <- matrix(rep(0,num.records*length(beh)), nrow = 1, ncol = num.records*length(beh))
s = 0

for (k in stage){
  for (g in group){
    for (j in MRS.trans){
      NOI <- paste0(ROI.names,j)
      validity <- bumedata$QC
      dat <- bumedata[validity==TRUE & bumedata$group4plot == g & bumedata$stage4plot == k, ]
      covariates <- cbind(bumedata$age, bumedata$sex, bumedata$cutoffDQ75IQ70)
      for (b in beh){
        print(s)
        s = s + 1
        lookupindex[k,g,ROI.names,j,b] <- s 
        rownames(stats.corr)[s] <- paste(k,g,ROI.names,j,b,sep = ",")
        datanow <- cbind(dat[[NOI]],dat[[b]], covariates)
        if (length(which(rowSums(is.na(datanow))>0)) > 0){
          datanow <- datanow[-which(rowSums(is.na(datanow))>0),]
        }
        r <- pcor.test(datanow[,1],datanow[,2],datanow[,c(3:5)], method = "pearson")
      }
      stats.corr$r.p[s] = r$p.value
      stats.corr$r[s] = r$estimate
      stats.corr$n[s] = r$n
      if (r$p.value < 0.05) {
        significant[s] <- 1
      }
    }
  }
}

print(stats.corr)
# make a table
write.xlsx(stats.corr, file="output_rct_behav_mrs.xls", sheetName = "crosssectionalbehcorr", append = T)


# 2) delta correlations between  baseline and follow-up 
set.seed(300)
ROI.names = c("ic")
dat <- data.frame()
lookupindex <- array(rep(0,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)), dim = c(length(group),length(ROI.names),length(MRS.trans),length(beh)), 
                     dimnames=list(group, ROI.names, MRS.trans, beh))
stats.corr <- data.frame(estimate = rep(1,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)), 
                         pvalue = rep(1,length(group)*length(ROI.names)*length(MRS.trans)*length(beh))) 
significant <- matrix(rep(0,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)), 
                      nrow = 1, ncol = length(group)*length(ROI.names)*length(MRS.trans)*length(beh))
s = 0

for (g in group){
  for (j in MRS.trans){
    NOI <- paste0(ROI.names,j)
    
    validity <- bumedata$QC.paired
    dat <- bumedata[which(validity==TRUE),]
    
    baseline.mrs <- dat[[NOI]][which(dat$stage4plot == "baseline" & dat$group4plot == g)]
    delta.mrs <- baseline.mrs - dat[[NOI]][which(dat$stage4plot == "3 months" & dat$group4plot == g)]
    covariates <- cbind(dat$age, dat$sex, dat$cutoffDQ75IQ70)[which(dat$stage4plot == "baseline" & dat$group4plot == g),]

    for  (b in beh){
      s  = s + 1
      lookupindex[g,ROI.names,j,b] <- s 
      rownames(stats.corr)[s] <- paste(g,ROI.names,j,b,sep = ".")
      baseline.beh <- dat[[b]][which(dat$stage4plot=="baseline" & dat$group4plot == g)]
      delta.beh <- baseline.beh - dat[[b]][which(dat$stage4plot == "3 months" & dat$group4plot == g)]
      
      datanow <- cbind(baseline.mrs, delta.mrs, baseline.beh, delta.beh, covariates)
      if (length(which(rowSums(is.na(datanow))>0)) > 0){
        datanow <- datanow[-which(rowSums(is.na(datanow))>0),]
      }
      colnames(datanow)[c(5:7)] <- c('age','sex','DIQ')
      datanow <- as.data.frame(datanow)

      fit <- lmp(delta.beh ~ delta.mrs + baseline.beh + baseline.mrs + age + sex + DIQ, data=datanow)
      fit2 <- lm(delta.beh ~ delta.mrs + baseline.beh + baseline.mrs + age + sex + DIQ, data=datanow)
      
      stats.corr$pvalue[s] = summary(fit)$coefficients['delta.mrs', 'Pr(Prob)']
      stats.corr$estimate[s] = summary(fit)$coefficients['delta.mrs', 'Estimate'] 
      stats.corr$tvalue[s]=summary(fit2)$coefficients['delta.mrs', 't value'] 
    }
  }
}

print(stats.corr)
write.xlsx(stats.corr, file="output_rct_behav_mrs.xls", sheetName = "deltabehcorr", append = T)

# for two significant subscales
corr.fdr.p = p.adjust(stats.corr[c('bumetanide.ic._gaba_corr.CARS_i1', 'bumetanide.ic._gaba_corr.CARS_i13'),]$pvalue, method = 'BH')
print(corr.fdr.p)

# puting three together
corr.fdr.p = p.adjust(stats.corr[c('bumetanide.ic._gaba_corr.CARS_tot', 'bumetanide.ic._gaba_corr.CARS_i1', 'bumetanide.ic._gaba_corr.CARS_i13'),]$pvalue, method = 'BH')
print(corr.fdr.p)

#################
# Figure 2C
#################

# FIgure of delta correlations between  baseline and follow-up 
bumedata.baseline <- bumedata[which(bumedata$stage==0),]
bumedata.3mont <- bumedata[which(bumedata$stage==1),]
jointtable <- merge(bumedata.baseline, bumedata.3mont, by = "id")
ic_gaba_corr_delta = jointtable$ic_gaba_corr.x - jointtable$ic_gaba_corr.y
cars_tot_delta = jointtable$cars_tot.x - jointtable$cars_tot.y
corr <- data.frame(id = jointtable$id, group = jointtable$group4plot.x, 
                   ic_gaba_corr_delta =  ic_gaba_corr_delta, 
                   cars_tot_delta=cars_tot_delta,
                   QC.paired = jointtable$QC.paired.x)
corr <- corr[which(corr$QC.paired==1),]

p2 <-ggplot(corr, aes(x=cars_tot_delta, y=ic_gaba_corr_delta, color=group, shape=group)) +
  geom_point() + 
  geom_smooth(method=lm, aes(fill=group))
p2+labs(x="Decrease of CARS total scores from baseline to month 3",y="Decrease of GABA/NAA ratio from baseline to month 3")+
  theme_classic()
ggsave("Figure2C.pdf")
dev.off()
