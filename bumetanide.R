# bumetanide repeated measures mixed effect model
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


#######################
#######################
## Step 1: read data ##
#######################
#######################

## 1.1 read clinical data: in the long format
clinicalAssessment <- read.csv("bumetanide_openlabel_clinical_2.csv", header = TRUE)
colnames(clinicalAssessment)[1] <- "ID"
colnames(clinicalAssessment)[44:48] <- c('CGI_severity', 'CGI_total', 
                                         'CGI_score', 'CGI_side', 'CGI_index')

## 1.2 read MRS data: in the long format
MRS.INS <- read.xlsx("Total_MRS_20190415.xlsx", sheetName = "Insula")
MRS.OFC <- read.xlsx("Total_MRS_20190415.xlsx", sheetName = "OFC")
MRS.VC <- read.xlsx("Total_MRS_20190415.xlsx", sheetName= "VC")


##################################
# MANUALLY CHECKING IS NECESSARY #
##################################

## 1.3 clean the data

# 1) extra columns with NA 
MRS.INS <- MRS.INS[,-38]
MRS.VC <- MRS.VC[,-38]




#  2) Some odd variable names after reading into the R 
colnames(MRS.INS)
colnames(MRS.INS) <- sub('..SD','.SD', colnames(MRS.INS))
colnames(MRS.INS)<- sub("..tissue.corrected.",".tissue.corrected", colnames(MRS.INS))
colnames(MRS.INS)<- sub("tissue_corrected.","tissue.corrected", colnames(MRS.INS))
colnames(MRS.INS)<- sub("X","INS", colnames(MRS.INS))
colnames(MRS.INS) <- sub("GABA.Cr..tissue_corrected.Gannet.", "GABA.Cr.tissue.corrected.Gannet", colnames(MRS.INS))
colnames(MRS.INS)

colnames(MRS.OFC)
colnames(MRS.OFC) <- sub('..SD','.SD', colnames(MRS.OFC))
colnames(MRS.OFC)<- sub("..tissue.corrected.",".tissue.corrected", colnames(MRS.OFC))
colnames(MRS.OFC)<- sub("tissue_corrected.","tissue.corrected", colnames(MRS.OFC))
colnames(MRS.OFC) <- sub("GABA.Cr..tissue_corrected.Gannet.", "GABA.Cr.tissue.corrected.Gannet", colnames(MRS.OFC))
colnames(MRS.OFC) <- sub("X","OFC", colnames(MRS.OFC))
colnames(MRS.OFC)

colnames(MRS.VC)
colnames(MRS.VC) <- sub('..SD','.SD', colnames(MRS.VC))
colnames(MRS.VC)<- sub("..tissue.corrected.",".tissue.corrected", colnames(MRS.VC))
colnames(MRS.VC)<- sub("tissue_corrected.","tissue.corrected", colnames(MRS.VC))
colnames(MRS.VC) <- sub("GABA.Cr..tissue_corrected.Gannet.", "GABA.Cr.tissue.corrected.Gannet", colnames(MRS.VC))
colnames(MRS.VC) <- sub("X","VC", colnames(MRS.VC))
colnames(MRS.VC)


# 3) inconsistent in subject names between baseline and follow-up
levels(MRS.INS$Subjects) <- sub('-','_', levels(MRS.INS$Subjects))
levels(MRS.OFC$Subjects) <- sub('-','_', levels(MRS.OFC$Subjects))
levels(MRS.VC$Subjects) <- sub('-','_', levels(MRS.VC$Subjects))



## 1.4 join to the tables into one big table
# make sure the index for join is exactly the same in both tables before merging

# 1) join the MRS tables
MRS <- merge(MRS.INS, MRS.OFC, by=c("ID","Subjects"))
MRS <- merge(MRS, MRS.VC, by=c("ID","Subjects"))


MRS$stage <- 0
MRS$stage[grep("3M", MRS$Subjects)] = 3

MRS$ID <- sub('0_', '', MRS$ID)

# 2) join the clinical table
bumedata <- merge(clinicalAssessment, MRS, by=c("ID", "stage"), all.x=TRUE) 



# 3) check whether the variables are numeric or factor. If it is factor, we need to know which levels it has; 
# we only allow subject names and sex to be factor in this case.
# BE CAREFUL, when converting factor to numeric. Replacing the NaNwith NA in the level set before converting 
for (i in c(1:length(bumedata))){
  if(is.factor(bumedata[,i])){
    print(colnames(bumedata)[i])
    print(levels(bumedata[,i]))
  }
}
#levels(bumedata$cars_i3) <- sub('_','', levels(bumedata$cars_i3)) # it had been addressed in the raw table
#bumedata$cars_i3 <- as.numeric(as.character(bumedata$cars_i3))
#bumedata$CRI_total[34] <- 2  #It had been addressed in the raw table
#bumedata$CRI_total<- droplevels(bumedata$CRI_total)
#bumedata$sex <- droplevels(bumedata$sex) # drop the extra level of ""
# check the datatype again after the convertion
#for (i in c(1:length(bumedata))){
#  if(!is.numeric(bumedata[,i])){
#    print(colnames(bumedata)[i])
#    print(levels(bumedata[,i]))
#  }
#}

# some missing values were 999  # No signal
# df <- bumedata %>% replace_with_na_at(.vars = c(colnames(bumedata)[50:130]), condition = ~.x == 999)
# bumedata <- df
# rm(df)

# recode sex as numeric variable
bumedata$ismale <- ifelse(bumedata$sex=='M',1,0)
# write the data after cleanning 
write.xlsx(bumedata, file = "bumedata20190517.xlsx")

# remove the redundant variables
rm(MRS.INS)
rm(MRS.OFC)
rm(MRS.VC)
rm(MRS)
rm(clinicalAssessment)
#rm(bumedata.MRS)


####################################
####################################
## Step 2: Description Statistics ##
####################################
####################################


######################################################
# table 1. characteristics and background 
######################################################
# check the distribution of the data
# check if any outlier 
# CAUTION: it is not easily to define an outlier by statistic only. We need a really good reason to exclude
# some one from the following analysis. For example, this data point has failed in some quality check. Otherwise,
# we might be asked to report the results using the full data and do sensitivity analysis (maybe a bootstrap)



bumedata$group4plot <- as.factor(ifelse(bumedata$group==1,'bumetanide','control'))
bumedata$stage4plot <- as.factor(ifelse(bumedata$stage==0,'baseline','3 months'))
bumedata$stage4plot = with(bumedata, relevel(stage4plot, "baseline"))
bumedata$group4plot = with(bumedata, relevel(group4plot, "bumetanide"))
levels(bumedata$group4plot)
levels(bumedata$stage4plot)

bumedata.baseline <- bumedata[which(bumedata$stage == 0),]
bumedata.3mont <- bumedata[which(bumedata$stage == 3),]
for (i in bumedata.3mont$ID){
  bumedata.3mont$IQ[which(bumedata.3mont$ID==i)] <- bumedata.baseline$IQ[which(bumedata.baseline$ID==i)] 
}
bumedata$IQ[bumedata$stage==3] = bumedata.3mont$IQ

bumedata$time <- ifelse(bumedata$stage==3, 1,0)

# 2.1 baseline comparison
table(bumedata.baseline$group, bumedata.baseline$sex)

prop.test(table(bumedata.baseline$group, bumedata.baseline$sex), correct = FALSE)



desmat <- matrix(rep(0,16*39), ncol = 16, nrow = 39)
rownames(desmat) <- colnames(bumedata.baseline)[c(5:43)]
colnames(desmat) <- c("n","min", "max", "mean", "sd", 
                      "n","min", "max", "mean", "sd", 
                      "t.df", "t.t", "t.p", 
                      "kw.df", "kw.chi-squared", "kw.p")
for (i in 5:43){
  # descriptive table:: desmat
  description <- describeBy(bumedata.baseline[,i], group = bumedata.baseline$group, mat = TRUE, na.rm = TRUE )
  desmat[i-4,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")], 
                                         description[2, c("n","min", "max", "mean", "sd")]))
  # baseline comparison: testmat
  ttest <- t.test(bumedata.baseline[,i]~group, data = bumedata.baseline)
  desmat[i-4,11] <- ttest$parameter # degree-of-freedom
  desmat[i-4,12] <- ttest$statistic # t-statistic
  desmat[i-4,13] <- ttest$p.value #p value
  kwtest <- kruskal.test(bumedata.baseline[,i]~group, data = bumedata.baseline)
  desmat[i-4,14] <- kwtest$parameter # degree-of-freedom
  desmat[i-4,15] <- kwtest$statistic # t-statistic
  desmat[i-4,16] <- kwtest$p.value #p value
}
print(desmat)
write.xlsx(desmat, file="resultsoutputNew.xls", sheetName = "baselinebehaviour")



# 2.2 follow-up
desmat <- matrix(rep(0,16*28), ncol = 16, nrow = 28)
rownames(desmat) <- colnames(bumedata.3mont)[c(16:43)]
colnames(desmat) <- c("n","min", "max", "mean", "sd", 
                      "n","min", "max", "mean", "sd", 
                      "t.df", "t.t", "t.p", 
                      "kw.df", "kw.chi-squared", "kw.p")
for (i in 16:43){
  # descriptive table:: desmat
  description <- describeBy(bumedata.3mont[,i], group = bumedata.3mont$group, mat = TRUE, na.rm = TRUE )
  desmat[i-15,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")], 
                                         description[2, c("n","min", "max", "mean", "sd")]))
  # baseline comparison: testmat
  ttest <- t.test(bumedata.3mont[,i]~group, data = bumedata.3mont)
  desmat[i-15,11] <- ttest$parameter # degree-of-freedom
  desmat[i-15,12] <- ttest$statistic # t-statistic
  desmat[i-15,13] <- ttest$p.value #p value
  kwtest <- kruskal.test(bumedata.3mont[,i]~group, data = bumedata.3mont)
  desmat[i-15,14] <- kwtest$parameter # degree-of-freedom
  desmat[i-15,15] <- kwtest$statistic # t-statistic
  desmat[i-15,16] <- kwtest$p.value #p value
}
print(desmat)
write.xlsx(desmat, file="resultsoutputNew.xls", sheetName = "month3behaviour", append = T)


# distribution is not perfect!  but if we count for the whole data set, it is still normally distributed
bp <- list()
for (i in 16:48){
  temp <- bumedata.3mont[which(bumedata.3mont$group==0),i]
  bp[[i-15]] <- gghistogram(bumedata.3mont, colnames(bumedata.3mont)[i])
}

ggarrange(plotlist = bp, 
          ncol = 3, nrow = 11)

# individual change
bp <- list()
for (i in 16:43){
  if (i==16){
    bp[[i-15]] <- ggplot(data=bumedata, aes_string(x="stage", y=colnames(bumedata)[i], group="ID"))+
      geom_line(aes(linetype=group4plot)) + 
      geom_point(aes(shape=group4plot))  +
      labs(x = 'month') +
      theme(legend.position="top", legend.direction = 'horizontal', 
            legend.title = element_blank()) 
  }
  else {
    bp[[i-15]] <- ggplot(data=bumedata, aes_string(x="stage", y=colnames(bumedata)[i], group="ID"))+
      geom_line(aes(linetype=group4plot)) + 
      geom_point(aes(shape=group4plot)) + 
      theme(legend.position="none") +
      labs(x = '')
  }
}

ggarrange(plotlist = bp, 
          ncol = 4, nrow = 7)


######################################################
######################################################
## Step 3: Three approaches to test the drug effect ##
######################################################
######################################################

# for those with the shaptio or the levent test being significant, we need to do permutation 

# 3.1 ANCOVA: without nuisance variables : assuming a balanced design at the baseline 
stats <- data.frame(t = rep(0,28), t.p = rep(1,28), f=rep(0,28), f.p=rep(1,28))
rownames(stats) <- colnames(bumedata.3mont[16:43])
for (i in 16:43){
  fit <- lm(bumedata.3mont[,i]~ bumedata.3mont$group + bumedata.baseline[,i])
  stats$shapitestp[i-15] <- shapiro.test(residuals(fit))$p.value
  stats$leventestp[i-15] <- leveneTest(bumedata.3mont[,i]~ as.factor(bumedata.3mont$group))[1,3]
  stats$t[i-15] <- summary(fit)$coefficients["bumedata.3mont$group","t value"]
  stats$t.p[i-15] <- summary(fit)$coefficients["bumedata.3mont$group","Pr(>|t|)"]
  ftests <- anova(fit)
  stats$f[i-15] <- ftests["bumedata.3mont$group","F value"]
  stats$f.p[i-15] <- ftests["bumedata.3mont$group","Pr(>F)"]
}
print(stats)

write.xlsx(stats, file="resultsoutputNew.xls", sheetName = "drugeffectonbehaviourModel1", append = T)

# 3.2 ANCOVA: with nuisance variables 
stats2 <- data.frame(t = rep(0,28), t.p = rep(1,28), f=rep(0,28), f.p=rep(1,28))
rownames(stats2) <- colnames(bumedata.3mont[16:43])
bumedata.3mont$ismale <- ifelse(bumedata.3mont$sex=="M", 1,0)
for (i in 16:43){
  fit <- lm(bumedata.3mont[,i]~ bumedata.3mont$group + bumedata.3mont$age 
            + bumedata.3mont$ismale + bumedata.baseline$IQ + bumedata.baseline[,i])
  stats2$shapitestp[i-15] <- shapiro.test(residuals(fit))$p.value
  stats2$leventestp[i-15] <- leveneTest(bumedata.3mont[,i]~ as.factor(bumedata.3mont$group))[1,3]
  stats2$t[i-15] <- summary(fit)$coefficients["bumedata.3mont$group","t value"]
  stats2$t.p[i-15] <- summary(fit)$coefficients["bumedata.3mont$group","Pr(>|t|)"]
  ftests <- anova(fit)
  stats2$f[i-15] <- ftests["bumedata.3mont$group","F value"]
  stats2$f.p[i-15] <- ftests["bumedata.3mont$group","Pr(>F)"]
}
print(stats2)
write.xlsx(stats2, file="resultsoutputNew.xls", sheetName = "drugeffectonbehaviourModel2", append = T)


# 3.3 LMM: liner mixed effect model with nuisance variables
stats3 <- data.frame(t = rep(0,28), t.p = rep(1,28), f=rep(0,28), f.p=rep(1,28), f.df=rep(1,28))
rownames(stats3) <- colnames(bumedata[16:43])
for (i in 16:43){
  score <- bumedata[,i]
  fit <- lmer( score ~  time + group + group*time + ismale + age + IQ + (1|ID), data=bumedata)
  stats3$shapitestp[i-15] <- shapiro.test(residuals(fit))$p.value
  stats3$leventestp[i-15] <- leveneTest(score ~ as.factor(group*time), data=bumedata)[1,3]
  stats3$t[i-15] <- summary(fit)$coefficients["time:group","t value"]
  stats3$t.p[i-15] <- summary(fit)$coefficients["time:group","Pr(>|t|)"]
  stats3$t.df[i-15] <- summary(fit)$coefficients["time:group","df"]
  ftests <- anova(fit, type=1, ddf="Kenward-Roger")
  stats3$f[i-15] <- ftests["time:group","F value"]
  stats3$f.p[i-15] <- ftests["time:group","Pr(>F)"]
  stats3$f.df[i-15] <- ftests["time:group","DenDF"]
}
print(stats3)
write.xlsx(stats3, file="resultsoutputNew.xls", sheetName = "drugeffectonbehaviourModel3", append = T)

# cars_tot
i = 21
# rerun the above fit to identify the outlier by qqplot of the model residuals
qqnorm(residuals(fit))
outlier <- which(abs(residuals(fit))>2)

score <- bumedata[-outlier,i]
fit <- lmer( score ~  time + group + group*time + ismale + age + IQ + (1|ID), data=bumedata[-outlier,])
stats3$shapitestp[i-15] <- shapiro.test(residuals(fit))$p.value
stats3$leventestp[i-15] <- leveneTest(score ~ as.factor(group*time), data=bumedata[-outlier,])[1,3]
stats3$t[i-15] <- summary(fit)$coefficients["time:group","t value"]
stats3$t.p[i-15] <- summary(fit)$coefficients["time:group","Pr(>|t|)"]
stats3$t.df[i-15] <- summary(fit)$coefficients["time:group","df"]
ftests <- anova(fit, type=1, ddf="Kenward-Roger")
stats3$f[i-15] <- ftests["time:group","F value"]
stats3$f.p[i-15] <- ftests["time:group","Pr(>F)"]
stats3$f.df[i-15] <- ftests["time:group","DenDF"]
print(stats3[i-15,])




# permutation test for subscales
set.seed(300)
currentdata <- bumedata
perm.p.cars <- matrix(rep(0,17), ncol=1, nrow=17)
rownames(perm.p.cars) <- colnames(currentdata)[c(21:37)]

for (i in rownames(perm.p.cars)){
  score <- (currentdata[[i]])
  # linear model comparison
  fit1 <- lmer(score ~ time * group + age + ismale + IQ + (1|ID), data = currentdata)
  fit0 <- lmer(score ~ time + group + age + ismale + IQ + (1|ID), data = currentdata)
  
  perm.p.cars[[i,1]] <- permlmer(fit0, fit1, perms = 3000, ncore = 8, plot = FALSE)$"Perm-p"[2]
}  
print(perm.p.cars)
perm.p.cars.fdr <- perm.p.cars
perm.p.cars.fdr[c(3:17)] <- p.adjust(perm.p.cars[c(3:17)], method = "fdr")

write.xlsx(cbind(perm.p.cars, perm.p.cars.fdr), file="resultsoutputNew.xls", 
           sheetName = "drugbehModel3PermFDR", append = T)


##########
# Figure 2: significant behaviour
##########
sig.beh <- c("cars_tot", "item_s3")

# group comparison
bp <- list()
s <- 1
for (i in sig.beh){
  bp[[s]] <- ggplot(data=bumedata, aes_string(x = "group4plot", y = i, fill = "stage4plot")) + 
    geom_boxplot(position=position_dodge(1)) +  
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(1)) + 
    labs(x = '') + theme(legend.position = 'none') 
  if (s == 3){
    bp[[s]] = bp[[s]] + theme(legend.position = c(0.40,0.90), legend.title = element_blank(), 
                              legend.key = element_rect(fill = "transparent", colour = NA),
                              legend.background = element_rect(fill = "transparent", colour = NA),     # get rid of legend bg
                              legend.box.background = element_rect(fill = "transparent", colour = NA)) # get rid of legend panel bg)
  }
  if (s == 3){
    bp[[s]] <- bp[[s]] + labs(x = 'group') 
  }
  s <- s + 2
}
# individual change
s <- 2
for (i in sig.beh){
  bp[[s]] <- ggplot(data=bumedata, aes_string(x="stage", y=i, group="ID")) +
    geom_line(aes(linetype=group4plot)) + 
    geom_point(aes(shape=group4plot))  +
    labs(x = '') + labs(y = '') + theme(legend.position = 'none') 
  if (s == 4){
    bp[[s]] <- bp[[s]] + labs(x = 'month') 
  }
  if (s == 4){
    bp[[s]] <- bp[[s]] + theme(legend.position = c(0.70,0.95), 
                               legend.title = element_blank(),
                               #panel.background = element_rect(fill = "transparent", colour = NA), # bg of the panel
                               #plot.background = element_rect(fill = "transparent", colour = NA), # bg of the plot
                               #panel.grid.major = element_blank(), # get rid of major grid
                               #panel.grid.minor = element_blank(), # get rid of minor grid
                               legend.key = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),     # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent", colour = NA)) # get rid of legend panel bg
  }
  
  s <- s + 2
}

pdf("Figure2.pdf",width=6,height=4,paper='special')
ggarrange(plotlist = bp, 
          ncol = 2, nrow = 2)
dev.off()



# 3.4 compare CRI indeces

# count
table(bumedata.3mont$CRI_score, bumedata.3mont$group)
table(bumedata.3mont$CRI_index, bumedata.3mont$group)
table(bumedata.3mont$CRI_severity, bumedata.3mont$group)
table(bumedata.3mont$CRI_total, bumedata.3mont$group)
table(bumedata.3mont$CRI_side, bumedata.3mont$group)

# nonparametric test
stats.cri <- matrix(rep(0, 20), ncol = 4, nrow = 5)
rownames(stats.cri) <- colnames(bumedata.3mont)[c(44:48)]
colnames(stats.cri) <- c("degree of freedom", "Kruskal-Wallis chi-squared", "p-value", "fdr-p")
# wilcox.test(CRI_score~group, data=bumedata.3mont)  failed to converge
for (i in c(44:48)){
  kstest <- kruskal.test(bumedata.3mont[,i]~group, data=bumedata.3mont)  # significant
  stats.cri[i-43,c(1:3)] <- c(kstest$parameter, kstest$statistic, kstest$p.value)
}

stats.cri[,"fdr-p"] <- t(p.adjust(stats.cri[,"p-value"], method = "fdr"))
print(stats.cri)
write.xlsx(stats.cri, file="resultsoutputNew.xls", sheetName = "month3CGI", append = T)


# permutation to confirm
# library(coin)

# independence_test(CRI_score~group, data=bumedata.3mont) # sign
# independence_test(CRI_index~group, data=bumedata.3mont) # sign
# independence_test(CRI_total~group, data=bumedata.3mont) # sign

##########################
##########################
## Step 4: MRS analysis ##
##########################
##########################


# 4.0 quality control

# NAA+NAAG %SD < 20  &  FWHM < 0.05 & SNR > 15
bumedata$QC <- matrix(rep(0,3*dim(bumedata)[1]), nrow = dim(bumedata)[1], ncol = 3)
ROI.names = c("INS","OFC", "VC")
FWHM.names = c("FWHM.x", "FWHM.y", "FWHM")
SNR.names = c("SNR.x", "SNR.y", "SNR")
des.QC <- matrix(rep(0,6*3), ncol = 2*3, nrow = 3)
rownames(des.QC) <- ROI.names
colnames(des.QC) <- c("SD.mean", "SD.sd","FW.mean", "FW.sd","SNR.mean", "SNR.sd")
for (i in c(1:3)){
  NAAplusSD = paste0(ROI.names[i], '.NAA.NAAG.SD')
  bumedata$QC[,i] <- (bumedata[[NAAplusSD]] < 20 & bumedata[[FWHM.names[i]]] < 0.05 & bumedata[[SNR.names[i]]] > 15)
  des.QC[i,c(1,2)] <- as.matrix(describe(bumedata[[NAAplusSD]]))[c(3,4)]
  des.QC[i,c(3,4)] <- as.matrix(describe(bumedata[[FWHM.names[i]]]))[c(3,4)]
  des.QC[i,c(5,6)] <- as.matrix(describe(bumedata[[SNR.names[i]]]))[c(3,4)]
}
counts <- matrix(rep(0,3*3), nrow = 3, ncol = 3, dimnames = list(c("baseline", "3 month", 'paired'), ROI.names))
counts[1,] <- colSums(bumedata$QC[which(bumedata$stage==0),], na.rm =  T)
counts[2,] <- colSums(bumedata$QC[which(bumedata$stage==3),], na.rm =  T)

mrsqc <- data.frame(cbind(bumedata$ID, bumedata$QC, bumedata$stage))
colnames(mrsqc) <- c('ID', 'INS.QC', 'OFC.QC', 'VC.QC', 'stage')
mrsqc_paired <- merge(mrsqc[which(mrsqc$stage==0),], mrsqc[which(mrsqc$stage==3),], by = "ID")
mrsqc_paired$INS.QC <- ((mrsqc_paired$INS.QC.x == mrsqc_paired$INS.QC.y) & mrsqc_paired$INS.QC.x == 1)
mrsqc_paired$OFC.QC <- ((mrsqc_paired$OFC.QC.x == mrsqc_paired$OFC.QC.y) & mrsqc_paired$OFC.QC.x == 1)
mrsqc_paired$VC.QC <- ((mrsqc_paired$VC.QC.x == mrsqc_paired$VC.QC.y) & mrsqc_paired$VC.QC.x == 1)
bumedata$QC.paired <- matrix(rep(0,3*dim(bumedata)[1]), nrow = dim(bumedata)[1], ncol = 3)
for (i in c(1:dim(bumedata)[1])){
  bumedata$QC.paired[i,] <- as.matrix(mrsqc_paired[which(mrsqc_paired$ID == bumedata$ID[i]), c(10:12)])
}
counts[3,] <- colSums(bumedata$QC.paired, na.rm = T)/2
print(counts) 
write.xlsx(counts, file="resultsoutputNew.xls", sheetName = "MRSafterQC", append = T)



# 4.1 data description 
colnames(bumedata)[which(colnames(bumedata)=="GABA.NAA.tissue.corrected.x")] <- "INS.GABA.NAA.tissue.corrected"
colnames(bumedata)[which(colnames(bumedata)=="NAA.tissue.corrected.x")] <- "INS.NAA.tissue.corrected"
colnames(bumedata)[which(colnames(bumedata)=="GABA.Glx.tissue.corrected.x")] <- "INS.GABA.Glx.tissue.corrected"

colnames(bumedata)[which(colnames(bumedata)=="GABA.NAA.tissue.corrected.y")] <- "OFC.GABA.NAA.tissue.corrected"
colnames(bumedata)[which(colnames(bumedata)=="NAA.tissue.corrected.y")] <- "OFC.NAA.tissue.corrected"
colnames(bumedata)[which(colnames(bumedata)=="GABA.Glx.tissue.corrected.y")] <- "OFC.GABA.Glx.tissue.corrected"

colnames(bumedata)[which(colnames(bumedata)=="GABA.NAA.tissue.corrected")] <- "VC.GABA.NAA.tissue.corrected"
colnames(bumedata)[which(colnames(bumedata)=="NAA.tissue.corrected")] <- "VC.NAA.tissue.corrected"
colnames(bumedata)[which(colnames(bumedata)=="GABA.Glx.tissue.corrected")] <- "VC.GABA.Glx.tissue.corrected"

currentdata = bumedata[bumedata$stage==0,]
ROI.names = c("INS","VC")
MRS.trans = c(".GABA.NAA.tissue.corrected", ".NAA.tissue.corrected", ".GABA.Glx.tissue.corrected")
desmat <- matrix(rep(0,6*19), ncol = 19, nrow = 6)
colnames(desmat) <- c("n","min", "max", "mean", "sd", "n2","min2", "max2", 
                      "mean2", "sd2", 
                      "shapi.p", "levene.p", "t", "t.p", "F", "F.p", 
                      "F.Df", "ks.p", "perm.p")
rownames(desmat) <- c(1:6)
tableCount <- matrix(rep(0,length(MRS.trans)*length(ROI.names)), ncol=length(MRS.trans), nrow=length(ROI.names))
dimnames(tableCount) <- list(ROI.names, MRS.trans)
testmat <- list()
bp <- list()
s <- 0
index = matrix(c(1,3), ncol = 1, nrow = 2) # this was nedeed as we had droped OFC 
rownames(index) <- ROI.names
t <- tableCount
for (i in ROI.names){
  for (j in MRS.trans){ 

    s <- s + 1
    
    NOI <- paste0(i,j)
    validity <- (currentdata$QC[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
    tableCount[i,j] = sum(validity==TRUE, na.rm=T)
    score <- (currentdata[[NOI]][validity==1])
    # descriptive table:: desma
    description <- describeBy(score, group = currentdata$group[validity==1], mat = TRUE, na.rm = TRUE )
    rownames(desmat)[s] <- NOI
    desmat[s,c(1:10)] <- as.matrix(cbind(description[1, c("n","min", "max", "mean", "sd")], 
                                            description[2, c("n","min", "max", "mean", "sd")]))
    
    # baseline comparison: testmat
    
    # linear model comparison
    fit <- lm(score~group + age + ismale + IQ, data = currentdata[validity==1,])
    desmat[s,11] <- shapiro.test(residuals(fit))$p.value
    desmat[s,12]<- leveneTest(score ~ as.factor(group), data=currentdata[validity==1,])[1,3]
    desmat[s,13] <- summary(fit)$coefficients["group","t value"]
    t[i,j] <- desmat[s,13]
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
}

ggarrange(plotlist = bp, 
          ncol = length(MRS.trans), nrow = length(ROI.names))



# permutation
set.seed(300)
perm.count <- matrix(rep(0,length(MRS.trans)*length(ROI.names)), ncol=length(MRS.trans), nrow=length(ROI.names))
dimnames(perm.count) <- list(ROI.names, MRS.trans)
perm.t <- perm.count
perm.total <- 3000
for (nperm in c(1:perm.total)){
  perm.sample <- sample.int(dim(currentdata)[1])
  currentdata$perm.group <- currentdata$group[perm.sample]
  st <- 0
  for (i in ROI.names){
    st <- st + 1
    for (j in MRS.trans){
      NOI <- paste0(i,j)
      validity <- (currentdata$QC[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3)
      score <- currentdata[[NOI]][validity==1]
      # linear model comparison
      fit <- lm(score~perm.group + age + ismale + IQ, data = currentdata[validity==1,])
      perm.t[i,j] <- summary(fit)$coefficients["perm.group","t value"]
    }
  }
  for (i in ROI.names){
    for (j in MRS.trans){
      if (j == ".NAA.tissue.corrected"){
        if ((abs(perm.t[i,j])) > abs(t[i,j])){
          perm.count[i,j] <- perm.count[i,j] + 1
        }
      }
      else{
        if (max(abs(perm.t[,-which(MRS.trans == ".NAA.tissue.corrected")])) > abs(t[i,j])){
          perm.count[i,j] <- perm.count[i,j] + 1
        }
      }
    }
  }
}
print(perm.count/perm.total)
desmat[,19] <- unmatrix(perm.count/perm.total, byrow = T)
fdrcorrected <- p.adjust(desmat[-c(2,5),19], method = "fdr")


print(desmat)
write.xlsx(desmat, file="resultsoutputNew.xls", sheetName = "MRSbaseline", append = T)


# distribution is not perfect!  but if we count for the whole data set, it is still normally distributed


#function that takes in vector of data and a coefficient,
#returns boolean vector if a certain point is an outlier or not
check_outlier <- function(v, coef=1.5){
  quantiles <- quantile(v,probs=c(0.25,0.75))
  IQR <- quantiles[2]-quantiles[1]
  res <- v < (quantiles[1]-coef*IQR)|v > (quantiles[2]+coef*IQR)
  return(res)
}


bp <- list()
outlierstable <- list()
s = 0
st <- 0
for (i in ROI.names){
  st <- st + 1
  for (j in MRS.trans){
    NOI <- paste0(i,j)
    s = s + 1
    #outlierstable[[s]] <- bumedata$ID[which(abs(scale(bumedata[[NOI]])) >= 3)]
    #validity <- (bumedata[[paste0(i,j,'.SD')]]<20 & abs(scale(bumedata[[NOI]])) < 3)
    validity <- (bumedata$QC[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
    
    dat <- bumedata[which(validity==TRUE),]
    dat$label <- dat$Subjects
    levels(dat$label) <- c(levels(dat$label), "")
    for (m in levels(dat$group4plot)){
      for (n in levels(dat$stage4plot)){
        groupidx <- which(dat$group4plot == m & dat$stage4plot == n)
        groupdat = dat[[NOI]][groupidx]
        dat$label[groupidx[!check_outlier(groupdat)]] <- ""
      }
    }
    
    bp[[s]] <- ggplot(dat, aes_string(x = "group4plot", y = NOI, fill = "stage4plot")) + 
               geom_boxplot(position=position_dodge(1)) +  
               geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(1)) + 
               geom_text(aes(label=label), size=3) +
               theme(legend.position = c(0.85,0.85))
    if (s != 3){
      bp[[s]] = bp[[s]] + theme(legend.position="none") 
    }
      
  }
}

ggarrange(plotlist = bp, 
          ncol = length(MRS.trans), nrow = length(ROI.names)) 

# individual change

bp <- list()
s = 0
st <- 0
for (i in ROI.names){
  st <- st + 1
  for (j in MRS.trans){
    NOI <- paste0(i,j)
    s = s + 1
    #validity <- (bumedata[[paste0(i,j,'.SD')]]<20 & abs(scale(bumedata[[NOI]])) < 3)
    #validity[bumedata$stage==0] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
    #validity[bumedata$stage==3] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
    validity <- (bumedata$QC.paired[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
    dat <- bumedata[which(validity==TRUE),]
    
    bp[[s]] <- ggplot(data=dat, aes_string(x="stage", y=NOI, group="ID")) +
      geom_line(aes(linetype=group4plot)) + 
      geom_point(aes(shape=group4plot))  +
      theme(legend.position=c(0.7,0.9), #legend.direction = 'horizontal', 
            legend.background = element_rect(fill = NA, colour= NA),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.title = element_blank()) + labs(x = 'month')
    if (s!=3) {
      bp[[s]] = bp[[s]] + theme(legend.position="none") 
    }
    if (s!=5) {
      bp[[s]] = bp[[s]] + labs(x = '')
    }
  }
}
ggarrange(plotlist = bp, 
          ncol = length(MRS.trans), nrow = length(ROI.names))  

dev.off()

# 4.2 test drug effect on MRS
# ancova
# # 1) no covariates
# t <- matrix(rep(0, length(ROI.names) * length(MRS.trans)), nrow = length(ROI.names), ncol = length(MRS.trans))
# dimnames(t) <- list(ROI.names, MRS.trans)
# f <- t
# t.p <- t + 1
# f.p <- t.p
# 
# for (i in ROI.names){
#   for (j in MRS.trans){
#     NOI <- paste0(i,j)
#     #validity <- (bumedata[[paste0(i,j,'.SD')]]<20 & abs(scale(bumedata[[paste0(i,j)]])) < 3 &
#     #               !is.na(bumedata[[paste0(i,j)]]) & bumedata$sex == 'M')
#     #validity[bumedata$stage==0] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
#     #validity[bumedata$stage==3] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
#     validity <- (bumedata$QC.paired[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
# 
#     dat <- bumedata[validity==TRUE,]
#     fit <- lm((dat[[NOI]][dat$stage==3]) ~ dat$group[dat$stage==3] + dat[[NOI]][dat$stage==0])
#     t[i,j] <- summary(fit)$coefficients["dat$group[dat$stage == 3]","t value"]
#     t.p[i,j] <- summary(fit)$coefficients["dat$group[dat$stage == 3]","Pr(>|t|)"]
#     ftests <- anova(fit)
#     f[i,j] <- ftests["dat$group[dat$stage == 3]","F value"]
#     f.p[i,j] <- ftests["dat$group[dat$stage == 3]","Pr(>F)"]
#   }
# }
# stats.mrs <- data.frame(t = t, t.p = t.p, f=f, f.p=f.p)
# print(stats.mrs)
# write.xlsx(stats.mrs, file="resultsoutput.xls", sheetName = "drugeffectMRSmodel1", append = T)
# 
# # 2) covariates
# t <- matrix(rep(0, length(ROI.names) * length(MRS.trans)), nrow = length(ROI.names), ncol = length(MRS.trans))
# dimnames(t) <- list(ROI.names, MRS.trans)
# f <- t
# t.p <- t + 1
# f.p <- t.p
# for (i in ROI.names){
#   for (j in MRS.trans){
#     NOI <- paste0(i,j)
#     #validity <- (bumedata[[paste0(i,j,'.SD')]]<20 & abs(scale(bumedata[[paste0(i,j)]])) < 3 &
#     #               !is.na(bumedata[[paste0(i,j)]]) & bumedata$sex == "M" )
#     #validity[bumedata$stage==0] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
#     #validity[bumedata$stage==3] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
#     
#     validity <- (bumedata$QC.paired[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
#     
#     dat <- bumedata[validity==TRUE,]
#     fit <- lm((dat[[NOI]][dat$stage==3]) ~ dat$group[dat$stage==3] + dat[[NOI]][dat$stage==0] + 
#                  dat$age[dat$stage==0] + dat$IQ[dat$stage==0] )
#     t[i,j] <- summary(fit)$coefficients["dat$group[dat$stage == 3]","t value"]
#     t.p[i,j] <- summary(fit)$coefficients["dat$group[dat$stage == 3]","Pr(>|t|)"]
#     ftests <- anova(fit)
#     f[i,j] <- ftests["dat$group[dat$stage == 3]","F value"]
#     f.p[i,j] <- ftests["dat$group[dat$stage == 3]","Pr(>F)"]
#   }
# }
# stats2.mrs <- data.frame(t = t, t.p = t.p, f=f, f.p=f.p)
# print(stats2.mrs)
# write.xlsx(stats2.mrs, file="resultsoutput.xls", sheetName = "drugeffectMRSmodel2", append = T)

# 3) LMM
t <- matrix(rep(0, length(ROI.names) * length(MRS.trans)), nrow = length(ROI.names), ncol = length(MRS.trans))
dimnames(t) <- list(ROI.names, MRS.trans)
shapitestp <- t
leventestp <- t
f <- t
t.p <- t + 1
t.df <- t
f.p <- t.p
f.df <- t
for (i in ROI.names){
  for (j in MRS.trans){
    
    NOI <- paste0(i,j)
    print(NOI)
    #validity <- (bumedata[[paste0(i,j,'.SD')]]<20 & abs(scale(bumedata[[paste0(i,j)]])) < 3 &
    #               !is.na(bumedata[[paste0(i,j)]]) & bumedata$sex == 'M')
    #validity[bumedata$stage==0] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
    #validity[bumedata$stage==3] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
    
    validity <- (bumedata$QC.paired[,index[i,1]] == 1 & !is.na(bumedata$QC.paired[,index[i,1]])) #(abs(scale(currentdata[[NOI]])) < 3) 
   
    
    dat <- bumedata[validity==TRUE,]
    score <- (dat[[NOI]])
    baseline.score <- score
    baseline.score[dat$stage==3] <- baseline.score[dat$stage==0]
    fit <- lmer( score ~  time + group + group*time + baseline.score + age + ismale + IQ + (1|ID), data=dat)
    shapitestp[i,j] <- shapiro.test(residuals(fit))$p.value
    leventestp[i,j] <- leveneTest(score ~ as.factor(group), data=bumedata[validity==1,])[1,3]
    t[i,j] <- summary(fit)$coefficients["time:group","t value"]
    t.p[i,j] <- summary(fit)$coefficients["time:group","Pr(>|t|)"]
    t.df[i,j] <- summary(fit)$coefficients["time:group","df"]
    ftests <- anova(fit, type=1, ddf="Kenward-Roger")
    f[i,j] <- ftests["time:group","F value"]
    f.p[i,j] <- ftests["time:group","Pr(>F)"]
    f.df[i,j] <- ftests["time:group","DenDF"]
  }
}
#fdr.p <- matrix(p.adjust(unmatrix(f.p), method = "fdr"), ncol = 3, nrow = 2)
#fdr.p <- matrix(p.adjust(unmatrix(f.p[,c(1,3)]), method = "fdr"), ncol = 2, nrow = 2)



# library(predictmeans)
# Oats$nitro <- factor(Oats$nitro)
#fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
#fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
# permlmer(fm0, fm)

# permutation test
set.seed(300)
perm.p <- matrix(rep(0, length(ROI.names) * length(MRS.trans)), nrow = length(ROI.names), ncol = length(MRS.trans))
dimnames(perm.p) <- list(ROI.names, MRS.trans)
currentdata = bumedata
st <- 0
for (i in ROI.names){
  st <- st + 1
  for (j in MRS.trans){
    NOI <- paste0(i,j)
    print(NOI)
    validity <- (currentdata$QC.paired[,index[i,1]] == 1 & !is.na(currentdata$QC.paired[,index[i,1]])) #(abs(scale(currentdata[[NOI]])) < 3) 
    score <- (currentdata[[NOI]][validity==1])
    baseline.score <- score
    baseline.score[currentdata$stage[validity==1]==3] <- baseline.score[currentdata$stage[validity==1]==0]
    # linear model comparison
    fit1 <- lmer(score ~ time + group + time * group + baseline.score + age + ismale + IQ + (1|ID), data = currentdata[validity==1,])
    fit0 <- lmer(score ~ time + group + + baseline.score + age + ismale + IQ + (1|ID), data = currentdata[validity==1,])
    perm.p[i,j] <- permlmer(fit0, fit1, perms = 3000, ncore = 8, plot = FALSE)$"Perm-p"[2]
  }
}

  
#perm.fdr.p <- matrix(p.adjust(unmatrix(perm.p), method = "fdr"), ncol = 3, nrow = 2)
perm.fdr.p <- matrix(p.adjust(unmatrix(perm.p[,c(1,3)]), method = "fdr"), ncol = 2, nrow = 2)

stats3.mrs <- data.frame(t = t, t.p = t.p, f=f, f.p=f.p, f.df = f.df,
                         shapitestp = shapitestp, leventestp = leventestp, 
                         perm.p = perm.p, fdr.perm.p = perm.fdr.p)
print(stats3.mrs)
write.xlsx(stats3.mrs, file="resultsoutputNew.xls", sheetName = "DEMRSM3baselineCov", append = T)

#################
# Figure 3
#################
ylabel <- c("INS.GABA", "", "INS.GABA/Glx", "", "VC.GABA/Glx", "")
bp <- list()
outlierstable <- list()
s = 1
i <- ROI.names[1]
for (j in MRS.trans[c(1,3)]){
  NOI <- paste0(i,j)
  validity <- (bumedata$QC[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
  
  dat <- bumedata[which(validity==TRUE),]
  dat$label <- dat$Subjects
  levels(dat$label) <- c(levels(dat$label), "")
  for (m in levels(dat$group4plot)){
    for (n in levels(dat$stage4plot)){
      groupidx <- which(dat$group4plot == m & dat$stage4plot == n)
      groupdat = dat[[NOI]][groupidx]
      dat$label[groupidx[!check_outlier(groupdat)]] <- ""
    }
  }
  
  bp[[s]] <- ggplot(dat, aes_string(x = "group4plot", y = NOI, fill = "stage4plot")) + 
    geom_boxplot(position=position_dodge(1)) +  
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(1)) + 
    theme(legend.position="none") + labs(x = '', y = ylabel[s]) #geom_text(aes(label=label), size=3) + 
  if (s == 5){
    bp[[s]] = bp[[s]] + labs(x = 'group') +
      theme(legend.position = c(0.80,0.99), legend.title = element_blank(), 
            legend.key = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),     # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent", colour = NA)) # get rid of legend panel bg)
    
  }
  s <- s + 2
}


# individual change

s = 2
for (j in MRS.trans[c(1,3)]){
  NOI <- paste0(i,j)
  validity <- (bumedata$QC.paired[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
  dat <- bumedata[which(validity==TRUE),]
  
  bp[[s]] <- ggplot(data=dat, aes_string(x="stage", y=NOI, group="ID")) +
    geom_line(aes(linetype=group4plot)) + 
    geom_point(aes(shape=group4plot)) + labs(x = '', y = '') + theme(legend.position="none") 
  if (s==6) {
    bp[[s]] = bp[[s]] + labs(x = 'month') +
      theme(legend.position = c(0.70,0.99), legend.title = element_blank(), 
            legend.key = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),     # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.spacing.y = unit(0.1, 'cm')) # get rid of legend panel bg)
  }
  s = s + 2
}

s = 5
i <- ROI.names[1]
for (j in MRS.trans[3]){
  NOI <- paste0(i,j)
  validity <- (bumedata$QC[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
  
  dat <- bumedata[which(validity==TRUE),]
  dat$label <- dat$Subjects
  levels(dat$label) <- c(levels(dat$label), "")
  for (m in levels(dat$group4plot)){
    for (n in levels(dat$stage4plot)){
      groupidx <- which(dat$group4plot == m & dat$stage4plot == n)
      groupdat = dat[[NOI]][groupidx]
      dat$label[groupidx[!check_outlier(groupdat)]] <- ""
    }
  }
  
  bp[[s]] <- ggplot(dat, aes_string(x = "group4plot", y = NOI, fill = "stage4plot")) + 
    geom_boxplot(position=position_dodge(1)) +  
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(1)) + 
    theme(legend.position="none") + labs(x = '', y = ylabel[s]) #geom_text(aes(label=label), size=3) + 
  if (s == 5){
    bp[[s]] = bp[[s]] + labs(x = 'group') +
      theme(legend.position = c(0.70,0.99), legend.title = element_blank(), 
            legend.key = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),     # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent", colour = NA)) # get rid of legend panel bg)
    
  }
  s <- s + 2
}


# individual change

s = 6
for (j in MRS.trans[3]){
  NOI <- paste0(i,j)
  validity <- (bumedata$QC.paired[,index[i,1]] == 1) #(abs(scale(currentdata[[NOI]])) < 3) 
  dat <- bumedata[which(validity==TRUE),]
  
  bp[[s]] <- ggplot(data=dat, aes_string(x="stage", y=NOI, group="ID")) +
    geom_line(aes(linetype=group4plot)) + 
    geom_point(aes(shape=group4plot)) + labs(x = '', y = '') + theme(legend.position="none") 
  if (s==6) {
    bp[[s]] = bp[[s]] + labs(x = 'month') +
      theme(legend.position = c(0.70,0.99), legend.title = element_blank(), 
            legend.key = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),     # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.spacing.y = unit(0.05, 'mm')) # get rid of legend panel bg)
  }
  s = s + 2
}

pdf("Figure3.pdf",width=6,height=4,paper='special')
ggarrange(plotlist = bp, 
          ncol = 2, nrow = 3)
dev.off()
 

# 4.3 correlations between behaviour and MRS association

# 1) cross-sectional correlations at both baseline and follow-up 
stage = levels(bumedata$stage4plot)[c(2,1)]
group = levels(bumedata$group4plot)
beh = colnames(bumedata)[c(16:37)]
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
    for (i in ROI.names){
      for (j in MRS.trans){

        NOI <- paste0(i,j)
        
        #validity <- (bumedata[[paste0(i,j,'.SD')]]<20 & abs(scale(bumedata[[paste0(i,j)]])) < 3 &
        #               !is.na(bumedata[[paste0(i,j)]]))
        #validity[bumedata$stage==0] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
        #validity[bumedata$stage==3] <- ifelse(validity[bumedata$stage==0] & validity[bumedata$stage==3], TRUE, FALSE)
        validity <- bumedata$QC[,index[i,1]]
        dat <- bumedata[validity==TRUE & bumedata$group4plot == g & bumedata$stage4plot == k, ]
        covariates <- cbind(bumedata$age, bumedata$sex, bumedata$IQ)
        covariates <- covariates[validity==TRUE & bumedata$group4plot == g & bumedata$stage4plot == k, ]
        for (b in beh){
          #print(s)
          s = s + 1
          lookupindex[k,g,i,j,b] <- s 
          rownames(stats.corr)[s] <- paste(k,g,i,j,b,sep = ",")
          datanow <- cbind(dat[[NOI]],dat[[b]], covariates)
          if (length(which(rowSums(is.na(datanow))>0)) > 0){
            datanow <- datanow[-which(rowSums(is.na(datanow))>0),]
          }
          if (var(datanow[,4])==0){
            r <- pcor.test(datanow[,1],datanow[,2], datanow[,c(3,5)], method = "spearman")
            
          }else{
            r <- pcor.test(datanow[,1],datanow[,2], datanow[,c(3:5)], method = "spearman")
            
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
  }
}

# make a table
write.xlsx(stats.corr, file="resultsoutputNew.xls", sheetName = "crosssectionalbehcorr", append = T)





# 2) delta correlations between  baseline and follow-up 

#stage = levels(bumedata$stage4plot)
#group = levels(bumedata$group4plot)
#beh = colnames(bumedata)[c(16:43)]
dat <- data.frame()
lookupindex <- array(rep(0,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)), dim = c(length(group),length(ROI.names),length(MRS.trans),length(beh)), 
                     dimnames=list(group, ROI.names, MRS.trans, beh))
stats.corr <- data.frame(r = rep(0,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)), 
                         r.p = rep(1,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)), 
                         n = rep(0,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)))
significant <- matrix(rep(0,length(group)*length(ROI.names)*length(MRS.trans)*length(beh)), 
                      nrow = 1, ncol = length(group)*length(ROI.names)*length(MRS.trans)*length(beh))
s = 0
for (g in group){
  for (i in ROI.names){
    for (j in MRS.trans){
        NOI <- paste0(i,j)
        

        validity <- bumedata$QC.paired[,index[i,1]]
        dat <- bumedata[validity==TRUE,]
        
        baseline.mrs <- dat[[NOI]][dat$stage4plot == "baseline" & dat$group4plot == g]
        delta.mrs <- baseline.mrs - dat[[NOI]][dat$stage4plot == "3 months" & dat$group4plot == g]
        covariates <- cbind(dat$age, dat$sex, dat$IQ)
        covariates <- covariates[dat$stage4plot == "baseline" & dat$group4plot == g,]
        
        for  (b in beh){
          s  = s + 1
          lookupindex[g,i,j,b] <- s 
          rownames(stats.corr)[s] <- paste(g,i,j,b,sep = ",")
          baseline.beh <- dat[[b]][dat$stage4plot=="baseline" & dat$group4plot == g]
          delta.beh <- baseline.beh - dat[[b]][dat$stage4plot == "3 months" & dat$group4plot == g]
          
          datanow <- cbind(baseline.mrs, delta.mrs, baseline.beh, delta.beh, covariates)
          if (length(which(rowSums(is.na(datanow))>0)) > 0){
            datanow <- datanow[-which(rowSums(is.na(datanow))>0),]
          }
          if (var(datanow[,6])==0){
            r <- pcor.test(datanow[,2],datanow[,4], datanow[,c(1,3,5,7)], method = "spearman")
            
          }else{
            r <- pcor.test(datanow[,2],datanow[,4], datanow[,c(1,3,5:7)], method = "spearman")
            
          }
          
          
          stats.corr$r.p[s] = r$p.value
          stats.corr$r[s] = r$estimate
          stats.corr$n[s] = r$n
          if (!is.na(r$p.value)){
            if (r$p.value < 0.05) {
            significant[s] <- 1
            }
          }
        }
    }
  }
}
write.xlsx(stats.corr, file="resultsoutputNew.xls", sheetName = "deltabehcorr", append = T)


# scatter plots if any interesting correlation has been identified
#scatterplot(delta.cars~delta.ei, data=mydata[group==1,], smooth = FALSE, 
#            ylab = paste('delta-', colnames(bumedata)[31+2*i-1]),
#            xlab = paste('delta-', 'OFC:GABA/Glu+Gln'),
#            main = 'bumetanide')
         
        
# if anything interseting we can do the group comparison on these correlations: library(cocor)  

# permutation
set.seed(300)
perm.total <- 3000
g = group[1]
i = ROI.names[1]
j = MRS.trans[3]
s = 0
NOI <- paste0(i,j)
validity <- bumedata$QC.paired[,index[i,1]]
dat <- bumedata[validity==TRUE,]
perm.count <- matrix(rep(0, 17), nrow = 1, ncol = 17)
stats.corr0 <- stats.corr[c(50:66), "r"]
for (nperm in c(1:perm.total)){
  baseline.mrs <- dat[[NOI]][which(dat$stage4plot == "baseline" & dat$group4plot == g)]
  delta.mrs <- baseline.mrs - dat[[NOI]][which(dat$stage4plot == "3 months" & dat$group4plot == g)]
  perm.sample <- sample.int(length(delta.mrs))
  
  delta.mrs.perm <- delta.mrs[perm.sample]
  covariates <- cbind(dat$age, dat$ismale, dat$IQ)
  covariates <- covariates[which(dat$stage4plot == "baseline" & dat$group4plot == g),]
  stats.corr.perm <- perm.count
  s <- 0
  for  (b in beh[c(6:22)]){
    s = s + 1
    baseline.beh <- dat[[b]][which(dat$stage4plot=="baseline" & dat$group4plot == g)]
    delta.beh <- baseline.beh - dat[[b]][which(dat$stage4plot == "3 months" & dat$group4plot == g)]
    
    datanow <- cbind(baseline.mrs, delta.mrs.perm, baseline.beh, delta.beh, covariates)
    if (length(which(rowSums(is.na(datanow))>0)) > 0){
      datanow <- datanow[-which(rowSums(is.na(datanow))>0),]
    }
    
    r <- pcor.test(datanow[,2],datanow[,4], datanow[,c(1,3,5:7)], method = "spearman")
    stats.corr.perm[s] <- r$estimate
  }
  
  perm.count[c(1,2)] <- perm.count[c(1,2)] + (max(abs(stats.corr.perm[c(1,2)])) > abs(stats.corr0[c(1,2)]))
  perm.count[-c(1:2,16)] <- perm.count[-c(1:2,16)] + (max(abs(stats.corr.perm[-c(1:2,16)]),na.rm = TRUE) > abs(stats.corr0[-c(1:2,16)]))
  
}
write.xlsx(perm.count/3000, file="resultsoutputNew.xls", sheetName = "deltabehcorrperm", append = T)

library(FSA)   
bumedata.baseline <- bumedata[bumedata$stage==0,]
bumedata.3mont <- bumedata[bumedata$stage==3,]

# behaviour grouping at baseline in the bumetanide group
beh.index <- (bumedata.baseline$item_s3 > median(bumedata.baseline$item_s3) & bumedata.baseline$group == 1)
bumedata.3mont$group.stritify.beh <- bumedata.3mont$group
for (i in bumedata.baseline$ID[beh.index]){
  bumedata.3mont$group.stritify.beh[which(bumedata.3mont$ID == i)] = 2
}
bumedata.3mont$group.stritify.beh <- as.factor(bumedata.3mont$group.stritify.beh)
levels(bumedata.3mont$group.stritify.beh) <- c("control", "bumetanide.less", "bumetanide.more")

# MRS grouping at baseline in the bumetanide group
mrs.index <- (bumedata.baseline$INS.GABA.Glx.tissue.corrected > 
                median(bumedata.baseline$INS.GABA.Glx.tissue.corrected, na.rm = T) & 
                bumedata.baseline$group == 1)
bumedata.3mont$group.stritify.mrs <- bumedata.3mont$group
for (i in bumedata.baseline$ID[mrs.index]){
  bumedata.3mont$group.stritify.mrs[which(bumedata.3mont$ID == i)] = 2
}
bumedata.3mont$group.stritify.mrs <- as.factor(bumedata.3mont$group.stritify.mrs)
levels(bumedata.3mont$group.stritify.mrs) <- c("control", "bumetanide.lower", "bumetanide.higher")

table(bumedata.3mont$group.stritify.mrs, bumedata.3mont$group.stritify.beh)

bumedata.3mont$group.stritify <- matrix(rep(0,dim(bumedata.3mont)[1]),ncol=1, nrow=dim(bumedata.3mont)[1])
for (i in c(1:dim(bumedata.3mont)[1])){
  if (bumedata.3mont$group.stritify.mrs[i] == "control" && bumedata.3mont$group.stritify.beh[i] == "control"){
    bumedata.3mont$group.stritify[i] = "control"
  }
  if (bumedata.3mont$group.stritify.mrs[i] == "bumetanide.lower" && bumedata.3mont$group.stritify.beh[i] == "bumetanide.less"){
    bumedata.3mont$group.stritify[i] = "bumetanide1"
  } 
  if (bumedata.3mont$group.stritify.mrs[i] == "bumetanide.lower" && bumedata.3mont$group.stritify.beh[i] == "bumetanide.more"){
    bumedata.3mont$group.stritify[i] = "bumetanide2"
  } 
  if (bumedata.3mont$group.stritify.mrs[i] == "bumetanide.higher" && bumedata.3mont$group.stritify.beh[i] == "bumetanide.less"){
    bumedata.3mont$group.stritify[i] = "bumetanide3"
  } 
  if (bumedata.3mont$group.stritify.mrs[i] == "bumetanide.higher" && bumedata.3mont$group.stritify.beh[i] == "bumetanide.more"){
    bumedata.3mont$group.stritify[i] = "bumetanide4"
  }
}
bumedata.3mont$group.stritify = as.factor(bumedata.3mont$group.stritify)
#levels(bumedata.3mont$group.stritify) <- c("placebo", "bumetanide1", "bumetanide2")

table(bumedata.3mont$CRI_score, bumedata.3mont$group.stritify.mrs)
table(bumedata.3mont$CRI_index, bumedata.3mont$group.stritify.mrs)
table(bumedata.3mont$CRI_severity, bumedata.3mont$group.stritify.mrs)
table(bumedata.3mont$CRI_total, bumedata.3mont$group.stritify.mrs)
table(bumedata.3mont$CRI_side, bumedata.3mont$group.stritify.mrs)


# nonparametric test
stats.cri <- matrix(rep(0, 4), ncol = 4, nrow = 4)
rownames(stats.cri) <- colnames(bumedata.baseline)[c(54,55,124,125)]
colnames(stats.cri) <- c("degree of freedom", "Kruskal-Wallis chi-squared", "p-value", "fdr-p")
# wilcox.test(CRI_score~group, data=bumedata.3mont)  failed to converge
temp.index <- c(1,1,3,3)
s <- 1
for (i in c(54,55,124,125)){
  validity <- (bumedata.3mont$QC[,temp.index[s]]==1)
  kstest <- kruskal.test(bumedata.baseline[validity,i]~bumedata.3mont[validity,]$group.stritify.beh, data=bumedata.baseline[validity,])  # significant
  stats.cri[s,c(1:3)] <- c(kstest$parameter, kstest$statistic, kstest$p.value)
  s <- s + 1
}

stats.cri[,"fdr-p"] <- t(p.adjust(stats.cri[,"p-value"], method = "fdr"))
print(stats.cri)


# comparison between MRS patients and no-MRS patients in the control group
withmrs.index <- (rowSums(is.na(bumedata.baseline$QC[,c(1,3)]))==0 & bumedata.baseline$group==0) # control & MRS

bumedata.3mont$group.stritify.withmrs <- bumedata.3mont$group
for (i in bumedata.baseline$ID[which(withmrs.index)]){
  bumedata.3mont$group.stritify.withmrs[which(bumedata.3mont$ID == i)] = 3
}
bumedata.3mont$group.stritify.withmrs <- as.factor(bumedata.3mont$group.stritify.withmrs)
levels(bumedata.3mont$group.stritify.withmrs) <- c("control", "bumetanide", "control.mrs")

table(bumedata.3mont$group.stritify.withmrs)

# nonparametric test
stats.cri.mrs <- matrix(rep(0, 4), ncol = 4, nrow = 39)
rownames(stats.cri.mrs) <- colnames(bumedata.baseline)[c(5:43)]
colnames(stats.cri.mrs) <- c("degree of freedom", "Kruskal-Wallis chi-squared", "p-value", "fdr-p")
# wilcox.test(CRI_score~group, data=bumedata.3mont)  failed to converge
s <- 1
for (i in c(5:43)){
  kstest <- kruskal.test(bumedata.baseline[,i]~bumedata.3mont$group.stritify.mrs, data=bumedata.baseline)  # significant
  stats.cri.mrs[s,c(1:3)] <- c(kstest$parameter, kstest$statistic, kstest$p.value)
  s <- s + 1
}

stats.cri.mrs[,"fdr-p"] <- t(p.adjust(stats.cri.mrs[,"p-value"], method = "fdr"))
print(stats.cri.mrs)








