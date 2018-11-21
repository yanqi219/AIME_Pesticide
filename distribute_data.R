library(dplyr)
library(xmsPANDA)

setwd("C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range")

load(file = "HILIC_feature_combat.rda")
load(file = "HILIC_class.rda")
load(file = "HILIC_sample_link.rda")

#############################
# Parameters
#############################

outdata_type = "residule"   # "raw", "transformed", "redisual"
pest = "OPs"   # "Pyrethroids", "Glyphosates", "Fungicides", "Neonicotinoids", "total"
outdata_loc = "C:/Users/QiYan/Dropbox/AIME/Pesticide/HILICpos/HILIC_OPs/PANDA_input/"
online = FALSE

#############################
#
#############################
{
# Get sample id & Covs

if(pest=="total"){
  sampleID <- HILIC_class[HILIC_class$casecontrol..Factor.1.=="0",]
  sampleID <- sampleID[sampleID$Exposure_category..Factor2.=="A"|sampleID$Exposure_category..Factor2.=="B"|
                         sampleID$Exposure_category..Factor2.=="C"|sampleID$Exposure_category..Factor2.=="D"|
                         sampleID$Exposure_category..Factor2.=="E"|sampleID$Exposure_category..Factor2.=="H",c(1,11:21)]
  sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
  sampleID <- sampleID[!is.na(sampleID$SampleID),]
  
  online_sampleID <- sampleID
  online_sampleID <- online_sampleID[order(online_sampleID$SampleID),]
} else if(pest == "OPs"){
  pest.code = "A"
  } else if(pest == "Pyrethroids"){
    pest.code = "B"
    } else if(pest == "Glyphosates"){
      pest.code = "C"
      } else if(pest == "Fungicides"){
        pest.code = "D"
        } else if(pest == "Neonicotinoids"){
          pest.code = "E"
        }

if(pest!="total"){
  sampleID <- HILIC_class[HILIC_class$casecontrol..Factor.1.=="0",]
  sampleID <- sampleID[sampleID$Exposure_category..Factor2.==pest.code|sampleID$Exposure_category..Factor2.=="H",c(1,11:21)]
  sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
  sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.==pest.code] <- "Exposed"
  sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="H"] <- "Unexposed"
  sampleID <- sampleID[!is.na(sampleID$SampleID),]
  
  online_sampleID <- sampleID
  online_sampleID <- online_sampleID[order(online_sampleID$SampleID),]
}

#get subjects' 3 replicates sample ID

rownumber <- nrow(sampleID)
sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
rownames(sampleID) <- c(1:nrow(sampleID))

tail <- c("_1","_3","_5")
tail <- rep(tail,times=rownumber)
V4 <- paste(sampleID$SampleID,tail,sep = "")

sampleID <- cbind(V4,sampleID)

sampleID <- sampleID[,-2]

colnames(sampleID) <- c("SampleID","factorcase","sex","birthyear","maternal_age","maternal_raceeth","maternal_edu"
                                             ,"lengthgestation","pregcompl","ttcbl","preterm","usborn")

#get the feature table

link <- sampleID[,c(1,2)]
colnames(link) <- c("Sample.ID","factorcase")
link <- merge(link,HILIC_sample_link,by="Sample.ID")
link <- link[,1:3]
link$V4 <- ".mzXML"
link$File.Name <- paste(link$File.Name,link$V4,sep = "")
addon <- data.frame(Sample.ID=c("mz","time"),factorcase=c("mz","time"),File.Name=c("mz","time"),V4=c("mz","time"))
link <- rbind(addon,link)

featurename <- as.vector(link$File.Name)
twomore <- c("mz","time")
featurename <- c(twomore,featurename)

feature <- HILIC_feature_combat[,colnames(HILIC_feature_combat) %in% featurename]    ##keep only wanted columns

colname_feature <- colnames(feature)

long_feature <- as.data.frame(colname_feature)
colnames(long_feature) <- "File.Name"
long_feature <- merge(long_feature,link,by="File.Name")

class <- long_feature[3:nrow(long_feature),2:3]
colnames(class) <- c("SampleID","factorcase")

long_feature <- long_feature[,1:2]
wide_feature <- as.data.frame(t(long_feature))

names(feature) <- lapply(wide_feature[2, ], as.character)

rownames(feature) <- paste("met_",c(1:nrow(feature)),sep="")

temp <- feature[,c(1,2)]
feature <- data.frame(t(feature[,-c(1,2)]))
feature <- feature[order(rownames(feature)),]
feature <- data.frame(t(feature))
feature <- cbind(temp,feature)   # feature is the input for downstream analysis
}

if(online){
  ready_for_regression<-data_preprocess(Xmat=feature,Ymat=NA,feature_table_file=NA,parentoutput_dir=outdata_loc,class_labels_file=NA,num_replicates=3,feat.filt.thresh=NA,
                                        summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0,group.missing.thresh=0,
                                        log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
                                        samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="zeros",featselmethod=NA)
  
  after.prepro.feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
  met <- rownames(after.prepro.feature)
  label <- online_sampleID$Exposure_category..Factor2.
  online_linkid <- cbind(met,after.prepro.feature[,c(1,2)])
  online_feature <- data.frame(t(after.prepro.feature[,-c(1,2)]))
  sample <- as.character(online_sampleID$SampleID)
  online_feature <- cbind(sample,label,online_feature)
  
  ################### save file here
}

# Start pre-processing

class_labels_file <- sampleID[,c(1,2)]

if(outdata_type=="raw"){
  ready_for_regression<-data_preprocess(Xmat=feature,Ymat=class_labels_file,feature_table_file=NA,parentoutput_dir=outdata_loc,class_labels_file=NA,num_replicates=3,feat.filt.thresh=NA,
                                        summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.5,group.missing.thresh=0.8,
                                        log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
                                        samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="halfdatamin",featselmethod=NA)

  after.prepro.feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
  na_count <-sapply(after.prepro.feature, function(y) sum(is.na(y)))
  summary(na_count)
  
  row.names(after.prepro.feature) <- c(paste("met_",1:nrow(after.prepro.feature),sep = ""))
  after.prepro.linkid <- after.prepro.feature[,1:2]
  after.prepro.feature <- t(after.prepro.feature[,-c(1:2)])
  after.prepro.feature <- after.prepro.feature[order(row.names(after.prepro.feature)), ]
  
  ##################save file here
  
} else{
  ready_for_regression<-data_preprocess(Xmat=feature,Ymat=class_labels_file,feature_table_file=NA,parentoutput_dir=outdata_loc,class_labels_file=NA,num_replicates=3,feat.filt.thresh=NA,
                                        summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.5,group.missing.thresh=0.8,
                                        log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
                                        samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="halfdatamin",featselmethod=NA)
  
  after.prepro.feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
  na_count <-sapply(after.prepro.feature, function(y) sum(is.na(y)))
  summary(na_count)
  
  row.names(after.prepro.feature) <- c(paste("met_",1:nrow(after.prepro.feature),sep = ""))
  after.prepro.linkid <- after.prepro.feature[,1:2]
  after.prepro.feature <- t(after.prepro.feature[,-c(1:2)])
  after.prepro.feature <- after.prepro.feature[order(row.names(after.prepro.feature)), ]
  
  
}
  










#############################
#Pesticides exposed controls and unexposed controls
#############################
split.pest <- function(pest){
  
  if(pest == "OPs"){
    pest.code = "A"
  } else if(pest == "Pyrethroids"){
    pest.code = "B"
  } else if(pest == "Glyphosates"){
    pest.code = "C"
  } else if(pest == "Fungicides"){
    pest.code = "D"
  } else if(pest == "Neonicotinoids"){
    pest.code = "E"
  } 
  
  ##extract cases and controls ID from class dataset
  sampleID <- HILIC_class[HILIC_class$casecontrol..Factor.1.=="0",]
  sampleID <- sampleID[sampleID$Exposure_category..Factor2.==pest.code|sampleID$Exposure_category..Factor2.=="H",c(1,11)]
  colnames(sampleID)
  sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
  sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.==pest.code] <- "Exposed"
  sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="H"] <- "Unexposed"
  sampleID <- sampleID[!is.na(sampleID$SampleID),]
  
  ##get subjects' 3 replicates sample ID
  
  rownumber <- nrow(sampleID)
  sampleID <- sampleID[rep(row.names(sampleID),3),]
  sampleID <- sampleID[order(sampleID$SampleID),]
  rownames(sampleID) <- c(1:nrow(sampleID))
  
  tail <- data.frame(V3=c("_1","_3","_5"))
  tail <- data.frame(rep(tail$V3,times=rownumber))
  colnames(tail) <- "V3"
  
  sampleID <- cbind(sampleID,tail)
  sampleID$V4 <- paste(sampleID$SampleID,sampleID$V3,sep = "")
  
  sampleID <- sampleID[,c("V4","Exposure_category..Factor2.")]
  
  colnames(sampleID) <- c("SampleID","factorcase")
  
  ##get the feature table
  
  link <- sampleID
  colnames(link) <- c("Sample.ID","factorcase")
  link <- merge(link,HILIC_sample_link,by="Sample.ID")
  link <- link[,1:3]
  link$V4 <- ".mzXML"
  link$File.Name <- paste(link$File.Name,link$V4,sep = "")
  addon <- data.frame(Sample.ID=c("mz","time"),factorcase=c("mz","time"),File.Name=c("mz","time"),V4=c("mz","time"))
  link <- rbind(addon,link)
  
  featurename <- as.vector(link$File.Name)
  twomore <- c("mz","time")
  featurename <- c(twomore,featurename)
  
  feature <- HILIC_feature_combat[,colnames(HILIC_feature_combat) %in% featurename]    ##keep only wanted columns
  
  colname_feature <- colnames(feature)
  
  long_feature <- as.data.frame(colname_feature)
  colnames(long_feature) <- "File.Name"
  long_feature <- merge(long_feature,link,by="File.Name")
  
  class <- long_feature[3:nrow(long_feature),2:3]
  colnames(class) <- c("SampleID","factorcase")
  
  long_feature <- long_feature[,1:2]
  wide_feature <- as.data.frame(t(long_feature))
  
  names(feature) <- lapply(wide_feature[2, ], as.character)
  
  ##save data file
  
  file.dir = paste("C:/Users/QiYan/Dropbox/AIME/Pesticide/HILICpos/HILIC_",pest,"/PANDA_input", sep = "")
  setwd(file.dir)
  
  write.table(class,file="HILIC_classlabels_control_expo_unexpo.txt",sep="\t",row.names = F,quote = F)
  write.table(feature,file="HILIC_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)
  
}

split.pest(pest = "OPs")
split.pest(pest = "Pyrethroids")
split.pest(pest = "Glyphosates")
split.pest(pest = "Fungicides")
split.pest(pest = "Neonicotinoids")

#############################
#All
#############################

##extract cases and controls ID from class dataset
sampleID <- HILIC_class[HILIC_class$Exposure_category..Factor2.=="G"|HILIC_class$Exposure_category..Factor2.=="H",c(1,10:21)]
colnames(sampleID)
sampleID$casecontrol..Factor.1.[sampleID$casecontrol..Factor.1.==1] <- "case"
sampleID$casecontrol..Factor.1.[sampleID$casecontrol..Factor.1.==0] <- "control"
sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="G"] <- "Exposed"
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="H"] <- "Unexposed"
sampleID <- sampleID[!is.na(sampleID$SampleID),]

##get subjects' 3 replicates sample ID

sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
rownames(sampleID) <- c(1:1011)

tail <- data.frame(V3=c("_1","_3","_5"))
tail <- data.frame(rep(tail$V3,times=337))
colnames(tail) <- "V3"

sampleID <- cbind(sampleID,tail)
sampleID$V4 <- paste(sampleID$SampleID,sampleID$V3,sep = "")

sampleID <- subset(sampleID,select=-c(V3,SampleID))
sampleID <- plyr::rename(sampleID,c('V4'='SampleID'))
colnames(sampleID)
sampleID <- sampleID[c("SampleID","casecontrol..Factor.1.","Exposure_category..Factor2.","sex","birthyear","maternal_age","maternal_raceeth","maternal_edu"
                       ,"lengthgestation","pregcompl","ttcbl","preterm","usborn")]
##get the feature table

link <- sampleID[,1:3]
colnames(link) <- c("Sample.ID","casecontrol..Factor.1.","Exposure_category..Factor2.")
link <- merge(link,HILIC_sample_link,by="Sample.ID")
link <- link[,1:4]
link$V4 <- ".mzXML"
link$File.Name <- paste(link$File.Name,link$V4,sep = "")
addon <- data.frame(Sample.ID=c("mz","time"),casecontrol..Factor.1.=c("mz","time"),Exposure_category..Factor2.=c("mz","time"),File.Name=c("mz","time"),V4=c("mz","time"))
link <- rbind(addon,link)

featurename <- as.vector(link$File.Name)
twomore <- c("mz","time")
featurename <- c(twomore,featurename)

# names(HILIC_feature_combat) <- lapply(HILIC_feature_combat[1, ], as.character)

feature <- HILIC_feature_combat[,colnames(HILIC_feature_combat) %in% featurename]    ##keep only wanted columns

colname_feature <- colnames(feature)

long_feature <- as.data.frame(colname_feature)
colnames(long_feature) <- "File.Name"
long_feature <- merge(long_feature,link,by="File.Name")

long_feature <- long_feature[,1:2]
wide_feature <- as.data.frame(t(long_feature))

names(feature) <- lapply(wide_feature[2, ], as.character)

transpose <- feature[,3:ncol(feature)]
transpose <- t(transpose)
transpose <- transpose[order(row.names(transpose)),]
transpose <- t(transpose)

feature <- cbind(feature[,1:2],transpose)

sample.label <- sampleID[,1:3]
sample.label$class <- paste(sample.label$casecontrol..Factor.1.,sample.label$Exposure_category..Factor2.,sep = "")
sample.label <- sample.label[,c(1,4)]

##save data file

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos")

write.table(sampleID,file="HILIC_classlabels_all.txt",sep="\t",row.names = F,quote = F)
write.table(sample.label,file="HILIC_classlabels_for_panda_all.txt",sep="\t",row.names = F,quote = F)
write.table(feature,file="HILIC_ftrsmzcalib_combat_ordered_all.txt",sep = "\t",row.names = F,quote = F)

#############################
#Get residual matrix for all
#############################

library(dplyr)
library(plyr)
library(xmsPANDA)
library(mixOmics)

class <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_classlabels_for_panda_all.txt"
feature <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_ftrsmzcalib_combat_ordered_all.txt"
outloc <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos"

ready_for_regression<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file=feature,parentoutput_dir=outloc,class_labels_file=class,num_replicates=3,feat.filt.thresh=NA,
                                      summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.5,group.missing.thresh=0.8,
                                      log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
                                      samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="halfdatamin",featselmethod=NA)

feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
na_count <-sapply(feature, function(y) sum(is.na(y)))
summary(na_count)

row.names(feature) <- c(paste("met_",1:nrow(feature),sep = ""))
after.prepro.linkid <- feature[,1:2]
after.prepro.feature <- t(feature[,-c(1:2)])
after.prepro.feature <- after.prepro.feature[order(row.names(after.prepro.feature)), ]

# get the class table
class<-read.table("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_classlabels_all.txt",header = T)
class$flag <- stringr::str_sub(class$SampleID,-2,-1)

class <- class[-which(class$flag!="_1"),-14]

# reformat feature table from wide to long, create unique metabolite id number, link with covariates table
save_feature <- as.data.frame(after.prepro.feature)
save_feature$SampleID <- row.names(save_feature)

# merge covariates data wirh feature data
feature_w_cov <- merge(class,save_feature, by = "SampleID")

fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,14:ncol(feature_w_cov)]) ~ as.factor(sex)+as.factor(maternal_age)+as.factor(maternal_edu)+as.factor(pregcompl)
                  +ttcbl+as.factor(maternal_raceeth), na.action = na.exclude)
#
residual_feature <- as.matrix(residuals(fit_feature),nrow = dim(feature_w_cov)[1],ncol = dim(save_feature)[1])
save_residual <- as.data.frame(residual_feature)
save_residual <- cbind(feature_w_cov$SampleID,save_residual)
save_residual <- save_residual[order(save_residual$`feature_w_cov$SampleID`),]
row.names(save_residual) <- c(1:nrow(save_residual))

# reformat residual from long to wide
n <- save_residual$`feature_w_cov$SampleID`

wide_save_residual <- as.data.frame(t(save_residual[,-1]))
colnames(wide_save_residual) <- n
id <- row.names(wide_save_residual)
wide_save_residual <- cbind(id,wide_save_residual)
names(wide_save_residual)[names(wide_save_residual) == 'id'] <- 'metablite'
row.names(wide_save_residual) <- c(1:nrow(wide_save_residual))

# link wide_save_residual with save_link, retrive m/z and time back
wide_save_residual <- cbind(after.prepro.linkid,wide_save_residual[,-1])
wide_save_residual <- wide_save_residual[order(wide_save_residual$mz,wide_save_residual$time),]
row.names(wide_save_residual) <- c(1:nrow(wide_save_residual))

# # replace na with 0
# wide_save_residual<-replace(wide_save_residual,is.na(wide_save_residual),0)

# remove na
wide_save_residual<-wide_save_residual[sapply(wide_save_residual, function(x) !any(is.na(x)))]
complete_sub <- row.names(t(wide_save_residual[,-c(1:2)]))
class <- subset(class, class$SampleID %in% complete_sub)

##save data file

save(class, wide_save_residual, file = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_all_residual.RData") ## for WGCNA