library(dplyr)
library(xmsPANDA)

setwd("C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range")

load(file = "HILIC_feature_combat.rda")
load(file = "HILIC_class.rda")
load(file = "HILIC_sample_link.rda")

#############################
# Parameters
#############################

outdata_type = "redisual"   # "raw", "transformed", "redisual"
pest = "OPs"   # "Pyrethroids", "Glyphosates", "Fungicides", "Neonicotinoids", "total"
outdata_loc = "C:/Users/QiYan/Dropbox/AIME/Pesticide/HILICpos/HILIC_OPs/PANDA_input/"
online = FALSE

setwd("C:/Users/QiYan/Dropbox/AIME/Pesticide/HILICpos/HILIC_OPs/PANDA_input/")

#############################
# Total/Each
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
  
  single_sampleID <- sampleID
  single_sampleID <- single_sampleID[order(single_sampleID$SampleID),]
  colnames(single_sampleID) <- c("SampleID","factorcase","sex","birthyear","maternal_age","maternal_raceeth","maternal_edu"
                                 ,"lengthgestation","pregcompl","ttcbl","preterm","usborn")
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
  
  single_sampleID <- sampleID
  single_sampleID <- single_sampleID[order(single_sampleID$SampleID),]
  colnames(single_sampleID) <- c("SampleID","factorcase","sex","birthyear","maternal_age","maternal_raceeth","maternal_edu"
                          ,"lengthgestation","pregcompl","ttcbl","preterm","usborn")
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

#############################
# Raw/Transform/Residual
#############################

{

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
  
  single_sampleID$SampleID <- row.names(after.prepro.feature)
  
  final_sample <- single_sampleID
  final_linkid <- after.prepro.linkid
  final_feature <- after.prepro.feature

}else{
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
  
  single_sampleID$SampleID <- row.names(after.prepro.feature)
  
  if(outdata_type=="transformed"){
    final_sample <- single_sampleID
    final_linkid <- after.prepro.linkid
    final_feature <- after.prepro.feature
  }else{
    
    # Recode some covariates due tp spase data
    single_sampleID$maternal_age[which(single_sampleID$maternal_age==5)] <- 4
    single_sampleID$maternal_raceeth[which(single_sampleID$maternal_raceeth==4)] <- 3
    single_sampleID$maternal_raceeth[which(single_sampleID$maternal_raceeth==5)] <- 3
    single_sampleID$maternal_edu[which(single_sampleID$maternal_edu==2)] <- 1
    single_sampleID$maternal_edu[which(single_sampleID$maternal_edu==3)] <- 2
    single_sampleID$maternal_edu[which(single_sampleID$maternal_edu==4)] <- 3
    single_sampleID$maternal_edu[which(single_sampleID$maternal_edu==5)] <- 4
    
    feature_w_cov <- cbind(single_sampleID,after.prepro.feature)
    
    # adjust for covariates on each metabolits and then get residuals
    fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,13:ncol(feature_w_cov)]) ~ as.factor(maternal_age)+as.factor(maternal_edu)+
                        as.factor(maternal_raceeth), na.action = na.exclude) # residual_3
    residual_feature <- as.matrix(residuals(fit_feature),nrow = dim(feature_w_cov)[1],ncol = dim(save_feature)[1])
    save_residual <- as.data.frame(residual_feature)
    save_residual <- cbind(feature_w_cov$SampleID,save_residual)
    save_residual <- save_residual[order(save_residual$`feature_w_cov$SampleID`),]
    row.names(save_residual) <- c(1:nrow(save_residual))
    
    row.names(save_residual) <- save_residual$`feature_w_cov$SampleID`
    save_residual <- save_residual[,-1]
    save_residual<-na.omit(save_residual)
    
    complete_sub <- row.names(save_residual)
    single_sampleID <- subset(single_sampleID, single_sampleID$SampleID %in% complete_sub)
    
    final_sample <- single_sampleID
    final_linkid <- after.prepro.linkid
    final_feature <- save_residual
    
  }
}

}

#############################
# Online/Offline
#############################
{
  
if(online){

  label <- final_sample$factorcase
  sample <- final_sample$SampleID
  final_feature <- cbind(sample,label,final_feature)
  
  save(final_feature,final_linkid,final_sample, file = paste("HILIC_Online_OPs_",outdata_type,".RData",sep = ""))
  write.table(final_feature,file=paste("HILIC_Online_OPs_",outdata_type,".txt",sep = ""),sep = "\t",row.names = F,quote = F)

}else{
  save(final_feature,final_linkid,final_sample, file = paste("HILIC_OPs_",outdata_type,".RData",sep = ""))
}
  
}
