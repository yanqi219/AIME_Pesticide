
library("xlsx")

setwd("C:/Users/Qi/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range")

HILIC_feature_combat<-read.table("HILICpos_mzcalibrated_featuretable_COMBAT.txt",header = T)
HILIC_feature_raw<-read.table("HILICpos_mzcalibrated_featuretable_RAW.txt",header = T)
HILIC_class<-read.xlsx("C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/AIME_to_Emory_586_111216_KU.xlsx", 1)
HILIC_sample_link<-read.table("PNS_sample_ID_mapping_HILICpos.txt",header = T)
HILIC_match_kegg<-read.table("DBmatches_KEGG.txt",header = T)
HILIC_match_hmdb<-read.table("DBmatches_HMDB.txt",header = T)
HILIC_match_lipidmaps<-read.table("DBmatches_LipidMaps.txt",header = T)
HILIC_replicates_corr<-read.table("Pairwise_Pearson_correlation_technical_replicates.txt",header = T)

save(HILIC_feature_combat, file = "HILIC_feature_combat.rda")
save(HILIC_feature_raw, file = "HILIC_feature_raw.rda")
save(HILIC_class, file = "HILIC_class.rda")
save(HILIC_sample_link, file = "HILIC_sample_link.rda")
save(HILIC_match_kegg, file = "HILIC_match_kegg.rda")
save(HILIC_match_hmdb, file = "HILIC_match_hmdb.rda")
save(HILIC_match_lipidmaps, file = "HILIC_match_lipidmaps.rda")
save(HILIC_replicates_corr, file = "HILIC_replicates_corr.rda")