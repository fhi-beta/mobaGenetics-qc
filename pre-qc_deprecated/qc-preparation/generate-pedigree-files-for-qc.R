# File preparation steps for files needed to run QC

library(dplyr)
options(scipen=999)

# load connection files from Bergen and Rotterdam
con.bergen <- read.table('/home/oyvind/hunt-cloud/mnt/archive/ROTTERDAM1/connection/connection_rotterdam1.txt.csv', header = T, stringsAsFactors = F, sep='\t')
con.rotterdam <- read.table('/home/oyvind/hunt-cloud/mnt/archive/ROTTERDAM1/qc/GSA2016_83_025_STS.csv', header = T, skip = 9, sep=',', stringsAsFactors = F)
fam <- read.table('/home/oyvind/hunt-cloud/mnt/archive/ROTTERDAM1/data/plink/PLINK_260318_0459/bedset/GSA2016_MOBA_025.fam', header = F, stringsAsFactors = F)

bergen.sub <- con.bergen %>% select(RetrievalDetail_ID, FamID, NumInTriade, Role, Gender)
rotterdam.sub <- con.rotterdam %>% select(Sample_ID, ArrayPicker)
fam.sub <- fam %>% select(V1, V2)
names(fam.sub) <- c('famfid','famiid')

# merge con info from bergen and rotterdam and then the famids present in the initial bedset
con.total.tmp <- merge(bergen.sub, rotterdam.sub, by.x='RetrievalDetail_ID', by.y='Sample_ID')
con.total <- merge(con.total.tmp, fam.sub, by.x='RetrievalDetail_ID', by.y='famiid')

#### generate plink file for updating family and sample id
recode.ids <- con.total %>% select(famfid, RetrievalDetail_ID, FamID, ArrayPicker)

# remove "ER00" prefix from famid
recode.ids$FamID <- gsub('ER0*', '',recode.ids$FamID)

#### generate plink file for updating parents
# separate roles
mother <- con.total %>% filter(Role=='Mother') %>% select(FamID, ArrayPicker) %>% rename(sentrixid_mother=ArrayPicker)
father <- con.total %>% filter(Role=='Father') %>% select(FamID, ArrayPicker) %>% rename(sentrixid_father=ArrayPicker)
child  <- con.total %>% filter(Role=='Child') %>% select(FamID, ArrayPicker) %>% rename(sentrixid_child=ArrayPicker)

# join on FamID
recode.parents = full_join(child,father, by='FamID') %>% full_join(mother, by='FamID')

# convert NAs to 0
recode.parents$sentrixid_child[is.na(recode.parents$sentrixid_child)] = 0
recode.parents$sentrixid_mother[is.na(recode.parents$sentrixid_mother)] = 0
recode.parents$sentrixid_father[is.na(recode.parents$sentrixid_father)] = 0
# tmp[is.na(tmp)] = 0  <-- possible but not ideal

# remove "ER00" prefix from famid
recode.parents$FamID <- gsub('ER0*', '',recode.parents$FamID)

#### generate file for updating gender
recode.sex <- con.total %>% select(FamID, ArrayPicker, Gender)

# recode gender to a plink accepted format
recode.sex[recode.sex$Gender=='X_MALE',]$Gender <- 'M'
recode.sex[recode.sex$Gender=='X_FEMALE',]$Gender <- 'F'

# remove "ER00" prefix from famid
recode.sex$FamID <- gsub('ER0*', '',recode.sex$FamID)

# One sample with pregID 300313 was not delivered in connection files from
# Bergen. In order to complete the recode lists the SentrixID for this sample
# was found by finding the remaining sentrixID not assigned from 
missing.sample <- rotterdam.sub[!rotterdam.sub$ArrayPicker %in% recode.ids$ArrayPicker,]

# Add missing sample to recode.ids
add.to.recode.ids <- data.frame(famfid=c(11161), RetrievalDetail_ID=c(300313), FamID=c(0), ArrayPicker=c('201904770001_R08C01'))
recode.ids <- rbind(recode.ids, add.to.recode.ids)

# Specify gender for missing sample (verified female by Jonas downstream)
add.to.recode.sex <- data.frame(FamID=c(0), ArrayPicker=c('201904770001_R08C01'), Gender=c('F'))
recode.sex <- rbind(recode.sex, add.to.recode.sex)

# write files
#out.path = c('/home/oyvind/hunt-cloud/mnt/archive/ROTTERDAM1/recode-files/')
out.path = c('/home/oyvind/harvest/media/local-disk2/helgeland/rotterdam1/recode-files/')

# OBS: Files are set to immutable (sudo chattr +i) at the server
write.table(recode.ids, file = file.path(out.path, 'recode-ids-rotterdam1.txt'), col.names = F, row.names = F, sep=' ', quote = F)
write.table(recode.parents, file = file.path(out.path, 'recode-parents-rotterdam1.txt'), col.names = F, row.names = F, sep=' ', quote = F)
write.table(recode.sex, file = file.path(out.path, 'recode-sex-rotterdam1.txt'), col.names = F, row.names = F, sep=' ', quote = F)

