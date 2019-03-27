### MAJOR-FUCKUP DETECTOR
### created for MoBa-Rotterdam2 QC
### date: 20190326-0330
### Jonas B.
#local: setwd("/Users/jb/Biostuff/ROTTERDAM2/")

library(dplyr)

####
#### arguments
####

args = commandArgs(trailingOnly=TRUE)

# variable inputs
fam_file = args[1] # ".fam"
ibd_file = args[2] # ".genome"
sex_file = args[3] # ".sex"
# static (permanent) inputs
del_file = args[4] # sample names to be ignored, one per line
rmp_file = args[5] # a-priori reconstruction file with columns FID IID PID MID SEX
# outputs
rep_fileRdat = args[6] # fckpclusters_20190xxx.RData 
upd_fam_file = args[7]
flg_fam_file = args[8]

### LOCAL ALTERNATIVE FOR DEBUGGING
# variable inputs
# fam_file = "fam.fam"
# ibd_file = "ibd.txt"
# sex_file = "sex.sex"
# # static (permanent) inputs
# del_file ="definitely_bad_samples.txt"
# rmp_file = "reconstructions.txt"
# # outputs
# rep_fileRdat = "experimental_fckpclusters_20190326_4.RData"
# upd_fam_file = "deleteme_updated_fam_file.txt"
# flg_fam_file = "deleteme_flag_file.txt"


####
#### load the files
####

fam = read.table(fam_file,h=F,stringsAsFactors = F) # head(fam); dim(fam)
if (any(!c("V1","V2","V3","V4","V5")  %in% colnames(fam))) warning("the needed fam-file columns do not exist!")
fam_orig = fam # save to RData for problemsolving if needed


# this could be unflagged
# CONS: in some cases detection of bad sample (a pair of undeclared relatedness) will be impossible
fam$V3[which(!fam$V3 %in% fam$V2)] = 0
fam$V4[which(!fam$V4 %in% fam$V2)] = 0




full = read.table(ibd_file,h=T,stringsAsFactors = F)# head(full)
if (any(!c("IID1","IID2","Z1","PI_HAT")  %in% colnames(full))) warning("the needed IBD-file columns do not exist!")
full = full[,c("IID1","IID2","Z1","PI_HAT")]

sex = read.table(sex_file,h=T,stringsAsFactors = F)
if (any(!c("FID","IID","PEDSEX","SNPSEX","STATUS","F","YCOUNT")  %in% colnames(sex))) warning("the needed sex-file columns do not exist!")
if (any(sex$IID!=fam$V2)) warning("the '.sexcheck' file IDs do not match '.fam' file IDs!")

# samples with clearly wrong DNA
del = read.table(del_file,stringsAsFactors = F)

# TEMPORARY MANUAL REMOVAL OF ONE VERY BAD CHIP *****
#chip_ids = unlist(lapply(fam$V2,function(x) unlist(strsplit(x,"_"))[1]))
#bad_samps = fam$V2[which(chip_ids  %in% c("203060680096","203060680146"))]
#del = data.frame(V1 = bad_samps,stringsAsFactors = F)

# blueprint of how to rearrange some families (index = IID)
upd = read.table(rmp_file,h=T,stringsAsFactors = F,sep="\t")


####
####  FLAG FILE
####

flg = fam[,c("V2","V5")]
colnames(flg) = c("IID","phenoOK")
flg$phenoOK = TRUE  # default
# will be continuously updated further
        

####
#### define thresholds to determine SEX problems
####

fem_y = median(sex$YCOUNT[which(sex$PEDSEX==2)]) # Yc center of female cluster
mal_y = median(sex$YCOUNT[which(sex$PEDSEX==1)]) # Yc center of   male cluster
fem_f = median(sex$F[which(sex$PEDSEX==2)]) # F center of female cluster
mal_f = median(sex$F[which(sex$PEDSEX==1)]) # F center of   male cluster
y_thr = fem_y + (mal_y-fem_y)*0.666 # threshold for Ycounts
f_thr = fem_f + (mal_f-fem_f)*0.666 # threshold for Fval (Xchr)
rm(fem_y,mal_y,fem_f,mal_f)



####
####  FAM FILE UPDATOR
####

## report problematic chips:
tbl = table(unlist(lapply(del$V1,function(x) unlist(strsplit(x,"_"))[1])))
problm_chips = data.frame(Array_Name=names(tbl),Problem_Count=as.numeric(tbl),stringsAsFactors = F)
problm_chips = problm_chips[order(problm_chips$Problem_Count,decreasing = T),]
#barplot(problm_chips$Problem_Count) # -> to report .Rdata

### 1) mark truly wrong samples for "deletion"

# update FID (aka family ID or PregID)
rix = which(fam$V2 %in% del$V1)
if (length(rix)>0) {

new_fids = unlist(lapply(seq(length(rix)),function(x) paste("prblm",paste(rep("0",3-nchar(x)),collapse=""),x,sep="")))
fam$V1[rix] = new_fids; rm(new_fids)

# FLAG FILE UPDATE: for wrong samples
flg$phenoOK[rix] = FALSE

# update SEX of these truly wrong samples
fam$V5[rix] = 0 # default, for unresolved cases
#hist(sex$YCOUNT,breaks=100,col="grey"); abline(v=y_thr,col="red")
fam$V5[rix][which((sex$YCOUNT[rix]>y_thr)&(sex$F[rix]>f_thr))] = 1 # genetic males
fam$V5[rix][which((sex$YCOUNT[rix]<y_thr)&(sex$F[rix]<f_thr))] = 2 # genetic females

# reset parents of these samples
fam$V3[rix] = 0
fam$V4[rix] = 0

}

# delete these samples from parental columns
fam$V3[which(fam$V3 %in% del$V1)] = 0
fam$V4[which(fam$V4 %in% del$V1)] = 0


### 2)   update relationships

if (nrow(upd)>0) {
for (i in 1:nrow(upd)) {
        rix = which(fam$V2==upd$IID[i])
        if (length(rix)==1) {
                fam$V1[rix] = upd$FID[i]
                fam$V3[rix] = upd$PID[i]
                fam$V4[rix] = upd$MID[i]
                fam$V5[rix] = upd$SEX[i] # sex update is also done!
        } else {
                if (length(rix)>1) {
                        warning("some lines in '.fam' file won't be updated!")
                        cat("\n reason: this particular IID in 'permanent_reconstruct-fam.txt' was found multiple times in '.fam': \n")
                        print(upd[i,])
                }
        }
        rm(rix)
}
} # end of IF

...  think about this & or |  at Fthr ...

# FLAG FILE UPDATE: remaining families (after hard-coded family rearrangements)
bad_boys_ix = which((fam$V5==1)&((sex$YCOUNT<y_thr)&(sex$F<f_thr)))  # declared males not males
bad_girl_ix = which((fam$V5==2)&((sex$YCOUNT>y_thr)&(sex$F>f_thr)))  # declared females not females
bad_indexes = unique(c(bad_boys_ix,bad_girl_ix))
#plot(y=sex$YCOUNT, x=sex$F)
#points(y=sex$YCOUNT[bad_boys_ix], x=sex$F[bad_boys_ix],pch=19,col="blue")
#points(y=sex$YCOUNT[bad_girl_ix], x=sex$F[bad_girl_ix],pch=19,col="red")
flg$phenoOK[bad_indexes] = FALSE
flg$phenoOK[bad_indexes] = FALSE 

# update SEX in fam file
fam$V5[bad_boys_ix] = 2
fam$V5[bad_girl_ix] = 1


# ...  this is problematic in fertile X0 kariotype females !!!  two instances found .. 
# ... not good for phasing ... 


# check whether updated sex in fam file does not conflict with parental status...
if(any(fam$V2[which(fam$V5==2)] %in% fam$V3)) warning("genetic females detected in V3 (dad's column)!")
if(any(fam$V2[which(fam$V5==1)] %in% fam$V4)) warning("genetic males detected in V4 (mom's column)!")


####
####  SEX FILE UPDATOR (used only for plotting)
####

sex_upd = sex
sex_upd$PEDSEX[which(sex_upd$IID %in% del$V1)] = 0 # only to remove from plotting
for (i in 1:nrow(upd)) {
        rix = which(sex_upd$IID==upd$IID[i])
        if (length(rix)==1) {
        sex_upd$PEDSEX[rix] = upd$SEX[i]
        } else {
                if (length(rix)>1) {
                        warning("some lines in '.sexcheck' file won't be updated!")
                        cat("\n reason: this particular IID in 'permanent_reconstruct-fam.txt' was found multiple times in '.fam': \n")
                        print(upd[i,])
                }
        }
        rm(rix)
}




############
############ DATA PREP
############

# what role individual plays (father, mother, child)
dad_ids =  fam$V3; dad_ids = sort(unique(dad_ids[which(dad_ids %in% fam$V2)])) # !***
mom_ids =  fam$V4; mom_ids = sort(unique(mom_ids[which(mom_ids %in% fam$V2)])) # !***

# what is genetic SEX of the individual
Yfemales = sort(unique(sex$IID[which(sex$YCOUNT< y_thr)]))
Ymales   = sort(unique(sex$IID[which(sex$YCOUNT>=y_thr)]))

# quick check for very likely pedigree mistakes:
#sum(dad_ids %in% Yfemales)
#sum(mom_ids %in% Ymales)
#sum(!mom_ids %in% c(Ymales,Yfemales))
#sum(!dad_ids %in% c(Ymales,Yfemales))

### IID-FID key (add-on for later)
sub2 = fam[,c("V2","V1")]
sub3 = fam[which(fam$V3!=0),c("V3","V1")]
sub4 = fam[which(fam$V4!=0),c("V4","V1")]
colnames(sub2) = c("IID","FID")
colnames(sub3) = c("IID","FID")
colnames(sub4) = c("IID","FID")
sub = rbind(sub2,sub3,sub4)
fun = function(FID) paste(sort(unique(FID)),collapse=",")
adon = group_by(sub,IID) %>% summarise(FID=fun(FID)) %>% ungroup()
adon = as.data.frame(adon)



### convert fam file into vertices

kid_dad = fam[which(fam$V3!=0),c("V2","V3")]
colnames(kid_dad) = c("IID1","IID2")
kid_dad$type = "kid_dad"

kid_mom = fam[which(fam$V4!=0),c("V2","V4")]
colnames(kid_mom) = c("IID1","IID2")
kid_mom$type = "kid_mom"

decl = rbind(kid_dad,kid_mom)

# sort values for comparability  - not sure whether it is useful yet
mm = t(apply(decl[,c("IID1","IID2")],1,function(x) sort(x)))
decl[,c("IID1","IID2")] = mm; rm(mm)
#head(decl)


### deal with the IBD file

gene = full[,c("IID1","IID2")]
mm = t(apply(gene,1,function(x) sort(x)))
gene = as.data.frame(mm,stringsAsFactors = F)
colnames(gene) = c("IID1","IID2")

# set the PI_HAT and Z1 thresholds

### Z1 threshold for PO (parent-offspring)
z1_thr = 0.8 #0.75  # Z1 values above are used for parent-offspring relationship to be called
### determine thresholds for FS (full-sibs)
fs_thr = c(0.35,0.65, 0.35,0.65) # for full siblings: PI_HAT low, PI_HAT high, Z1 low, Z1 high
### determine thresholds for TW (twins)
tw_thr = 0.8  # threshold used for twin relationship detection (PI_HAT)

#par(mfrow=c(1,1))
#plot(full$Z1,full$PI_HAT,pch=19,cex=0.5,col=rgb(0.3, 0.3, 0.3,0.3)); grid()
#polygon(x=c(fs_thr[3:4],fs_thr[4:3],fs_thr[3]),y = c(fs_thr[1],fs_thr[1],fs_thr[2],fs_thr[2],fs_thr[1]),col = NA) #FS area
#polygon(x=c(z1_thr,1.01,1.01,z1_thr,z1_thr),y = c(fs_thr[1],fs_thr[1],fs_thr[2],fs_thr[2],fs_thr[1]),col = NA) # PO
#polygon(x=c(-0.01,0.2,0.2,-0.01,-0.01),y = c(tw_thr,tw_thr,1.01,1.01,tw_thr),col = NA) # TW


gene$type = NA
gene$type[which(full$Z1>z1_thr)] = "PO"
gene$type[which(full$PI_HAT>tw_thr)] = "TW"
gene$type[which((full$Z1>=fs_thr[3])&(full$Z1<=fs_thr[4])&(full$PI_HAT>=fs_thr[1])&(full$PI_HAT<=fs_thr[2]))] = "FS"
gene$type[which( (full$PI_HAT>=fs_thr[1])&(is.na(gene$type)))] = "AMBG"
# reduce
gene = gene[which(full$PI_HAT>=fs_thr[1]),] # (PI_HAT low) 
#table(gene$type,useNA = "a")



### combo dataset

comb = rbind(decl,gene)#; head(comb)
comb$unqid = paste(comb$IID1,comb$IID2,sep="-")

fun = function(type) paste(type,collapse=",")
dtf = group_by(comb,unqid) %>% summarise(n=n(),new_type=fun(type)) %>% ungroup()
dtf = as.data.frame(dtf)


tmp1 = lapply(dtf[,"unqid"],function(x) unlist(strsplit(x,"-")))
tmp11 = do.call(rbind.data.frame, tmp1)
colnames(tmp11) = c("IID1","IID2")
origin = cbind(tmp11,type = dtf$new_type)
origin$IID1 = as.character(origin$IID1)
origin$IID2 = as.character(origin$IID2)
origin$type = as.character(origin$type)
#table(origin$type)

#####
##### THE GREAT DIGESTOR
##### [INDISCRIMINATELY extracts ALL family clusters determined by both genetics and pedigree]
##### [not to be confused with digestor script below, which only extracts suspicious clusters]
#####

lst_GREAT = list()      # the great (aka full) list of clusters
clust_szs_GREAT = NULL  # meta data for the great list of clusters
dat = origin            # "dat" will be gradually trimmed until dissapears
counter = 0             # counts number of clusters, assigns IDs
while(nrow(dat)>0) {
        counter = counter +1
       # print(counter)
        rixs = 1 # initiation index
        go = TRUE  # allows inner cycle to continue
        while(go) {
                l = length(rixs)
                iids = unique(as.character(as.matrix(dat[rixs,c("IID1","IID2")]))) #; length(iids)
                rixs = unique(which( (dat$IID1 %in% iids)|(dat$IID2 %in% iids) )) #; length(rixs)
                go = length(rixs)>l
        }
        clust_szs_GREAT = c(clust_szs_GREAT,length(rixs))
        lst_GREAT[[counter]] = dat[rixs,]  # save
        dat = dat[-rixs,] # remove from original
        rm(iids,rixs,l,go) # cleanup
}

#length(lst_GREAT)
#table(clust_szs_GREAT)


#####
##### THE SMART DIGESTOR
##### [ONLY extracts those family clusters, that have signs of problems]
##### [not to be confused with digestor script abowe, which extracts all clusters]
#####

# relationship types that are expected (i.e., both pedigree and genetics must agree):
ok_types = c("kid_dad,PO","kid_mom,PO","FS")
#  relationship types that are suspicious/problematic/imposible
susp_types = unique(origin$type[which(!origin$type %in% ok_types)]) # = "kid_dad","kid_mom","kid_mom,TW","PO","TW")  
if(length(susp_types)==0) susp_types = c("kid_dad","kid_mom","kid_mom,TW","PO","TW")

lst_SMART = list()      # the great (aka full) list of clusters
clust_szs_SMART = NULL  # meta data for the great list of clusters
dat = origin            # "dat" will be gradually trimmed until dissapears
counter = 0             # counts number of clusters, assigns IDs
while( sum(dat$type %in% susp_types)>0 ) {  # different than in script above****
        counter = counter +1
        #print(counter)
        rixs = which(dat$type %in% susp_types)[1] # initiation index
        go = TRUE  # allows inner cycle to continue
        while(go) {
                l = length(rixs)
                iids = unique(as.character(as.matrix(dat[rixs,c("IID1","IID2")]))) #; length(iids)
                rixs = unique(which( (dat$IID1 %in% iids)|(dat$IID2 %in% iids) )) #; length(rixs)
                go = length(rixs)>l
        }
        clust_szs_SMART = c(clust_szs_SMART,length(rixs))
        lst_SMART[[counter]] = dat[rixs,]  # save
        dat = dat[-rixs,] # remove from original
        rm(iids,rixs,l,go) # cleanup
}

#length(lst_SMART)
#table(clust_szs_SMART)

#####
##### EXPORT/SAVE for future use
#####

# export FLAG FILE
write.table(x = flg,file = flg_fam_file,row.names = F,col.names = T,sep="\t",quote=F)

# export UPDATED FAM FILE
write.table(x = fam,file = upd_fam_file,row.names = F,col.names = F,sep="\t",quote=F)

# export working objects
save(list = c("full","fam_orig","z1_thr","fs_thr","tw_thr","sex","sex_upd","fam","flg","upd","del",
              "y_thr","f_thr","problm_chips","lst_GREAT","lst_SMART","susp_types",
              "clust_szs_GREAT","clust_szs_SMART","adon","mom_ids","dad_ids","Yfemales","Ymales"),
     file = rep_fileRdat)

# now generate report using : "report_on_clusterFuckups.Rmd"


