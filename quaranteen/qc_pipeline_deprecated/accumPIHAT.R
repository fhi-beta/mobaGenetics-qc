args = commandArgs(trailingOnly=TRUE)


### Select individuals who are too-related to many others (accumulated PIHAT)
### Jonas 20170115
### Runtime: around 13 minutes

#### Description:
# 1) imports a FULL IBD pairwise comparisson file for ALL individuals (.geno)
#    also imports a .fam file to generate a FID,IID in the output
# 2) eliminates too-related pairs (edges), but keeps other edges of the same individuals
# 3) splits the input file to improve the speed (was optimised on full MoBa12)
# 4) estimates the accumulated PI_HAT for each individual
# 5) normalizes accumPIHAT for sample/cohort size
# 6) generates visual report
# 7) outputs all accumPIHAT values for reporting
# 8) outputs individuals-to-be-excluded based on hard-set threshold

### MANUALLY SET THE PARAMETERS AND ARGUMENTS
############################################################
# ###### input files
# ibd_file = ".geno"
# fam_file = "~/Biostuff/temp_export/deleteme_moba12mrg.fam"
# ###### output files
# out_file = ".txt"  # "remove file" to be used in PLINK
# report_dir = "reporting/", where pdf, txt and png reports will be stored
# ###### parameters
# pih_thr = 0.2   # PI_HAT threshold to ignore too-closely related pairs
# hard_thr = 3      # hard-set (non-SD) threshold for normalized accumPIHAT to determine outliers
############################################################


ibd_file        = args[1]
fam_file        = args[2]
out_file        = args[3]
report_dir      = args[4]  # pass no slash 
pih_thr         = as.numeric(args[5])
hard_thr        = as.numeric(args[6])
res_files       = args[7]   #  .png, .pdf. and .txt will be added to this trunk
		# Note that out_file and report_dir/res_files.txt should not be the same file ...

library(dplyr)
library(ggplot2)

# load fam file (to get famIDs)
fam = read.table(fam_file,h=F,stringsAsFactors = F)

plot_report = paste(report_dir, paste(res_files,".png",sep=""), sep="/")
pdf_report = paste(report_dir, paste(res_files,".pdf",sep=""), sep="/")
txt_report = paste(report_dir, paste(res_files,".txt2",sep=""), sep="/")

# load the genetic-relatedness matrix/dataframe
a = read.table(ibd_file,h=T,stringsAsFactors = F)
#load("~/Biostuff/temp_export/rez.RData")
#nrow(a); head(a)

### exclude all PAIRS that are too related
# (but leave all other pairs with those individuals)
a = a[which(a$PI_HAT<=pih_thr),] # remove real-family edges
#nrow(a); head(a)

# after loading it should take 11 min in total to run everything else

# extract unique sample IDs
unq_ids = unique(c(unique(a$IID1),unique(a$IID2)))
results = data.frame(IID=unq_ids,accum_PI_HAT = 0,stringsAsFactors = F) # here results will be stored
rm(unq_ids)
head(results)

# split the file into subfiles for faster search
n_grps = 12 # optimised for moba12 (should take 8 min)
grps = sample(1:n_grps,nrow(a),replace = T) # randomly assign to groups


for (gr in 1:n_grps) {
        print(gr)
        # isolate a subfile
        sub = a[which(grps==gr),]

        # restructure the file
        df1 = sub[,c("IID1","PI_HAT")]
        df2 = sub[,c("IID2","PI_HAT")]
        colnames(df1) = c("IID","PI_HAT")
        colnames(df2) = c("IID","PI_HAT")
        sub = rbind(df1,df2)

        # accumulated PI_HAT for each sample
        df = group_by(sub,IID) %>% summarise(s=sum(PI_HAT)) %>% ungroup() # n=n(),
        df = as.data.frame(df)

        # update the final results
        tmp = merge(results,df,by="IID",all=T)
        results = tmp[,1:2]
        results[,2] = apply(tmp[,2:3],1,sum)

        # cleanup
        rm(df,df1,df2,tmp,sub)
}

# spped optimisation (with n_grps = 5,10,25,50,100,200):
# 143 * 5 / 60    # 12 minutes to run
# 45 * 10 / 60    # 7 minutes to run
# 19 * 25 / 60    # 8 minutes to run
# 8.6 * 50 / 60   # 7 minutes to run
# 5.2 * 100 / 60  # 9 minutes to run
# 4.46 * 200 / 60 # 15 minutes to run


results$normalised = results$accum_PI_HAT / nrow(results) # max number of pairs for each ind
thr_accpihat = hard_thr # mean(results$normalised) + sd_thr * sd(results$normalised) # a bit arbitrary

# preview non-normalized
#hist(results$accum_PI_HAT,breaks=100,col="grey")
#hist(results$normalised,breaks=100,col="grey")

# visual reporting

# dirty hack: results$normalised with only NA will crash hist(). Only execture when there is at least one value
if (sum(!is.na(results$normalised)) > 0) {
   pdf(pdf_report)
   h = hist(results$normalised,breaks=100,col="grey",
       main="accumulated PI_HAT normalized to sample size",
       xlab="accumPI_HAT / N",ylab="individuals")
   abline(v=thr_accpihat,col="orange1",lwd=2,lty=2)
   x_max = max(results$normalised)
   n_out = sum(results$normalised > thr_accpihat)
   n_tot = length(results$normalised)
   text(x=x_max*0.6,y=max(h$counts)*0.9,paste("(PI_HAT>",pih_thr," edges were excluded)",sep=""))
   text(x=x_max*0.6,y=max(h$counts)*0.8,paste("threshold =",hard_thr,sep=" ")) # ,"SD"
   text(x=x_max*0.6,y=max(h$counts)*0.7,paste("N excluded = ",n_out,sep=""))
   text(x=x_max*0.6,y=max(h$counts)*0.6,paste("N total = ",n_tot,sep=""))
   dev.off()

   p = ggplot(results, aes(x=normalised)) + geom_histogram() +
      stat_bin(aes(label=ifelse(..count..==0, "", ..count..), y=..count..),
             geom="text", vjust=-0.5, size=3) +
   	     geom_vline(aes(xintercept=thr_accpihat), colour="red") +
	     theme_bw() + ggtitle(expression(paste("Accumulated ", hat(pi), ", normalized"))) + xlab(NULL)
      ggsave(plot_report, width=16, height=10, units="cm", plot=p)
} else {"Warning: Skipped ploting as results normalised had no values"} # end of dirty code to skip plotting

# output
bad_samples = results[which(results$normalised > thr_accpihat),"IID"]
out_bad_samples = fam[which(fam$V2 %in% bad_samples),1:2]
out_all_report = results

# write out the output
write.table(out_bad_samples,file = out_file,row.names = F,col.names = F,quote=F,sep="\t")
write.table(out_all_report,file = txt_report,row.names = F,col.names = T,quote=F,sep="\t")

# note 1:  be brave while excluding individuls, this is temporary exclusion!
