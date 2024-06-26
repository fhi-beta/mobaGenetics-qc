# Script to generate a rsID SNP conversion table
library(dplyr)
library(tidyr)

# load illumina provided rsid conversion list
illumina.rsid <- read.table('/home/oyvind/harvest/media/local-disk2/helgeland/rotterdam1/aux/GSA-24-v1.0/loci-name-to-rsid-conversion/GSA-24v1-0_A1_b144_rsids.txt', 
                            header = T, stringsAsFactors = F, sep='\t')

# retaining only the first rsID where there are multiple options
illumina.single.rsid <- illumina.rsid %>% separate(col = RsID, into = "RsID_illumina", sep=",")

# convert to use the original ID if dot notation in Illumina file
illumina.single.rsid$RsID_out_tmp <- ifelse(illumina.single.rsid$RsID=='.', illumina.single.rsid$Name, illumina.single.rsid$RsID)

# make IDs of triallelic markers unique by suffixes (to avoid PLINK problems)
illumina.single.rsid$RsID_out <- make.unique(illumina.single.rsid$RsID_out_tmp)

# generate final outfile
out.file <- illumina.single.rsid %>% select(Name, RsID_out)

# write file
write.table(out.file, file = '/home/oyvind/harvest/media/local-disk2/helgeland/rotterdam1/recode-files/recode-rsid.txt', 
            col.names = F, row.names = F, quote = F, sep=' ')
