##
#
# This script convers the chromosome naming of the strand files to the plink format. This needs to be run only once.
#
##

folder <- "qc-pipeline/resources/ab_conversion"

for (filename in list.files(folder)) {
  
  if (endsWith(filename, ".strand.gz")) {
    
    print(paste0("processing ", filename))
  
    strand_table <- read.table(
      file = file.path(folder, filename),
      header = F
    )
    
    strand_table$V2[strand_table$V2 == "X"] <- 23
    strand_table$V2[strand_table$V2 == "Y"] <- 24
    strand_table$V2[strand_table$V2 == "XY"] <- 25
    strand_table$V2[strand_table$V2 == "MT"] <- 26
    
    new_name <- paste0(substr(filename, 1, nchar(filename) - 10), "_chr.strand.gz")
    
    write.table(
      x = strand_table,
      file = gzfile(file.path(folder, new_name)),
      sep = "\t",
      col.names = F,
      row.names = F,
      quote = F
    )
    
  }
}

