debug <- T

if (debug){
    args <- c("/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_make_phasing_tables.best", "/mnt/work/oystein/tmp/phasing_test/report.md")
} else {
    args <- commandArgs(TRUE)
}

chrs <- sapply(1:22, function(i) paste0(args[1],"chr", i))

md_file <- args[1]
docs_folder <- dirname(md_file)
plots_folder <- paste0(docs_folder, "/plots/")


if(!dir.exists(docs_folder)){
    dir.create(docs_folder)
}


if(!dir.exists(pipeout_folder)){
    dir.create(pipeout_folder)
}

if(!dir.exists(plots_folder)){
    dir.create(plots_folder)
}


c1 <- read.table(chrs[1], header = T)
c1_0<- subset(c1, parents_in_batch == 0)
c1_1<- subset(c1, parents_in_batch == 1)
c1_2<- subset(c1, parents_in_batch == 2)