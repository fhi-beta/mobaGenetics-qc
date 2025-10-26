debug <- F

if (debug){
    args <- c(
        "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_make_phasing_tables.best", 
        "/mnt/work/oystein/tmp/phasing_test/report.md",
        "Phasing report debug")
} else {
    args <- commandArgs(TRUE)
}

input_files <- sapply(1:22, function(i) paste0(args[1],".chr", i))

md_file <- args[2]
title <- args[3]
docs_folder <- dirname(md_file)
plots_folder <- paste0(docs_folder, "/plot/")


if(!dir.exists(docs_folder)){
    dir.create(docs_folder)
}

if(!dir.exists(plots_folder)){
    dir.create(plots_folder)
}


write(
  x = paste("#", title),
  file = md_file,
  append = F
)



rates <- list(
  c("r_phasing_hom", "Phasing error rates"),
  c("r_mendel", "Mendelian error rates")
)

plot_density <- function(numbers, title, xlab, absolute_filepath, relative_filepath) {
    numbers <- na.omit(numbers)
    png(absolute_filepath)
    breaks_seq <- seq(0, 1, length.out = 100)
    hist(numbers,
    main = title,
    xlab = xlab,
    ylab = "Frequency",
    xlim = c(0, 1),
    col = "lightblue",
    border = "black",
    breaks = breaks_seq)

    dev.off()
    write(
        x = paste0("![](",relative_filepath,")"),
        file = md_file,
        append = T
    )
}

chromosome_tables <- lapply(input_files, function(file) read.table(file, header = TRUE))

write_rates_table_row <- function(p, rates){
    rates <- na.omit(rates)
    write(
        x = paste("|", p,"|", signif(min(rates), digits=2), "|", signif(max(rates), digits=2), "|", signif(median(rates), digits=2), "|", signif(mean(rates), digits=2), "|"),
        file = md_file,
        append = T
    )
}

for (rate in rates){
    write(
        x = paste0("## ", rate[2], " per chromosome and number of parents in same batch as child"),
        file = md_file,
        append = T
    )
    for (chr in 1:22) {
        tab <- chromosome_tables[[chr]]
        write(
            x = paste0("### ", rate[2], ", chromosome ", chr),
            file = md_file,
            append = T
        )
        write(
            x = "| Parents in batch | min | max | median | mean |\n| --- | --- | --- | --- | --- |",
            file = md_file,
            append = T
        )
        for (p in 0:2){
            tab_p <- subset(tab, parents_in_batch == p)
            write_rates_table_row(p, tab_p[[rate[1]]])
        }
        write(
            x = "\n",
            file = md_file,
            append = T
        )
        for (p in 0:2){
            plur <- "s"
            if (p==1){
                plur <- ""
            }
            tab_p <- subset(tab, parents_in_batch == p)
            filename <- paste0(gsub(" ", "_", tolower(rate[2])),"_chr",chr,"_",p,"parents.png")
            title <- paste0("Chromosome ", chr,", ",p, " parent", plur, " in batch")
            aboslute_path <- paste0(plots_folder, filename)
            relative_path <- paste0("plots/", filename)
            plot_density(tab_p[[rate[1]]], title, rate[2], aboslute_path, relative_path)

        }
    }
}
