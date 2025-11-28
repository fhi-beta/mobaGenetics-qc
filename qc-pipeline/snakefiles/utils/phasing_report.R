library(dplyr)
debug <- F
if (debug){
    args <- c(
        "/mnt/archive3/phasing_test/phase_related_orig_imputation/mod7_phase_check.chr21",
        "/mnt/archive3/phasing_test/phase_related_orig_imputation/expected_all_relations",
        "/mnt/archive3/phasing_test/debug/phase_report/phase_report.chr21.md", 
        "Phasing report, pre-phase related samples",
        "phasing_hom",
        "Phasing error rates",
        "21")
} else {
    args <- commandArgs(TRUE)
}


chr_input <- args[7]

if(chr_input != "0"){
    single_file <- TRUE
    chr <- chr_input
} else{
    single_file <- FALSE
}


if(single_file){
    input_file <- args[1]
} else{
    input_files <- sapply(1:22, function(i) paste0(args[1],".chr", i))
}


relations_file <- args[2]
md_file <- args[3]

title <- args[4]
rate <- c(args[5], args[6])
docs_folder <- dirname(md_file)
plots_folder <- paste0(docs_folder, "/phasing_plots/")



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

if (single_file){
    chromosome_table <- read.table(input_file, header = TRUE)
    rel <- read.table(relations_file, header = TRUE)
    chromosome_table <- chromosome_table %>% left_join(rel %>% select(iid, shared_chips), by = "iid")
    # if (filter_psam != "0")
    # {
    #     psam <- read.table(filter_psam)
    #     ptrios <- subset(psam, V3 != "0" & V4 != "0") %>% select(iid = V2, pat = V3, mat = V4)
    #     chromosome_table <- chromosome_table %>% semi_join(ptrios, by = c("iid", "pat", "mat"))

    # }
} else {chromosome_tables <- lapply(input_files, function(file) read.table(file, header = TRUE))
rel <- read.table(relations_file, header = TRUE)
chromosome_tables <- lapply(chromosome_tables, function(tab) tab %>% left_join(rel %>% select(iid, shared_chips), by = "iid"))
}


write_rates_table_row <- function(p, rates){
    rates <- na.omit(rates)
    write(
        x = paste("|", p,"|", signif(min(rates), digits=2), "|", signif(max(rates), digits=2), "|", signif(median(rates), digits=2), "|", signif(mean(rates), digits=2), "|"),
        file = md_file,
        append = T
    )
}

write_combined_rates_table_row <- function(p, trios, sites, errors){
    write(
        x = paste("|", p,"|", trios, "|", sites, "|", errors, "|", signif(errors/sites, digits=2), "|"),
        file = md_file,
        append = T
    )
}


# write(
#     x = paste0("## ", rate[2], " per chromosome and number of parents in same batch as child"),
#     file = md_file,
#     append = T
# )

write_chromosome_table <- function(tab, rate, chr){
    write(
            x = paste0("### ", rate[2], ", chromosome ", chr),
            file = md_file,
            append = T
        )
        write(
            x = paste0("#### Combined error rates:\n - Trios = Number of trios\n - Sites = Total number of sites where error can occur (summed over all trios)\n - Errors = Total number of errors\n"),
            file = md_file,
            append = T
        )
        
        write(
            x = "| Shared chips | Trios | Sites | Errors | Error rate |\n|:----|:----|:----|:----|:----|",
            file = md_file,
            append = T
        )
        for (p in 0:2){
            tab_p <- subset(tab, shared_chips == p)
            trios <- nrow(tab_p)
            sites <- sum(tab_p[[paste0("n_",rate[1])]])
            errors <- sum(tab_p[[paste0("e_",rate[1])]])
            write_combined_rates_table_row(p, trios, sites, errors)
        }
        tab_p <- tab
        trios <- nrow(tab_p)
        sites <- sum(tab_p[[paste0("n_",rate[1])]])
        errors <- sum(tab_p[[paste0("e_",rate[1])]])
        write_combined_rates_table_row("Total", trios, sites, errors)
        write(
            x = "\n",
            file = md_file,
            append = T
        )
        write(
            x = paste0("#### Individual error rates:"),
            file = md_file,
            append = T
        )
        write(
            x = "| Shared chips | min | max | median | mean |\n|:----|:----|:----|:----|:----|",
            file = md_file,
            append = T
        )
        for (p in 0:2){
            tab_p <- subset(tab, shared_chips == p)
            write_rates_table_row(p, tab_p[[paste0("r_",rate[1])]])
        }
        write_rates_table_row("Total", tab[[paste0("r_",rate[1])]])
        write(
            x = "\n",
            file = md_file,
            append = T
        )
        
        write(
            x = paste0("#### Histogram of individual error rates:"),
            file = md_file,
            append = T
        )
        for (p in 0:2){
            plur <- "s"
            if (p==1){
                plur <- ""
            }
            tab_p <- subset(tab, shared_chips == p)
            filename <- paste0(gsub(" ", "_", tolower(rate[2])),"_chr",chr,"_",p,"shared_chips.png")
            title <- paste0("Chromosome ", chr,", ",p, " shared chip", plur)
            aboslute_path <- paste0(plots_folder, filename)
            relative_path <- paste0("phasing_plots/", filename)
            plot_density(tab_p[[paste0("r_",rate[1])]], title, rate[2], aboslute_path, relative_path)

        }
        filename <- paste0(gsub(" ", "_", tolower(rate[2])),"_chr",chr,"_total.png")
        title <- paste0("Chromosome ", chr,", all samples")
        aboslute_path <- paste0(plots_folder, filename)
        relative_path <- paste0("phasing_plots/", filename)
        plot_density(tab[[paste0("r_",rate[1])]], title, rate[2], aboslute_path, relative_path)
    }

if(single_file){
    write_chromosome_table(chromosome_table, rate, chr)
} else{
for (chr in 1:22) {
    tab <- chromosome_tables[[chr]]
    write_chromosome_table(tab, rate, chr)
}
}


