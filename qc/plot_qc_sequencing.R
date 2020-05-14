plot_depth_by_amplicon <- function(df, outname)
{
    ggplot(df, aes(x=sample, y=mean_coverage)) + 
        geom_point() +
        scale_y_log10() + 
        facet_wrap(. ~ amplicon_id) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

        ggsave(outname, width=15, height=10)
}

plot_depth_by_amplicon_and_ct <- function(df, metadata, outname)
{
    # merge df with metadata   
    merged = dplyr::inner_join(df, metadata, by = "sample")

    summary <- merged %>% group_by(sample, amplicon_id, ct)  %>% summarise(mean_coverage=mean(depth))

    ggplot(summary, aes(x=ct, y=mean_coverage)) + 
        geom_point() +
        scale_y_log10() + 
        facet_wrap(. ~ amplicon_id) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

        ggsave(outname, width=15, height=10)
}

plot_fraction_covered_by_amplicon <- function(df, outname)
{
    ggplot(df, aes(x=sample, y=fraction_covered)) + 
        geom_point() +
        facet_wrap(. ~ amplicon_id) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

        ggsave(outname, width=15, height=10)
}

plot_depth_per_base <- function(df, metadata, outname)
{
    # merge df with metadata   
    merged = dplyr::full_join(df, metadata, by = "sample")

    # construct new column containing the facet label
    merged$label = paste(merged$sample, " Ct: ", merged$ct, sep="")

    num_samples = length(unique(merged$sample))
    plots_per_page = 8

    pdf(outname, width=8, height=16)

    for(i in seq(1, ceiling(num_samples / plots_per_page))) {
        # bedtools coverage reports the position along the interval, adding start gives
        # us the position along the genome
        p <- ggplot(merged, aes(x=start + position, y=depth)) + 
            geom_line() +
            scale_y_log10() + 
            facet_wrap_paginate(. ~ label, nrow=8, ncol=1, page=i) + 
            xlab("Genome position (nt)") +
            theme_bw()
        print(p)
    }
    dev.off()
}

plot_alt_frequency <- function(df, outname)
{
    pdf(outname, width=8, height=24)
    num_samples = length(unique(df$sample))
    plots_per_page = 8

    for(i in seq(1, ceiling(num_samples / plots_per_page))) {
        p <- ggplot(df, aes(x=position, y=alt_frequency, color=depth)) + 
            geom_line() +
            facet_wrap_paginate(. ~ sample, nrow=plots_per_page, ncol=1, page=i) + 
            #scale_colour_gradient(low = "white", high = "black") +
            xlab("Genome position (nt)") +
            ylab("Alt Frequency") +
            theme_bw()
        print(p)
    }
    dev.off()
}

plot_genome_completeness_by_ct <- function(df, metadata, outname)
{
    # merge df with metadata   
    merged = dplyr::inner_join(df, metadata, by = "sample")

    ggplot(merged, aes(x=ct, y=pct_covered_bases)) +
        geom_point() +
        xlab("Ct") +
        ylab("Genome Completeness (%)") +
        ylim(0, 100) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        ggsave(outname, width=12, height=8)
}

# load all the files matching this name into one data frame
# from: https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once
read_glob <- function(dir, pattern) {
    files = list.files(dir, pattern=pattern)
    read_sample_file <- function(x) {
        fn <- paste(dir, x, sep="/")
        sample_name <- strsplit(x, ".", fixed=T)[[1]][1]
        df <- read.table(fn, stringsAsFactors = F, header=T)
        df$sample <- sample_name
        return(df)
    }

    df = do.call(rbind, lapply(files, read_sample_file))
    return(df)
}

# read artic-nextflow's QC file
read_qc_csv <- function(fn) {
    df <- read.csv(fn)
    # standardize sample field
    df$sample <- df$sample_name
    return(df)
}

library(ggplot2)
library(dplyr)
library(ggforce)

if(!interactive()) {
    #mean_coverage <- read_glob("qc_sequencing", "*.mean_coverage.bed")
    #plot_depth_by_amplicon(mean_coverage)

    args = commandArgs(trailingOnly=TRUE)
    if(length(args) == 0) {
        q()
    }

    command = args[1]
    prefix = args[2]

    # metadata is optional so we make a default table here
    metadata = data.frame(sample=character(),
                          ct=double())

    if(length(args) == 3) {
        metadata.raw <- read.table(args[3], header=T)
        
        # clean metadata of unknown Cts
        metadata = subset(metadata.raw, ct != "unknown")
    
        # convert Cts to numeric
        metadata$ct = as.numeric(as.character(metadata$ct))
    }

    if(command == "depth_by_position") {
        
        df <- read_glob("qc_sequencing", "*.per_base_coverage.bed")
        outname = sprintf("plots/%s_depth_by_position.pdf", prefix)
        plot_depth_per_base(df, metadata, outname)

    } else if(command == "amplicon_depth_by_ct") {
        
        df <- read_glob("qc_sequencing", "*.per_base_coverage.bed")
        outname = sprintf("plots/%s_amplicon_depth_by_ct.pdf", prefix)
        plot_depth_by_amplicon_and_ct(df, metadata, outname)

    } else if(command == "amplicon_covered_fraction") {

        df <- read_glob("qc_sequencing", "*.amplicon_coverage.bed")
        outname = sprintf("plots/%s_amplicon_covered_fraction.pdf", prefix)
        plot_fraction_covered_by_amplicon(df, outname)

    } else if(command == "alt_allele_frequency") {
 
        df <- read_glob("qc_sequencing", "*.fpileup.tsv")
        outname = sprintf("plots/%s_alt_frequency.pdf", prefix)
        plot_alt_frequency(df, outname)

    } else if(command == "genome_completeness_by_ct") {
        
        df <- read_qc_csv("qc_analysis/merged.qc.csv")
        outname = sprintf("plots/%s_genome_completeness_by_ct.pdf", prefix)
        plot_genome_completeness_by_ct(df, metadata, outname)

    }
}
