plot_depth_by_amplicon <- function(df)
{
    ggplot(df, aes(x=sample, y=mean_coverage)) + 
        geom_point() +
        scale_y_log10() + 
        facet_wrap(. ~ amplicon_id) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

        ggsave("qc_sequencing/depth_by_amplicon_sample.pdf", width=15, height=10)
}

plot_depth_by_amplicon_and_ct <- function(df)
{
    ggplot(df, aes(x=ct, y=mean_coverage)) + 
        geom_point() +
        scale_y_log10() + 
        facet_wrap(. ~ amplicon_id) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

        ggsave("qc_sequencing/depth_by_amplicon_ct.pdf", width=15, height=10)
}

plot_fraction_covered_by_amplicon <- function(df)
{
    print(head(df))
    ggplot(df, aes(x=sample, y=fraction_covered)) + 
        geom_point() +
        facet_wrap(. ~ amplicon_id) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

        ggsave("qc_sequencing/covered_by_amplicon.pdf", width=15, height=10)
}

plot_depth_per_base <- function(df)
{
    # bedtools coverage reports the position along the interval, adding start gives
    # us the position along the genome
    ggplot(df, aes(x=start + position, y=depth)) + 
        geom_line() +
        scale_y_log10() + 
        facet_wrap(. ~ sample, ncol=1) + 
        xlab("Genome position (nt)") +
        theme_bw()
        
        ggsave("qc_sequencing/depth_by_position.pdf", width=8, height=16)
}

plot_alt_frequency <- function(df)
{
    # bedtools coverage reports the position along the interval, adding start gives
    # us the position along the genome
    ggplot(df, aes(x=position, y=alt_frequency, color=depth)) + 
        geom_line() +
        facet_wrap(. ~ sample, ncol=1) + 
        #scale_colour_gradient(low = "white", high = "black") +
        xlab("Genome position (nt)") +
        ylab("Alt Frequency") +
        theme_bw()
        
        ggsave("qc_sequencing/alt_frequency.pdf", width=8, height=24)
}

plot_genome_completeness <- function(df)
{
    ggplot(df, aes(x=sample_name, y=pct_covered_bases)) +
        geom_col() +
        xlab("Sample") +
        ylab("Genome Completeness (%)") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        ggsave("qc_analysis/genome_completeness.pdf", width=12, height=8)
}

# load all the files matching this name into one data frame
# from: https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once
read_glob <- function(dir, pattern) {
    files = list.files(dir, pattern=pattern)
    print(files)
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

library(ggplot2)
library(dplyr)

if(!interactive()) {
    #mean_coverage <- read_glob("qc_sequencing", "*.mean_coverage.bed")
    #plot_depth_by_amplicon(mean_coverage)

    args = commandArgs(trailingOnly=TRUE)
    if(length(args) == 0) {
        q()
    }

    if(args[1] == "depth_by_position") {
        df <- read_glob("qc_sequencing", "*.per_base_coverage.bed")
        plot_depth_per_base(df)
    } else if(args[1] == "amplicon_covered_fraction") {
        df <- read_glob("qc_sequencing", "*.amplicon_coverage.bed")
        plot_fraction_covered_by_amplicon(df)
    } else if(args[1] == "alt_allele_frequency") {
        df <- read_glob("qc_sequencing", "*.fpileup.tsv")
        plot_alt_frequency(df)
    } else if(args[1] == "genome_completeness") {
        df <- read.csv("qc_analysis/merged.qc.csv")
        plot_genome_completeness(df)
    }

}
