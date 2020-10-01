
plot_depth_by_amplicon_and_ct <- function(df, metadata, outname)
{
    # merge df with metadata
    merged = dplyr::inner_join(df, metadata, by = "sample")

    ggplot(merged, aes(x=ct, y=mean_depth)) +
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
    merged = dplyr::left_join(df, metadata, by = "sample")

    # construct new column containing the facet label
    merged$label = paste(merged$sample, " Ct: ", merged$ct, sep="")

    num_samples = length(unique(merged$sample))
    plots_per_page = min(8, num_samples)

    pdf(outname, width=8, height=2*plots_per_page)

    for(i in seq(1, ceiling(num_samples / plots_per_page))) {
        # bedtools coverage reports the position along the interval, adding start gives
        # us the position along the genome
        p <- ggplot(merged, aes(x=start + position, y=depth)) +
            geom_line() +
            scale_y_log10() +
            facet_wrap_paginate(. ~ label, nrow=plots_per_page, ncol=1, page=i) +
            xlab("Genome position (nt)") +
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
library(argparse)

if(!interactive()) {

# Main
    parser <- ArgumentParser(description="Generate various QC plots")

    parser$add_argument('--directory', '-d', default="qc_sequencing",
                        help="directory containing the input files")
    parser$add_argument('--plot-type', '-t', default="depth_by_position",
                        help="the type of plot to generate")
    parser$add_argument("--output", "-o", default="",
                        help="name of the output plot")
    parser$add_argument("--metadata", "-m", default="",
                        help="optional metadata file to annotate the plot with")
    args <- parser$parse_args()
    
    if(args$output == "") {
        q()
    }

    # metadata is optional so we make a default table here
    metadata = data.frame(sample=character(),
                          ct=double())

    if(args$metadata != "") {
        metadata.raw <- read.table(args$metadata, header=T, sep="\t")

        # clean metadata of unknown Cts
        metadata = subset(metadata.raw, ct != "unknown")

        # convert Cts to numeric
        metadata$ct = as.numeric(as.character(metadata$ct))
    }

    if(args$plot_type == "depth_by_position") {

        df <- read_glob("qc_sequencing", "*.per_base_coverage.bed")
        plot_depth_per_base(df, metadata, args$output)

    } else if(args$plot_type == "negative_control_depth_by_position") {

        df <- read_glob("qc_sequencing_negative_control", "*.per_base_coverage.bed")
        plot_depth_per_base(df, metadata, args$output)

    } else if(args$plot_type == "amplicon_depth_by_ct") {

        df <- read_glob("qc_sequencing", "*.amplicon_depth.bed")
        plot_depth_by_amplicon_and_ct(df, metadata, args$output)

    } else if(args$plot_type == "amplicon_covered_fraction") {

        df <- read_glob("qc_sequencing", "*.amplicon_coverage.bed")
        plot_fraction_covered_by_amplicon(df, args$output)

    } else if(args$plot_type == "genome_completeness_by_ct") {

        df <- read_qc_csv("qc_analysis/merged.qc.csv")
        plot_genome_completeness_by_ct(df, metadata, args$output)
    }
}
