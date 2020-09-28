library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(argparse)

default.heatmap.theme <- function(base_size = 24, base_family = 'Helvetica', xaxis.angle = 90) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text.x = element_text(angle = xaxis.angle, size = rel(0.8), hjust = 1)
    )
}


get_variant_annotation_files <- function(path='.', pattern='multianno.txt') {
  if (is.null(path)) stop("Mandatory argument path is missing")

  files <- list.files(
    path=path,
    pattern=pattern,
    full.names=TRUE,
    recursive=TRUE
  )

  return(files)
}


get_variant_annotation_data <- function(files) {
  if (is.null(files)) stop("Mandatory argument files is missing")

  data <- data.frame()
  for (i in 1:length(files)) {
    if (file.info(files[i])$size > 0) {
      tmp_data <- read.table(
        file=files[i],
        header=TRUE,
        as.is=TRUE,
        sep='\t',
        quote='\"')
    } else {
      next
    }
    data = rbind(data, tmp_data)
  }
  colnames(data) <- c('chr', 'start', 'end', 'ref', 'alt', 'func', 'gene', 'gene_detail', 'exonic_func', 'aa_change', 'sample')
  return(data)
}


get_aa_change <- function(variant=NULL) {
  if (is.null(variant)) stop("Mandatory argument variant is missing")

  if (variant == '') {
    return(NA)
  }

  var_vector <- strsplit(variant, ',')
  aa <- strsplit(var_vector[[1]], ':')[[1]][5]
  aa <- gsub(pattern='^p.', replacement='', x=aa)
  gene <- strsplit(var_vector[[1]], ':')[[1]][1]
  return(paste(sep='-', gene, aa))
}


get_gene_name <- function(variant=NULL) {
  if (is.null(variant)) stop("Mandatory argument variant is missing")
  var_vector <- strsplit(variant, ',')
  gene_name <- strsplit(var_vector[[1]], ':')[[1]][1]
  return(gene_name)
}


add_aa_to_data <- function(data=NULL) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data$aa <- apply(X=as.matrix(data$aa_change), MARGIN=1, FUN=get_aa_change)
  return(data)
}


create_recurrent_table <- function(data=NULL) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  df <- as.data.frame.matrix(table(data$aa, data$sample))
  df$rec <- apply(X=df, MARGIN=1, FUN=sum)
  df <- df[order(df$rec, decreasing=TRUE),]
  aa <- rownames(df)
  df <- cbind(aa, df)
  return(df)
}


# create a lookup table associating amino acid changes with a specific
# numeric value for the functional consequence
create_aa_lookup <- function() {
  func_cons <- c(
    'nonsynonymous SNV',
    'synonymous SNV',
    'frameshift deletion',
    'nonframeshift deletion',
    'frameshift insertion',
    'nonframeshift insertion',
    'stopgain',
    'stoploss',
    ''
    )
  func_cons_id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  lookup <- data.frame(matrix(ncol=2, nrow=length(func_cons)))
  colnames(lookup) <- c('func', 'id')
  lookup$func <- func_cons
  lookup$id <- func_cons_id

  return(lookup)
}


create_aa_table_count <- function(data=NULL) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data_tbl <- as.data.frame.matrix(table(data$aa, data$sample))
  return(data_tbl)
}


get_func_consequence_levels <- function() {
  return(c(
    'nonsynonymous SNV',
    'synonymous SNV',
    'frameshift deletion',
    'nonframeshift deletion',
    'frameshift insertion',
    'nonframeshift insertion',
    'stopgain',
    'stoploss'))
}


create_aa_dataframe <- function(path='.', threshold=2) {
  files <- get_variant_annotation_files(path=path)
  data <- get_variant_annotation_data(files=files)
  aa_df <- add_aa_to_data(data=data)
  # cleanup the dataframe by removing NAs and using ''
  #aa_df[is.na(aa_df$aa),]$aa <- 'none'
  # remove mutations that are not deleterious
  aa_df <- aa_df[!is.na(aa_df$aa),]

  # create a dataframe containing only those amino acid changes
  # with their corersponding functional consequence, this will
  # be used as a lookup table
  aa_lookup <- unique(aa_df[,c('aa', 'exonic_func')])

  # initialize the colour column
  func_ids <- create_aa_lookup()
  data_tbl <- create_recurrent_table(data=aa_df)

  # subset the data_tbl for aa and recurrence, then reduce recurrence
  # to greater than 4 samples
  data_df <- data.frame(aa=rownames(data_tbl), recurrence=data_tbl$rec)
  data_df <- data_df[data_df$recurrence >= threshold,]
  aa_levels <- as.character(data_df$aa)

  # create a geom_tile plottable dataframe containing aa, samples, exonic_func
  aa_df_del <- aa_df %>% filter(aa_df$aa %in% data_df$aa)
  aa_df_del <- aa_df_del %>% select(sample, aa, exonic_func)
  aa_df_del$aa <- factor(aa_df_del$aa, levels=rev(aa_levels))
  aa_df_del$exonic_func <- factor(aa_df_del$exonic_func, levels=get_func_consequence_levels())

  colnames(aa_df_del) <- c('sample', 'aa', 'Consequence')
  return(aa_df_del)
}


create_indel_heatmap <- function(data=NULL, base_size=12,
                                 type=c('frameshift insertion', 'frameshift deletion', 'nonframeshift insertion', 'nonframeshift deletion')) {
  if (is.null(data)) stop("Mandatory argument data is midding")

  data <- data %>% dplyr::filter(Consequence %in% type)
  plot <- ggplot(data, aes(sample, aa)) + geom_tile(aes(fill=Consequence)) + default.heatmap.theme(base_size=base_size) + xlab('') + ylab('')
  return(plot)
}


create_snv_heatmap <- function(data=NULL, base_size=12,
                               type=c('nonsynonymous SNV', 'synonymous SNV', 'stopgain', 'stoploss')) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data <- data %>% dplyr::filter(Consequence %in% type)
  plot <- ggplot(data, aes(sample, aa)) + geom_tile(aes(fill=Consequence)) + default.heatmap.theme(base_size=base_size) + xlab('') + ylab('')
}



create_mutation_heatmap <- function(data=NULL, base_size=12) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  plot <- ggplot(data, aes(sample, aa)) + geom_tile(aes(fill=Consequence)) + default.heatmap.theme(base_size=base_size) + xlab('') + ylab('')
}


# Main
parser <- ArgumentParser(description="Plot recurrent amino acid changes as a heatmap")
parser$add_argument('--path', '-p', default="./",
                    help="path to the *_multianno.txt files from ANOVAR")
parser$add_argument("--output", "-o", default="aa_mutation_heatmap.pdf",
                    help="path and filename to write heatmap to (default: aa_mutation_heatmap.pdf")
parser$add_argument("--threshold", "-t", default=2,
                    help="threshold value for recurrent mutations to plot (default: 2)")
parser$add_argument("--base_size", "-b", default=10,
                    help="base font size for the plot (default: 10)")
args <- parser$parse_args()

aa_df <- create_aa_dataframe(path=args$path, threshold=args$threshold)
num_samples <- length(unique(aa_df$sample))
num_aa <- length(unique(aa_df$aa))
pdf_height <- 0.075 * num_aa
pdf_width <- 0.075 * num_samples
pdf_height <- max(11, pdf_height)
pdf_width <- max(8.5, pdf_width)
aa_heatmap <- create_mutation_heatmap(data=aa_df, base_size=10) + scale_fill_brewer(palette="Set1")
ggsave(filename=args$output, plot=aa_heatmap, device='pdf', units="in", width=pdf_width, height=pdf_height)
