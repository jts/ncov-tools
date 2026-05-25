library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(argparse)

get_func_consequence_levels <- function() {
  return(c(
    'chromosome_number_variation',
    'exon_loss_variant',
    'frameshift_variant',
    'rare_amino_acid_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'start_lost',
    'stop_gained',
    'transcript_ablation',
    '3_prime_UTR_truncation',
    'exon_loss',
    '5_prime_UTR_truncation',
    'coding_sequence_variant',
    'conservative_inframe_deletion',
    'conservative_inframe_insertion',
    'disruptive_inframe_deletion',
    'disruptive_inframe_insertion',
    'missense_variant',
    'regulatory_region_ablation',
    'splice_region_variant',
    'TFBS_ablation',
    '5_prime_UTR_premature_start_codon_gain_variant',
    'initiator_codon_variant',
    'start_retained',
    'stop_retained_variant',
    'synonymous_variant',
    '3_prime_UTR_variant',
    '5_prime_UTR_variant',
    'conserved_intergenic_variant',
    'conserved_intron_variant',
    'downstream_gene_variant',
    'exon_variant',
    'feature_elongation',
    'feature_truncation',
    'gene_variant',
    'intergenic_region',
    'intragenic_variant',
    'intron_variant',
    'mature_miRNA_variant',
    'miRNA',
    'NMD_transcript_variant',
    'non_coding_transcript_exon_variant',
    'non_coding_transcript_variant',
    'regulatory_region_amplification',
    'regulatory_region_variant',
    'TF_binding_site_variant',
    'TFBS_amplification',
    'transcript_amplification',
    'transcript_variant',
    'upstream_gene_variant'))
}

get_func_consequence_levels_annovar <- function() {
  return(c(
    'frameshift insertion',
    'frameshift deletion',
    'frameshift block substitution',
    'stopgain',
    'stoploss',
    'nonframeshift insertion',
    'nonframeshift deletion',
    'nonframeshift block substitution',
    'nonsynonymous SNV',
    'synonymous SNV',
    'unknown')
  )
}


default.heatmap.theme <- function(base_size = 24, base_family = 'Helvetica', xaxis.angle = 90) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text.x = element_text(angle = xaxis.angle, size = rel(0.8), hjust = 1)
    )
}

create_mutation_heatmap <- function(data=NULL, base_size=12) {
  if (is.null(data)) stop("Mandatory argument data is missing")
  
  plot <- ggplot(data, aes(sample, aa)) + geom_tile(aes(fill=Consequence)) + default.heatmap.theme(base_size=base_size) + xlab('') + ylab('')
}

import_var_data <- function(path=NULL, pattern='aa_table.tsv', recursive=TRUE) {
  files = list.files(path=path, pattern=pattern, full.names=TRUE, recursive=recursive)
  data = data.frame()
  for(i in 1:length(files)) {
    tmp_data = read.table(file=files[i], header=TRUE, as.is=TRUE, sep='\t', quote='\"')
    data <- rbind(data, tmp_data)
  }
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

create_aa_dataframe <- function(path='.', threshold=2) {
  aa_df <- import_var_data(path=path)
  aa_df <- aa_df[!is.na(aa_df$aa),]
  
  aa_lookup <- unique(aa_df[,c('aa', 'Consequence')])
  data_tbl <- create_recurrent_table(data=aa_df)
  
  data_df <- data.frame(aa=rownames(data_tbl), recurrence=data_tbl$rec)
  data_df <- data_df[data_df$recurrence >= threshold,]
  aa_levels <- as.character(data_df$aa)
  
  aa_df_del <- aa_df %>% filter(aa_df$aa %in% data_df$aa)
  aa_df_del <- aa_df_del %>% select(sample, aa, Consequence)
  aa_df_del$aa <- factor(aa_df_del$aa, levels=rev(aa_levels))
#  aa_df_del$Consequence <- factor(aa_df_del$Consequence, levels=get_func_consequence_levels_annovar)
  aa_df_del$exonic_func <- factor(aa_df_del$Consequence, levels=get_func_consequence_levels())

  #colnames(aa_df_del) <- c('sample', 'aa', 'Consequence')
  return(aa_df_del)
}


#
# Main
#
description <- 'Plot recurrent amino acid changes as a heatmap'
parser <- ArgumentParser(description=description)
parser$add_argument('--path', '-p', default='./',
                    help='path to the annotated variant files')
parser$add_argument('--output', '-o', default='aa_mutation_heatmap.pdf')
parser$add_argument('--threshold', '-t', default=2,
                    help='threshold value for recurrent mutations to plot (default: 2)')
parser$add_argument('--base_size', '-b', default=10,
                    help='base font size for the plot (default: 10)')
args <- parser$parse_args()

data <- create_aa_dataframe(path=args$path, threshold=args$threshold)
num_samples <- length(unique(data$sample))
num_aa <- length(unique(data$aa))
pdf_height <- 0.075 * num_aa
pdf_width <- 0.075 * num_samples
pdf_height <- max(11, pdf_height)
pdf_width <- max(8.5, pdf_width)

aa_heatmap <- create_mutation_heatmap(data=data, base_size=8) + scale_fill_brewer(palette='Set1')
ggsave(filename=args$output, plot=aa_heatmap, device='pdf', units='in', width=pdf_width, height=pdf_height)
