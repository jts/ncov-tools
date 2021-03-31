library(tidyr)
library(dplyr)
library(ggplot2)
library(argparse)

init_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument("-f", "--file", help="Amplicon coverage file to process")
    parser$add_argument("-t", "--threshold", default=20.0, type='double',
                        help="Threshold cutoff for coverage (default: 20.0)")
    parser$add_argument("-o", "--output", help="Filename for PDF")
    return(parser$parse_args())
}

default.boxplot.theme <- function (base_size = 24, base_family = "Helvetica", xaxis.angle = 90) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text.x = element_text(angle = xaxis.angle, 
                                     size = rel(0.6), hjust = 1))
}

args <- init_args()
data <- read.table(file=args$file, header=TRUE, as.is=TRUE, sep='\t', quote='\"')
amps <- colnames(data[,2:ncol(data)])
data_df <- tidyr::gather(data[,2:ncol(data)])
data_df$key <- factor(data_df$key, levels=amps)
colnames(data_df) <- c('amplicon', 'coverage')

df <- data.frame(matrix(ncol=3, nrow=ncol(data)-1))
colnames(df) <- c('amplicon', 'coverage', 'low_coverage')
df$amplicon <- colnames(data)[2:ncol(data)]
for(i in 2:ncol(data)) {
  df$coverage[i-1] <- mean(data[,i])
}
amps <- df$amplicon
df$amplicon <- factor(df$amplicon, levels=amps)

cutoff <- args$threshold
for(i in 1:nrow(df)) {
  if (df$coverage[i] < cutoff) {
    df$low_coverage[i] <- TRUE
  } else {
    df$low_coverage[i] <- FALSE
  }
}
df <- df[,c('amplicon', 'low_coverage')]
data_df <- dplyr::left_join(x=data_df, y=df, by=c('amplicon'))

legend_title <- bquote(mu ~ covg ~ "<" ~ .(cutoff))

pdf_height <- 7
pdf_width <- 10
plot <- data_df %>% ggplot(aes(x=amplicon, y=coverage, fill=low_coverage)) +
  geom_boxplot() +
  scale_y_log10() +
  xlab('') +
  ylab(bquote(mu ~ 'Coverage per Sample')) +
  scale_fill_manual(legend_title, values=c("#377EB8", "#E41A1C")) +
  default.boxplot.theme(base_size=12)
ggsave(filename=args$output, plot=plot, device='pdf', width=pdf_width, height=pdf_height)
