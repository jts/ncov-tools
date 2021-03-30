library(dplyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-f", "--file", help="Amplicon coverage file to process")
parser$add_argument("-t", "--threshold", default=10.0,
                    help="Threshold cutoff for coverage (default: 10.0)")
parser$add_argument("-o", "--output", help="Filename for PDF")
args <- parser$parse_args()

default.barplot.theme <- function (base_size = 24, base_family = "Helvetica", xaxis.angle = 90) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text.x = element_text(angle = xaxis.angle, 
                                     size = rel(0.6), hjust = 1))
}

file = args$file
data <- read.table(file=file, header=TRUE, as.is=TRUE, sep='\t', quote='\"')
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

legend_title <- paste(sep=' ', "Coverage <", cutoff)
pdf_height <- 7
pdf_width <- 10
plot <- df %>% ggplot(aes(x=amplicon, y=coverage, fill=low_coverage)) +
  geom_bar(stat='identity') +
  default.barplot.theme(base_size=12) +
  xlab('') +
  ylab('Average Coverage') +
  scale_y_log10() +
  scale_fill_manual(legend_title, values=c("#377EB8", "#E41A1C"))
ggsave(filename=args$output, plot=plot, device='pdf', width=pdf_width, height=pdf_height)