### plot_amplicon_coverage_heatmap.R ##############################################################
#  A tool to plot the amplicon coverage across samples as part of the QC analysis for the
#  nCoV project.

### HISTORY #######################################################################################
# Version           Date            Developer               Comments
# 0.01              2020-06-16      Richard J. de Borja     Initial development
# 0.02              2020-07-09      Richard J. de Borja     Adjusted PDF page height for large
#                                                           sample numbers
# 0.03              2020-08-27      Richard J. de Borja     Added output for table of amplicon
#                                                           coverage
### PREAMBLE ######################################################################################
library(argparse)

usage <- function() {
  usage.text <- '\nUsage: script.R --path <full path to coverage files>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("-p", "--path", default=".", help="full path to the coverage files")
parser$add_argument("-o", "--output", default="amplicon_coverage_heatmap.pdf",
                    help="filename to use for plot output")
parser$add_argument("-t", "--table", default="amplicon_coverage_table.tsv",
                    help="filename to use for table output")
args <- parser$parse_args()

# Set the environment
library(lattice)
library(dplyr)
library(tidyr)


### FUNCTIONS #####################################################################################
import_amplicon_coverage <- function(file=NULL, sample=NULL) {
  if (is.null(file)) stop("Mandatory argument file is missing")
  if (is.null(sample)) stop("Mandatory argument sample is missing")

  data <- read.table(
    file=file,
    header=TRUE,
    as.is=TRUE,
    sep='\t',
    quote='\"'
  )
  data$sample <- sample

  return(data)
}


get_sample_amplicon_coverage <- function(data=NULL, sample_name=NULL, feature_column=c('amplicon_id', 'read_count')) {
  if (is.null(data)) stop("Mandatory argument data is missing")
  if (is.null(sample_name)) stop("Mandatory argument sample_name is missing")
  x <- as.data.frame(data[,c(feature_column)])
  colnames(x) <- c(sample_name)
  return(x)
}


get_amplicon_coverage_files <- function(path='.', pattern='.amplicon_coverage.bed') {
  amp_cov_files <- list.files(
    path=path,
    pattern=pattern,
    recursive=TRUE,
    full.names=TRUE
  )
  amp_cov_df <- data.frame(matrix(ncol=2, nrow=length(amp_cov_files)))
  colnames(amp_cov_df) <- c('sample', 'file')
  for(i in 1:length(amp_cov_files)) {
    sample_name <- gsub(pattern=pattern, replacement='', x=basename(amp_cov_files[i]))
    amp_cov_df$sample[i] <- sample_name
    amp_cov_df$file[i] <- amp_cov_files[i]
  }
  return(amp_cov_df)
}

get_amplicon_depth_files <- function(path='.', pattern='.amplicon_depth.bed') {
  files <- list.files(
    path=path,
    pattern=pattern,
    recursive=TRUE,
    full.names=TRUE
  )
  df <- data.frame(matrix(ncol=2, nrow=length(files)))
  colnames(df) <- c('sample', 'file')
  for(i in 1:length(files)) {
    sample_name <- gsub(pattern=pattern, replacement='', x=basename(files[i]))
    df$sample[i] <- sample_name
    df$file[i] <- files[i]
  }
  return(df)
}

### GET DATA ######################################################################################
path <- args$path
amp_cov_files <- get_amplicon_depth_files(path=path)
cov_data <- data.frame()
for(i in 1:nrow(amp_cov_files)) {
  tmp_data <- import_amplicon_coverage(file=amp_cov_files$file[i], sample=amp_cov_files$sample[i])
  tmp_data_df <- tmp_data %>% dplyr::select(amplicon_id, sample, mean_depth)
  tmp_data_df$amplicon_id <- paste(
    sep='', 'amp',
    gsub(x=tmp_data_df$amplicon_id, pattern='nCoV-2019_', replacement=''))
  cov_data <- rbind(cov_data, tmp_data_df)
}

### PROCESS DATA ##################################################################################
cov_data_df <- as.data.frame(cov_data %>% pivot_wider(names_from=amplicon_id, values_from=mean_depth))
write.table(
  x=cov_data_df,
  file=args$table,
  quote=FALSE,
  row.names=FALSE,
  sep='\t'
)
rownames(cov_data_df) <- cov_data_df$sample
cov_data_df <- cov_data_df[,2:ncol(cov_data_df)]
num_samples <- nrow(cov_data_df)
print(paste(sep=' ', 'number of samples: ', num_samples))

# add a pseudo count of 0.01 to fix the log10 0 issue
cov_data_df <- cov_data_df + 0.01
cov_data_df <- log10(cov_data_df)

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################
cols <- colorRampPalette(c('#E41A1C', 'white', '#377EB8'))
pdf_height <- 0.075 * num_samples
pdf_height <- max(8, pdf_height)
pdf(file=args$output, height=pdf_height, width=10)
lattice::levelplot(
  x = as.matrix(t(cov_data_df)),
  col.regions=cols,
  scales=list(
    x=list(rot=90),
    cex=c(0.45,0.45)
  ),
  colorkey=list(space="top"),
  xlab="",
  ylab=""
)
dev.off()

