#
# Plot a phylogenetic tree with associated mutations
#

#
# this method of aligning the tree and tiles is from:
# https://thackl.github.io/ggtree-composite-plots
#

# fix the alignment of the tree panel
scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}

# 
plot_tree_with_snps <- function(tree, alleles, lineage)
{
  require(ggtree)
  require(ggplot2)
  library(patchwork)
  
  # get order of tips in the tree
  # from: https://groups.google.com/forum/#!topic/bioc-ggtree/LqRDK78m3U4
  d = fortify(tree)
  d = subset(d, isTip)
  tip.order = with(d, label[order(y, decreasing=F)])
  
  # order SNPs according to tip
  alleles$name = factor(alleles$name, levels=tip.order)
  alleles$pos = factor(alleles$pos)
  alleles$alt_allele = factor(alleles$alt_allele, c("A", "C", "G", "T", "N"))
  
  # basic tree plot
  # uncomment to show tip labels, to make sure they are aligned correctly
  g <- ggtree(tree) + 
    #geom_tiplab(align=TRUE) + 
    #xlim(0, 0.00050) +
    scale_y_tree() +
    theme_tree2()
  
  # draw panel with variants
  cols <- c("blue", "red", "green3", "purple3", "lightgrey")
  panel.snps <- ggplot(alleles, aes(x=pos, y=name)) + 
    geom_tile(aes(fill=alt_allele), color="white") +
    ylim(tip.order) +
    theme_bw() + 
    theme(axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
                       axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                       axis.text.x = element_text(angle = 90, hjust = 1),
                       legend.position = "top") +
    scale_fill_manual(name="Variant", values=cols)
  
  lineage <- lineage[c("taxon", "lineage")]
  colnames(lineage) <- c("name", "lineage")
  lineage$name <- factor(lineage$name, levels=tip.order)
  lineage$pos <- 1
  lineage$lineage <- as.factor(lineage$lineage)

  panel.cov <- ggplot(lineage, aes(x=pos, y=name)) +
    geom_tile(aes(fill=lineage), color='white') +
    ylim(tip.order) +
    theme_bw() +
    theme(axis.line = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
          axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

  r <- g + panel.cov + panel.snps + plot_layout(widths = c(1, 0.1, 3), guides="collect")
  return(r)
}

require(ggtree)
if(!interactive()) {
    args = commandArgs(trailingOnly=TRUE)

    # default to nextstrain structure
    data_dir <- "results"
    tree_fn <- "tree.nwk"
    lineage_fn <- "lineage_report.csv"

    if(length(args) == 1) {
        data_dir <- args[1]
    }

    tree_path = paste(data_dir, tree_fn, sep="/")
    tree <- read.tree(tree_path)
    
    alleles_path = paste(data_dir, "alleles.tsv", sep="/")
    alleles <- read.table(alleles_path, header=T)

    lineage <- read.csv(
        file=lineage_fn,
        header=TRUE,
        as.is=TRUE,
        quote="\""
    )

    p <- plot_tree_with_snps(tree, alleles, lineage)

    plot_path = paste("tree_snps_lineage.pdf", sep="/")
    
    # count number of samples, for scaling the plot
    d = fortify(tree)
    d = subset(d, isTip)
    num_samples = nrow(d)
    pdf_height = 0.125 * num_samples
    pdf_height = max(8, pdf_height)
    ggsave(plot_path, p, height=pdf_height, width=20)
}
