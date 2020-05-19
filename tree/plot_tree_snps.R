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
plot_tree_with_snps <- function(tree, alleles, lineage_path)
{
  require(ggtree)
  require(ggplot2)
  require(patchwork)
  
  # get order of tips in the tree
  # from: https://groups.google.com/forum/#!topic/bioc-ggtree/LqRDK78m3U4
  d = fortify(tree)
  d = subset(d, isTip)
  tip.order = with(d, label[order(y, decreasing=F)])
  
  # order SNPs according to tip
  alleles$name = factor(alleles$name, levels=tip.order)
  alleles$pos = factor(alleles$pos)

  bases = c("A", "C", "G", "T", "N", "X")
  alleles$alt_allele = factor(alleles$alt_allele, bases)
  ambiguous_bases = !(alleles$alt_allele %in% bases)
  alleles$alt_allele[ambiguous_bases] = "X"
  
  # basic tree plot
  # uncomment to show tip labels, to make sure they are aligned correctly
  g <- ggtree(tree) + 
    #geom_tiplab(align=TRUE) + 
    #xlim(0, 0.00050) +
    scale_y_tree() +
    theme_tree2()
  
  # draw panel with variants
  cols <- c("blue", "red", "green3", "purple3", "lightgrey", "black")
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
    scale_fill_manual(name="Variant", values=cols, drop=FALSE)
  
  # if a lineage file exists, annotate lineages on the plot
  if(file.exists(lineage_path)) { 
     lineage <- read.csv(
          file=lineage_path,
          header=TRUE,
          as.is=TRUE,
          quote="\"")

      lineage <- lineage[c("taxon", "lineage")]
      colnames(lineage) <- c("name", "lineage")
      lineage$name <- factor(lineage$name, levels=tip.order)
      lineage$pos <- 1
      lineage$lineage <- as.factor(lineage$lineage)

      panel.cov <- ggplot(lineage, aes(x=pos, y=name)) +
        geom_tile(aes(fill=lineage), color='white') +
        geom_text(aes(label=lineage), size=3) +
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
  } else {
    r <- g + panel.snps + plot_layout(widths = c(1, 3))
  }
  return(r)
}

require(ggtree)
if(!interactive()) {
    args = commandArgs(trailingOnly=TRUE)

    if(length(args) < 3) {
        q()
    }

    plot_path <- args[1]
    tree_path <- args[2]
    alleles_path <- args[3]
    lineage_path = ""
    if(length(args) == 4) {
        lineage_path = args[4]
    }

    tree <- read.tree(tree_path)
    alleles <- read.table(alleles_path, header=T)
    p <- plot_tree_with_snps(tree, alleles, lineage_path)
    
    # count number of samples, for scaling the plot
    d = fortify(tree)
    d = subset(d, isTip)
    num_samples = nrow(d)
    pdf_height = 0.125 * num_samples
    pdf_height = max(8, pdf_height)
    ggsave(plot_path, p, height=pdf_height, width=20, limitsize=FALSE)
}
