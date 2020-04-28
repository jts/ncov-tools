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
plot_tree_with_snps <- function(tree, alleles)
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
  
  r <- g + panel.snps + plot_layout(widths = c(1, 3))
  return(r)
}

require(ggtree)
if(!interactive()) {
    args = commandArgs(trailingOnly=TRUE)

    # default to nextstrain structure
    data_dir <- "results"
    tree_fn <- "tree.nwk"

    if(length(args) == 1) {
        data_dir <- args[1]
        tree_fn <- "tree_raw.nwk"
    }

    tree_path = paste(data_dir, tree_fn, sep="/")
    tree <- read.tree(tree_path)

    alleles_path = paste(data_dir, "alleles.tsv", sep="/")
    alleles <- read.table(alleles_path, header=T)
    p <- plot_tree_with_snps(tree, alleles)

    plot_path = paste(data_dir, "tree_snps.pdf", sep="/")
    ggsave(plot_path, p, height=8, width=20)
}
