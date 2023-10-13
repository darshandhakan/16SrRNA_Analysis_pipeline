#!/usr/bin/env Rscript
################################################################################
### R script to compare two different conditions with count file (generated using corset and salmon) and DESeq2 package
### Aditya Narayan Sarangi
### Designed to be executed with bulkRNASeqPIPE
################################################################################

rm(list=ls())                                        # remove all the objects from the R session
suppressMessages(library(MicrobiotaProcess))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(gghalves))
suppressMessages(library(ggtreeExtra))
suppressMessages(library(ggtree))
suppressMessages(library(optparse))
suppressMessages(library(tibble))
suppressMessages(library(scales))
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))
suppressMessages(library(microbiomeMarker))
suppressMessages(library(RColorBrewer))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))


# to run the script in command lines

# options list with associated default value.
option_list <- list( 
  
  make_option(c("-p", "--project"),
              dest="project",
              default="Figures_and_Tables",
              help="Project Name to store images and tables [default: %default]."),
  
  make_option(c("-m", "--metadata"),
              dest="metadata",
              help="path to the QIIME2 mapping file [default: %default]."),
  
  make_option(c("-o", "--abundance"),
              dest="abundance",
              help="path to the QIIME OTU table in biom format [default: %default]."),
  
  make_option(c("-t", "--tree"),
              dest="tree",
              help="path to the tree.nwk file [default: %default]."),
  
  make_option(c("-g", "--groups"),
              dest="groups",
              help="Give the names of the groups to be compared")
)	


# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list, 
                       description="Genetate Images using microeco.",
                       epilogue="For comments, bug reports etc... please contact Aditya Narayan Sarangi <aditya.sarangi@basesolve.com>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options


#Check mandetory inputs 
if ( is.null(opt$metadata) ) {
  stop("--metadata tsv file in QIIME2 format must be provided. See script usage (--help)")
}

if ( is.null(opt$abundance) ) {
  stop("--abundance QIIME2 filtered .qza file must be provided. See script usage (--help)")
}

if ( is.null(opt$tree) ) {
  stop("--taxonmy path to silva taxonomy .qza file must be provided. See script usage (--help)")
}

# get options and arguments
workDir <- getwd()
metadata <- opt$metadata  
biom <- opt$abundance  
tree <- opt$tree
projectName <- opt$project

################################################################################
###                             running script                               ###
################################################################################



Plot_16S_rRNA <- function(biom, tree, metadata, projectName)
{
  
  ## Load the package:
  library(phyloseq)
  suppressMessages(library(MicrobiotaProcess))
  suppressMessages(library(ggplot2))
  suppressMessages(library(dplyr))
  suppressMessages(library(gghalves))
  suppressMessages(library(ggtreeExtra))
  suppressMessages(library(ggtree))
  suppressMessages(library(optparse))
  suppressMessages(library(tibble))
  suppressMessages(library(scales))
  suppressMessages(library(xlsx))
  suppressMessages(library(ampvis2))
  
  
  # - titleStyle : style object to use for title
  xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
    rows <-createRow(sheet,rowIndex=rowIndex)
    sheetTitle <-createCell(rows, colIndex=1)
    setCellValue(sheetTitle[[1,1]], title)
    setCellStyle(sheetTitle[[1,1]], titleStyle)
  }
  
  
  ## import the Qiime output:
  dataset <- import_biom(BIOMfilename = biom, treefilename = tree)
  SAMPLE_TABLE <- read.csv(file = metadata, header = TRUE, sep = "\t", row.names = 1)
  sam1 <- sample_data(SAMPLE_TABLE)
  dataset <- merge_phyloseq(dataset, sam1)
  colnames(tax_table(dataset)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  #dataset@taxatree@phylo$node.label <- gsub(x = dataset@taxatree@phylo$node.label, pattern = "^(.)__(un)_(.)__(.*)", replacement = "\\1__\\2_\\3__\\4_Unknown")
  
  ## Generate Ampvis data structure:
  
  metadata_table <- read.delim(file = metadata, header = T, sep = "\t")
  metadata_table$sample.name <- factor(metadata_table$sample.name, levels = metadata_table$sample.name)
  
  data <- amp_load(otutable = biom, metadata = metadata_table)
  
  
  
  ########## DIVERSITY ANALYSIS ############
  
  mpse <- dataset %>% as.MPSE()
  mpse %<>% mp_rrarefy()
  mpse %<>% mp_cal_rarecurve(.abundance=RareAbundance, chunks=10, action="add")
  mpse %<>% mp_cal_alpha(.abundance=RareAbundance)
  mpse2 <- mpse
  
  plot_types <- c("Bar_plot", "Alpha_rarefaction_curve", "Hierarchical_cluster", "Alpha_Diversity", "Beta_Diversity", "Heatmap")
  for (plot in plot_types){
    dir.create(file.path(projectName,"figures",plot),showWarnings = FALSE, recursive = TRUE)
  }
  
  dir.create(file.path(projectName, "tables"), showWarnings = FALSE, recursive = TRUE)
  
  
  #### Plot Rarefaction Curves ####
  
  p1 <- mpse %>% 
    mp_plot_rarecurve(
      .rare = RareAbundanceRarecurve, 
      .alpha = Observe,
      .group = sample.name
    ) +
    theme(axis.text=element_text(angle = 90, vjust = 0.5, hjust=1, size=8), panel.grid=element_blank(),
          strip.background = element_rect(colour=NA,fill="grey"),
          strip.text.x = element_text(face="bold"))
  
  rareplot_group <- p1
  rownames(mpse2@colData) <- SAMPLE_TABLE$Group
  ggsave(file=paste0(projectName,"/figures/","Alpha_rarefaction_curve/","rarefaction_curve_group.svg"), rareplot_group, width = 20, dpi = 300, units = "cm", device='svg')
  ggsave(file=paste0(projectName,"/figures/","Alpha_rarefaction_curve/","rarefaction_curve_group.png"), rareplot_group, width = 20, dpi = 300, units = "cm", device='png')
  
  ####### Plot the Alpha Diversity Metrics ########
  
  rownames(mpse2@colData) <- SAMPLE_TABLE$Group
  p2 <- mp_plot_alpha(.data = mpse2, .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)) +
  ylab(label = "ALPHA INDEX VALUES") + theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), strip.text.x.top = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"))
  ggsave(file=paste0(projectName,"/figures/","Alpha_Diversity/","Individual_alpha_diversity.svg"), p2, width = 20, height = 15, dpi = 300, units = "cm", device='svg')
  ggsave(file=paste0(projectName,"/figures/","Alpha_Diversity/","Individual_alpha_diversity.png"), p2, width = 20, height = 15, dpi = 300, units = "cm", device='png')
  
  
  ##### Plot the Beta Diversity Plots ######
  B1 <- mpse %<>% mp_cal_dist(.abundance=Abundance, distmethod="bray")
  betadiv <- B1 %>% mp_plot_dist(.distmethod = bray, .group = sample.name, group.test=TRUE, textsize=2) + 
  theme(axis.text.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"))
  ggsave(file=paste0(projectName,"/figures/","Beta_Diversity/","beta_diversity.svg"), betadiv, width = 20, dpi = 300, units = "cm", device='svg')
  ggsave(file=paste0(projectName,"/figures/","Beta_Diversity/","beta_diversity.png"), betadiv, width = 20, dpi = 300, units = "cm", device='png')
  
  
  ####### Plotted the Alpha Diversity Metrics #########
   
  ############################## TABLE ##############################
  ## Generate the Tables of Alpha Diversity and Beta Diversity ##
  
  wb <- createWorkbook(type = "xlsx")
  TITLE_STYLE3 <- CellStyle(wb)+ Font(wb,  heightInPoints=16,
                                       color="dodgerblue4", isBold=TRUE, underline=1)
  SUB_TITLE_STYLE3 <- CellStyle(wb) +
    Font(wb,  heightInPoints=14,
         isItalic=TRUE, isBold=FALSE)
  # Styles for the data table row/column names
  TABLE_ROWNAMES_STYLE3 <- CellStyle(wb) + Font(wb, isBold=TRUE)
  TABLE_COLNAMES_STYLE3 <- CellStyle(wb) + Font(wb, isBold=TRUE) +
    Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
    Border(color="black", position=c("TOP", "BOTTOM"),
           pen=c("BORDER_THIN", "BORDER_THICK"))
  
  ## Generate the Tables of Taxa Abundance
  
  wb2<-createWorkbook(type="xlsx")
  TITLE_STYLE2 <- CellStyle(wb2)+ Font(wb2,  heightInPoints=16,
                                       color="dodgerblue4", isBold=TRUE, underline=1)
  SUB_TITLE_STYLE2 <- CellStyle(wb2) +
    Font(wb2,  heightInPoints=14,
         isItalic=TRUE, isBold=FALSE)
  # Styles for the data table row/column names
  TABLE_ROWNAMES_STYLE2 <- CellStyle(wb2) + Font(wb2, isBold=TRUE)
  TABLE_COLNAMES_STYLE2 <- CellStyle(wb2) + Font(wb2, isBold=TRUE) +
    Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
    Border(color="black", position=c("TOP", "BOTTOM"),
           pen=c("BORDER_THIN", "BORDER_THICK"))
  
  #######################################################################
  
  
  
  ALPHA_DIV_TABLE <-  data.frame(mpse@colData@listData$sample.name, mpse@colData@listData$sample.name, mpse@colData@listData$Observe, mpse@colData@listData$Chao1, mpse@colData@listData$ACE, mpse@colData@listData$Shannon, mpse@colData@listData$Simpson, mpse@colData@listData$Pielou)
  colnames(ALPHA_DIV_TABLE) <- c("sample.name", "Group", "Observe", "Chao1", "ACE", "Shannon", "Simpson", "Pielou")
  sheet <- createSheet(wb, sheetName = paste0("ALPHA_DIVERSITY_TABLE"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("ALPHA DIVERSITY METRICS TABLE"), titleStyle = TITLE_STYLE3)
  addDataFrame(ALPHA_DIV_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE3, rownamesStyle = TABLE_ROWNAMES_STYLE3)
  
  ## Bray Curtis Table ##
  BRAY_TABLE <- data.frame()
  BRAY_STATS <- data.frame()
  BETWEEN_SAMPLE_DISTANCE <- c()
  WITHIN_SAMPLE_DISTANCE <- c()
  
  for (b in seq_along(B1@colData@listData$bray))
  {
    BRAY <- unlist(c(rep(NA, length(B1@colData@listData$bray) - nrow(B1@colData@listData$bray[[b]])), B1@colData@listData$bray[[b]][,"bray"]))
    BRAY_TABLE <- rbind(BRAY_TABLE, BRAY)
    BETWEEN_SAMPLE_DISTANCE <- c(BETWEEN_SAMPLE_DISTANCE, BRAY[B1@colData@listData$Group[b] != B1@colData@listData$Group])
    WITHIN_SAMPLE_DISTANCE <- c(WITHIN_SAMPLE_DISTANCE, BRAY[B1@colData@listData$Group[b] == B1@colData@listData$Group])
  }
  colnames(BRAY_TABLE) <- B1@colData@listData$bray[[1]]$braySampley
  rownames(BRAY_TABLE) <- B1@colData@listData$bray[[1]]$braySampley
  sheet <- createSheet(wb, sheetName = paste0("BETA_DIVERSITY_MATRIX"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("BETA DIVERSITY (BRAY CURTIS) DISTANCES ACROSS SAMPLES"), titleStyle = TITLE_STYLE3)
  addDataFrame(BRAY_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE3, rownamesStyle = TABLE_ROWNAMES_STYLE3)
  rownames(BRAY_TABLE) <- SAMPLE_TABLE$sample.name
  colnames(BRAY_TABLE) <- SAMPLE_TABLE$sample.name
  ## Save the Workbook of Diversity Metrics ##
  saveWorkbook(wb, paste0(projectName, "/tables/DIVERSITY_ANALYSIS_TABLE.xlsx"))
  
  cat ("##########\n\nCompleted Alpha and Beta Diversity Analysis\n\n##########")
  
  #### Plot PCoA plot of the samples ######
  
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(ncol(BRAY_TABLE))
  res.pca <- PCA(BRAY_TABLE,  graph = FALSE)
  
  PCA_PLOT <- fviz_pca_ind(res.pca, col.ind = colnames(BRAY_TABLE),repel = TRUE, pointsize = 3, labelsize = 6) + scale_color_manual(name = "Samples", values = mycolors) + scale_shape(name = "Samples") + 
    labs(title = "PCoA Plot of samples", x = paste0("PC1 [",signif(res.pca$eig[1,2], digits = 4), "] % Variation explained"), y = paste0("PC2 [", signif(res.pca$eig[2,2], digits = 4), "] % Variation explained")) + 
    theme(axis.text =  element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"))
  
  ggsave(file=paste0(projectName,"/figures/","Beta_Diversity/","PCOA.png"), PCA_PLOT, width = 30, height = 15, dpi = 300, units = "cm", device='png', bg = "white")
  ggsave(file=paste0(projectName,"/figures/","Beta_Diversity/","PCOA.svg"), PCA_PLOT, width = 30, height = 15, dpi = 300, units = "cm", device='svg', bg = "white")
  
  cat ("##########\n\nCompleted PCoA Analysis##########")
  
  ### Perform Hierarchical clustering of the samples ####
  
  #The hierarchical cluster result of samples
  Clust <- mpse %<>%
    mp_cal_clust(
      .abundance = Abundance, 
      distmethod = "bray",
      hclustmethod = "average", # (UPGAE)
      action = "add" # action is used to control which result will be returned
    )
  
  sample.clust <- Clust %>% mp_extract_internal_attr(name='SampleClust')
  
  library(ggtree)
  p <- ggtree(sample.clust) + 
    geom_tippoint(aes(color=Group), size = 4) +
    geom_tiplab(aes(label = sample.name), size = 20, as_ylab = TRUE) +
    ggplot2::scale_x_continuous(expand=c(0, 0.01)) + theme(axis.text = element_text(size = 12, face = "bold"), legend.position = "none")
  
  ggsave(file=paste0(projectName,"/figures/","Hierarchical_cluster/","hierarchical_cluster_samples.png"), p, width = 30, dpi = 300, units = "cm", device='png')
  ggsave(file=paste0(projectName,"/figures/","Hierarchical_cluster/","hierarchical_cluster_samples.svg"), p, width = 30, dpi = 300, units = "cm", device='svg')
  
  cat ("##########\n\nCompleted Sample Clustering\n\n##########")
  
  ####
  #--------------------------------------------------------@
  
  ######### PLOT TAXA-WISE BARPLOTS ##################
  
  RareAbund <- mpse %<>%
    mp_cal_abundance( 
      .abundance = RareAbundance
    ) %>%
    mp_cal_abundance( 
      .abundance=RareAbundance,
      .group=Group
    )
  
  my_colors <- c(microshades_palette("micro_cvd_purple", 5, lightest = TRUE), microshades_palette("micro_cvd_blue", 5, lightest = TRUE), microshades_palette("micro_cvd_green", 5, lightest = TRUE), microshades_palette("micro_cvd_gray", 5, lightest = TRUE), microshades_palette("micro_cvd_orange", 5, lightest = TRUE))
  COL <- sample(x = my_colors, size = 25, replace = FALSE)
  
  ##### Phylum ######
  
  ## Get the Phylum  
  phyla.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Phylum)
  bar <- phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phylum="label")
  BAR <- bar[c("Phylum","Sample","RareAbundance")]
  BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
  sheet <- createSheet(wb2, sheetName = paste0("Phylum_Abundance_Table"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("Phylum Abundance Table"), titleStyle = TITLE_STYLE2)
  addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
  
  phyla.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Phylum, topn=20)
  phyla.tb$label <- gsub(x = phyla.tb$label, pattern = ".*([a-z]__)", "\\1")
  bar <- phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phylum="label")
  p1 <- ggplot(bar, aes(x=paste0(Sample), y=RelRareAbundanceBySample, fill = Phylum)) + geom_bar(stat = "identity", position = "stack", width = 0.8) + coord_flip() + scale_fill_manual(values = COL) +  ylab(label = "Relative Abundance of Phylum") + xlab(label = "Samples") + 
    theme(axis.text.y = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 12, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  
  p1 <- p + geom_fruit(data=bar, geom=geom_col, mapping = aes(x = RelRareAbundanceBySample, y = Sample, fill = Phylum),
                       orientation = "y", pwidth = 3, 
                       axis.params = list(axis = "x", text.size = 2, limits = c(0,120), vjust = 1),
                       grid.params = list()) + scale_fill_manual(values = mycolors)
  p1 <- p1 + theme(legend.position = "bottom", axis.text.y = element_text(size = 12, face = "bold"), legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14, face = "bold"))  
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Phylum_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Phylum_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
  
  cat ("##########\n\nPhylum Bar plots completed\n\n##########")
  
  
  ####### Class ########
  
  class.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Class)
  bar <- class.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Class="label")
  ## Add taxa Table
  BAR <- bar[c("Class","Sample","RareAbundance")]
  BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
  sheet <- createSheet(wb2, sheetName = paste0("Class_Abundance_Table"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("Class Abundance Table"), titleStyle = TITLE_STYLE2)
  addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
  
  # Plot the Class Taxa.
  class.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Class, topn = 20)
  class.tb$label <- gsub(x = class.tb$label, pattern = ".*([a-z]__)", "\\1")
  bar <- class.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Class="label")
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(21)
  p1 <- p + geom_fruit(data=bar, geom=geom_col, mapping = aes(x = RelRareAbundanceBySample, y = Sample, fill = Class),
                       orientation = "y", pwidth = 3, 
                       axis.params = list(axis = "x", text.size = 2, limits = c(0,120), vjust = 1),
                       grid.params = list()) + scale_fill_manual(values = mycolors)
  p1 <- p1 + theme(legend.position = "bottom", axis.text.y = element_text(size = 12, face = "bold"), legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14, face = "bold"))  
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Class_barplots.png"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='png')
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Class_barplots.svg"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='svg')

  cat ("##########\n\nClass Bar Plots Completed\n\n##########")
  
  ########## Order ############
  
  order.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Order)
  bar <- order.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Order="label")
  ## Add taxa Table
  BAR <- bar[c("Order","Sample","RareAbundance")]
  BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
  sheet <- createSheet(wb2, sheetName = paste0("Order_Abundance_Table"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("Order Abundance Table"), titleStyle = TITLE_STYLE2)
  addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
  
  # Plot the Class Taxa.
  order.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Order, topn = 20)
  order.tb$label <- gsub(x = order.tb$label, pattern = ".*([a-z]__)", "\\1")
  bar <- order.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Order="label")
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(21)
  p1 <- p + geom_fruit(data=bar, geom=geom_col, mapping = aes(x = RelRareAbundanceBySample, y = Sample, fill = Order),
                       orientation = "y", pwidth = 3, 
                       axis.params = list(axis = "x", text.size = 2, limits = c(0,120), vjust = 1),
                       grid.params = list()) + scale_fill_manual(values = mycolors)
  p1 <- p1 + theme(legend.position = "bottom", axis.text.y = element_text(size = 12, face = "bold"), legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14, face = "bold"))  
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Order_barplots.png"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='png')
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Order_barplots.svg"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='svg')
  
  cat ("##########\n\nOrder Bar Plots Completed\n\n##########")
  
  ######### Family ############
  
  family.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Family)
  bar <- family.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Family="label")
  ## Add taxa Table
  BAR <- bar[c("Family","Sample","RareAbundance")]
  BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
  sheet <- createSheet(wb2, sheetName = paste0("Family_Abundance_Table"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("Family Abundance Table"), titleStyle = TITLE_STYLE2)
  addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
  
  # Plot the Family Taxa.
  family.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Family)
  strings <- c("o__", "c__", "p__")
  family.tb <- family.tb[!str_detect(string = family.tb$label, pattern = paste(strings, collapse = '|')), ]
  family.tb$label <- gsub(x = family.tb$label, pattern = ".*([a-z]__)", "\\1")
  bar <- family.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Family="label")
  bar <- bar[rev(order(bar$RareAbundance)),]
  Family <- unique(bar$Family)[1:20]
  bar <- bar[bar$Family %in% Family,]
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(21)
  p1 <- p + geom_fruit(data=bar, geom=geom_col, mapping = aes(x = RelRareAbundanceBySample, y = Sample, fill = Family),
                       orientation = "y", pwidth = 3, 
                       axis.params = list(axis = "x", text.size = 2, limits = c(0,120), vjust = 1),
                       grid.params = list()) + scale_fill_manual(values = mycolors)
  p1 <- p1 + theme(legend.position = "bottom", axis.text.y = element_text(size = 12, face = "bold"), legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14, face = "bold"))
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Family_barplots.png"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='png')
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Family_barplots.svg"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='svg')
  
  cat ("##########\n\nFamily Bar Plots completed\n\n##########")
  
  ############# Genus #################
  
  genus.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Genus)
  bar <- genus.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Genus="label")
  ## Add taxa Table
  BAR <- bar[c("Genus","Sample","RareAbundance")]
  BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
  sheet <- createSheet(wb2, sheetName = paste0("Genus_Abundance_Table"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("Genus Abundance Table"), titleStyle = TITLE_STYLE2)
  addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
  
  # Plot the Class Taxa.
  genus.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Genus)
  strings <- c("f__", "o__", "c__", "p__")
  genus.tb <- genus.tb[!str_detect(string = genus.tb$label, pattern = paste(strings, collapse = '|')), ]
  genus.tb$label <- gsub(x = genus.tb$label, pattern = ".*([a-z]__)", "\\1")
  bar <- genus.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Genus="label")
  bar <- bar[rev(order(bar$RareAbundance)),]
  Genus <- unique(bar$Genus)[1:20]
  bar <- bar[bar$Genus %in% Genus,]
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(21)
  p1 <- p + geom_fruit(data=bar, geom=geom_col, mapping = aes(x = RelRareAbundanceBySample, y = Sample, fill = Genus),
                       orientation = "y", pwidth = 3, 
                       axis.params = list(axis = "x", text.size = 2, limits = c(0,120), vjust = 1),
                       grid.params = list()) + scale_fill_manual(values = mycolors)
  p1 <- p1 + theme(legend.position = "bottom", axis.text.y = element_text(size = 12, face = "bold"), legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14, face = "bold"))  
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Genus_barplots.png"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='png')
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Genus_barplots.svg"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='svg')
  
  cat ("##########\n\nGenus Bar Plots completed\n\n##########")
  
  ############# Species #################
  
  species.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Species)
  bar <- species.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Species="label")
  ## Add taxa Table
  BAR <- bar[c("Species","Sample","RareAbundance")]
  BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
  sheet <- createSheet(wb2, sheetName = paste0("Species_Abundance_Table"))
  xlsx.addTitle(sheet, rowIndex=1, title=paste0("Species Abundance Table"), titleStyle = TITLE_STYLE2)
  addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
  
  # Plot the Species Taxa #
  species.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Species)
  strings <- c("f__", "o__", "c__", "p__", "g__")
  species.tb <- species.tb[!str_detect(string = species.tb$label, pattern = paste(strings, collapse = '|')), ]
  species.tb$label <- gsub(x = species.tb$label, pattern = ".*([a-z]__)", "\\1")
  bar <- species.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Species="label")
  bar <- bar[rev(order(bar$RareAbundance)),]
  Species <- unique(bar$Species)[1:20]
  bar <- bar[bar$Species %in% Species,]
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(21)
  p1 <- p + geom_fruit(data=bar, geom=geom_col, mapping = aes(x = RelRareAbundanceBySample, y = Sample, fill = Species),
                       orientation = "y", pwidth = 3, 
                       axis.params = list(axis = "x", text.size = 2, limits = c(0,120), vjust = 1),
                       grid.params = list()) + scale_fill_manual(values = mycolors)
  p1 <- p1 + theme(legend.position = "bottom", axis.text.y = element_text(size = 12, face = "bold"), legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20, face = "bold"))  
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Species_barplots.png"), p1, width = 100, height = 60, dpi = 300, units = "cm", device='png')
  ggsave(file=paste0(projectName,"/figures/","Bar_plot/","Species_barplots.svg"), p1, width = 100, height = 60, dpi = 300, units = "cm", device='svg')

  cat ("##########\n\nSpecies Bar Plots Completed\n\n##########")
  
  cat ("##########\n\nGenerate the Excel Tables\n\n##########")
  saveWorkbook(wb2, paste0(projectName, "/tables/TAXA_ABUNDANCE_TABLE.xlsx"))
  
  ## Generate the Heatmaps using ampvis2
  
  ## Heatmap ##
  lineage <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (t in seq_along(lineage))
  {
    plot_heatmap <- amp_heatmap(data, tax_aggregate = lineage[t], group_by = "Group", tax_show = 50, min_abundance = 0.001,
                                facet="sample.name", round=5, tax_empty = "remove", plot_values_size = 3, normalise = TRUE, plot_values = TRUE) +
      theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1), axis.text.y = element_text(size=10), legend.position="right")
    
    ggsave(file=paste0(projectName,"/figures/","Heatmap/",lineage[t],"_heatmap.svg"), plot_heatmap, width = 15, dpi = 300, units = "cm", device='svg')
    ggsave(file=paste0(projectName,"/figures/","Heatmap/",lineage[t],"_heatmap.png"), plot_heatmap, width = 15, dpi = 300, units = "cm", device='png')
    ggsave(file=paste0(projectName,"/figures/","Heatmap/",lineage[t],"_heatmap.pdf"), plot_heatmap, width = 15, dpi = 300, units = "cm", device='pdf')
  }
  
  
  cat ("##########----------##########\n\n16SrRNA Pipeline Executed\n\n##########----------##########")
}

suppressWarnings(Plot_16S_rRNA(biom = biom, tree = tree, metadata = metadata, projectName = projectName))
