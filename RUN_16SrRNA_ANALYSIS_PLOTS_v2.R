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

  ## import the Qiime output:
  dataset <- import_biom(BIOMfilename = biom, treefilename = tree)
  SAMPLE_TABLE <- read.csv(file = metadata, header = TRUE, sep = "\t", row.names = 1)
  sam1 <- sample_data(SAMPLE_TABLE)
  dataset <- merge_phyloseq(dataset, sam1)
  colnames(tax_table(dataset)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  #dataset@taxatree@phylo$node.label <- gsub(x = dataset@taxatree@phylo$node.label, pattern = "^(.)__(un)_(.)__(.*)", replacement = "\\1__\\2_\\3__\\4_Unknown")
  
  
  ## First perform the individual Alpha Rarefaction plote
  
  mpse <- dataset %>% as.MPSE()
  mpse %<>% mp_rrarefy()
  mpse %<>% mp_cal_rarecurve(.abundance=RareAbundance, chunks=10, action="add")
  mpse %<>% mp_cal_alpha(.abundance=RareAbundance)
  
  
  ## Generate tables:
  # - titleStyle : style object to use for title
  
  xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
    rows <-createRow(sheet,rowIndex=rowIndex)
    sheetTitle <-createCell(rows, colIndex=1)
    setCellValue(sheetTitle[[1,1]], title)
    setCellStyle(sheetTitle[[1,1]], titleStyle)
  }
  
  
  plot_types <- c("Bar_plot", "Alpha_rarefaction_curve", "Hierarchical_cluster", "Alpha_Diversity", "PCOA1")
  for (plot in plot_types){
    dir.create(file.path(projectName,"figures",plot),showWarnings = FALSE, recursive = TRUE)
  }
  
  p1 <- mpse %>% 
    mp_plot_rarecurve(
      .rare = RareAbundanceRarecurve, 
      .alpha = Observe,
      .group = Group
    ) +
    theme(axis.text=element_text(angle = 90, vjust = 0.5, hjust=1, size=8), panel.grid=element_blank(),
          strip.background = element_rect(colour=NA,fill="grey"),
          strip.text.x = element_text(face="bold"))
  
  rareplot_group <- p1
  
  ggsave(file=paste0(projectName,"/figures/","Alpha_rarefaction_curve/","rarefaction_curve_group.svg"), rareplot_group, width = 20, dpi = 300, units = "cm", device='svg')
  ggsave(file=paste0(projectName,"/figures/","Alpha_rarefaction_curve/","rarefaction_curve_group.png"), rareplot_group, width = 20, dpi = 300, units = "cm", device='png')
  
  assign("mpse", value = mpse, envir = .GlobalEnv)
  assign("dataset", value = dataset, envir = .GlobalEnv)
  
  #### Group-wise Analysis and Comparisons ####
  GROUPS <- unique(SAMPLE_TABLE$Group)
  GP <- unique(GROUPS)
  print (GP)
  
  for (i in seq_along(GP[-length(GP)]))
  {
    GP2 <- GP[-(1:i)]
    print (GP2)
    for (j in seq_along(GP2))
    {
      if (GP[i] == GP2[j])
      {
        next
      }
      else {
        
        wb<-createWorkbook(type="xlsx")
        wb2<-createWorkbook(type="xlsx")
        # Title and sub title styles
        TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16,
                                           color="dodgerblue4", isBold=TRUE, underline=1)
        SUB_TITLE_STYLE <- CellStyle(wb) +
          Font(wb,  heightInPoints=14,
               isItalic=TRUE, isBold=FALSE)
        # Styles for the data table row/column names
        TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
        TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
          Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
          Border(color="black", position=c("TOP", "BOTTOM"),
                 pen=c("BORDER_THIN", "BORDER_THICK"))
        
        
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
        
        print (GP[i])
        print (GP2[j])
        print (i)
        print (j)
        assign ("GP", value = GP, envir = .GlobalEnv)
        assign ("i", value = i, envir = .GlobalEnv)
        assign ("GP2", value = GP2, envir = .GlobalEnv)
        assign ("j", value = j, envir = .GlobalEnv)
        ph1 <- subset_samples(physeq = dataset, Group %in% c(GP[i], GP2[j]))
        print (ph1)
        assign("ph1", value = ph1, envir = .GlobalEnv)
        mpse <- ph1 %>% as.MPSE()
        
        Comp <- paste0(GP[i], "_vs_", GP2[j])
        dir.create(file.path(projectName, "GROUPWISE",Comp))
        dir.create(file.path(projectName,"GROUPWISE", Comp,"tables"),  showWarnings = FALSE,recursive = TRUE)
        plot_types <- c("Bar_plot", "Alpha_rarefaction_curve", "Hierarchical_cluster", "Alpha_Diversity", "PCOA1", "Box_plot", "Heatmap")
        for (plot in plot_types){
          dir.create(file.path(projectName,"GROUPWISE",Comp,"figures",plot),showWarnings = FALSE, recursive = TRUE)
        }
        mpse %<>% mp_rrarefy()
        mpse %<>% mp_cal_rarecurve(.abundance=RareAbundance, chunks=1000, action="add")
        mpse %<>% mp_cal_alpha(.abundance=RareAbundance)
        
        ## Plot Rarefaction Curve ###
        
        p2 <- mpse %>% 
          mp_plot_rarecurve(
            .rare = RareAbundanceRarecurve, 
            .alpha = "Observe", 
            .group = Group, 
            plot.group = TRUE
          ) +
          scale_color_manual(values=c("orangered", "dodgerblue")) +
          scale_fill_manual(values=c("orangered", "dodgerblue"),guide="none") + 
          theme(axis.text=element_text(angle = 90, vjust = 0.5, hjust=1, size=8), panel.grid=element_blank(),
                strip.background = element_rect(colour=NA,fill="grey"),
                strip.text.x = element_text(face="bold"))
        
        rareplot_group <- p2
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Alpha_rarefaction_curve/","rarefaction_curve_group.svg"), rareplot_group, width = 20, dpi = 300, units = "cm", device='svg')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Alpha_rarefaction_curve/","rarefaction_curve_group.png"), rareplot_group, width = 20, dpi = 300, units = "cm", device='png')
        
        
        ## Plot Alpha Diversity Metrics ##
        
        f1 <- mpse %<>% 
          mp_cal_alpha(.abundance=RareAbundance)
        f1
        
        ## Plot the Alpha Diversity groupwise
        f2 <- f1 %>% 
          mp_plot_alpha(
            .group=Group, 
            .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
          ) +
          scale_fill_manual(values=c("orangered", "dodgerblue"), guide="none") +
          scale_color_manual(values=c("orangered", "dodgerblue"), guide="none") +
          theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), strip.text.x.top = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"))
        
        # Plot the Alpha Diversity Individual Plots
        
        f3 <- mp_plot_alpha(.data = mpse, .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)) +
          ylab(label = "ALPHA INDEX VALUES") + theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), strip.text.x.top = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"))
        
        alphadiv <- f2
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Alpha_Diversity/","alpha_diversity.svg"), alphadiv, width = 20, height = 15, dpi = 300, units = "cm", device='svg')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Alpha_Diversity/","alpha_diversity.png"), alphadiv, width = 20, height = 15, dpi = 300, units = "cm", device='png')
        
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Alpha_Diversity/","Individual_alpha_diversity.svg"), f3, width = 20, height = 15, dpi = 300, units = "cm", device='svg')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Alpha_Diversity/","Individual_alpha_diversity.png"), f3, width = 20, height = 15, dpi = 300, units = "cm", device='png')
        
        ## Calculate Alpha Diversity Stats ##
        
        ## Generate the Table sof Alpha Diversity and Beta Diversity
        ############################## TABLE ##############################
        wb3 <- createWorkbook(type = "xlsx")
        TITLE_STYLE3 <- CellStyle(wb3)+ Font(wb3,  heightInPoints=16,
                                           color="dodgerblue4", isBold=TRUE, underline=1)
        SUB_TITLE_STYLE3 <- CellStyle(wb3) +
          Font(wb3,  heightInPoints=14,
               isItalic=TRUE, isBold=FALSE)
        # Styles for the data table row/column names
        TABLE_ROWNAMES_STYLE3 <- CellStyle(wb3) + Font(wb3, isBold=TRUE)
        TABLE_COLNAMES_STYLE3 <- CellStyle(wb3) + Font(wb3, isBold=TRUE) +
          Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
          Border(color="black", position=c("TOP", "BOTTOM"),
                 pen=c("BORDER_THIN", "BORDER_THICK"))
        #######################################################################
        
        ALPHA_DIV_TABLE <-  data.frame(mpse@colData@listData$sample.name, mpse@colData@listData$Group, mpse@colData@listData$Observe, mpse@colData@listData$Chao1, mpse@colData@listData$ACE, mpse@colData@listData$Shannon, mpse@colData@listData$Simpson, mpse@colData@listData$Pielou)
        colnames(ALPHA_DIV_TABLE) <- c("sample.name", "Group", "Observe", "Chao1", "ACE", "Shannon", "Simpson", "Pielou")
        sheet <- createSheet(wb3, sheetName = paste0("ALPHA_DIVERSITY_TABLE"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("ALPHA DIVERSITY METRICS TABLE"), titleStyle = TITLE_STYLE3)
        addDataFrame(ALPHA_DIV_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE3, rownamesStyle = TABLE_ROWNAMES_STYLE3)
        
        ######### Calculate the Statsitics of Alpha Diversity ############
        
        ALPHA_STATS_TABLE <- data.frame()
        ALPHA_VAL <- c("Observe", "Chao1", "ACE", "Shannon", "Simpson", "Pielou")
        for (al in seq_along(ALPHA_VAL))
        {
          MEAN_GP1 <- mean(ALPHA_DIV_TABLE[ALPHA_DIV_TABLE$Group == GP[i], ALPHA_VAL[al]])
          MEAN_GP2 <- mean(ALPHA_DIV_TABLE[ALPHA_DIV_TABLE$Group == GP2[j], ALPHA_VAL[al]])
          SD_GP1 <- sd(ALPHA_DIV_TABLE[ALPHA_DIV_TABLE$Group == GP[i], ALPHA_VAL[al]])
          SD_GP2 <- sd(ALPHA_DIV_TABLE[ALPHA_DIV_TABLE$Group == GP2[j], ALPHA_VAL[al]])
          SEM_GP1 <- SD_GP1/sqrt(length(ALPHA_DIV_TABLE[ALPHA_DIV_TABLE$Group == GP[i], ALPHA_VAL[al]]))
          SEM_GP2 <- SD_GP2/sqrt(length(ALPHA_DIV_TABLE[ALPHA_DIV_TABLE$Group == GP2[j], ALPHA_VAL[al]]))
          TAB2 <- ALPHA_DIV_TABLE[ALPHA_DIV_TABLE$Group %in% c(GP[i], GP2[j]),]
          
          TEST <- pairwise.wilcox.test(x = TAB2[,ALPHA_VAL[al]], g = TAB2[,"Group"], p.adjust.method = "fdr", paired = FALSE)
          P_VAL <- TEST$p.value
          T_TEST <- pairwise.t.test(x = TAB2[,ALPHA_VAL[al]], g = TAB2[,"Group"], p.adjust.method = "fdr", paired = FALSE)
          TP_VAL <- T_TEST$p.value
          RES <- c(paste0(ALPHA_VAL[al]), MEAN_GP1, MEAN_GP2, SD_GP1, SD_GP2, SEM_GP1, SEM_GP2, TP_VAL, P_VAL)
          ALPHA_STATS_TABLE <- rbind(ALPHA_STATS_TABLE, RES)
        }
        colnames(ALPHA_STATS_TABLE) <- c("Alpha Metric", paste0(GP[i], " [Mean]"), paste0(GP2[j], " [Mean]"), paste0(GP[i], " [Standard Deviation]"), paste0(GP[j], " [Standard Deviation]"), paste0(GP[i], " [SEM]"), paste0(GP[j], " [SEM]"), "FDR corrected p-value [Student's t-test]" ,"FDR corrected p-value [Wilcoxon's rank sum test]")
        sheet <- createSheet(wb3, sheetName = paste0("ALPHA_DIVERSITY_STATISTICS"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("ALPHA DIVERSITY STATISTICAL ANALYSIS"), titleStyle = TITLE_STYLE3)
        addDataFrame(ALPHA_STATS_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE3, rownamesStyle = TABLE_ROWNAMES_STYLE3)
        
        
        #### Beta_Diversity Plots #####
        
        B1 <- mpse %<>% mp_cal_dist(.abundance=Abundance, distmethod="bray")
        
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
        sheet <- createSheet(wb3, sheetName = paste0("BETA_DIVERSITY_MATRIX"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("BETA DIVERSITY DISTANCES ACROSS SAMPLES"), titleStyle = TITLE_STYLE3)
        addDataFrame(BRAY_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE3, rownamesStyle = TABLE_ROWNAMES_STYLE3)
        
        # Calculate Stats of Beta Diversity ##
        BETWEEN_SAMPLE_DISTANCE <- BETWEEN_SAMPLE_DISTANCE[!is.na(BETWEEN_SAMPLE_DISTANCE)]
        WITHIN_SAMPLE_DISTANCE <- WITHIN_SAMPLE_DISTANCE[!is.na(WITHIN_SAMPLE_DISTANCE)]
        Distances <- c(BETWEEN_SAMPLE_DISTANCE, WITHIN_SAMPLE_DISTANCE)
        Type_Distance <- c(rep("BETWEEN_SAMPLE_DISTANCE", length(BETWEEN_SAMPLE_DISTANCE)), rep("WITHIN_SAMPLE_DISTANCE", length(WITHIN_SAMPLE_DISTANCE)))
        BRAY_STATS <- data.frame(Distances, Type_Distance)
        MEAN_BTWN <- mean(BRAY_STATS[BRAY_STATS$Type_Distance == "BETWEEN_SAMPLE_DISTANCE", "Distances"], na.rm = TRUE)
        MEAN_WITHN <- mean(BRAY_STATS[BRAY_STATS$Type_Distance == "WITHIN_SAMPLE_DISTANCE", "Distances"], na.rm = TRUE)
        TEST <- pairwise.wilcox.test(x = BRAY_STATS$Distances, g = BRAY_STATS$Type_Distance, p.adjust.method = "fdr", paired = FALSE)
        BETA_DIVERSITY_STATS <- data.frame(MEAN_BTWN, MEAN_WITHN, TEST$p.value)
        colnames(BETA_DIVERSITY_STATS) <- c("Between Group Samples Distances (Mean)", "Within Group Samples Distance (Mean)", "Wilcoxon's test (P-value)")
        sheet <- createSheet(wb3, sheetName = paste0("BETA_DIVERSITY_STATS"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("BETA DIVERSITY DISTANCES STATISTICAL ANALYSIS"), titleStyle = TITLE_STYLE3)
        addDataFrame(BETA_DIVERSITY_STATS, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE3, rownamesStyle = TABLE_ROWNAMES_STYLE3)
        
        ## Plot the Beta Diversity ANalysis ##
        
        betadiv <- B1 %>% mp_plot_dist(.distmethod = bray, .group = Group, group.test=TRUE, textsize=2)
        
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","PCOA1/","beta_diversity.svg"), betadiv, width = 20, dpi = 300, units = "cm", device='svg')
        ggsave(file=paste0(projectName,"/GROUPWISE/", Comp,"/figures/","PCOA1/","beta_diversity.png"), betadiv, width = 20, dpi = 300, units = "cm", device='png')
        
        #PCOA
        
        P <- mpse %<>% 
          mp_cal_pcoa(.abundance=Abundance, distmethod="bray")
        
        p1 <- P %>%
          mp_plot_ord(
            .ord = pcoa, 
            .group = Group, 
            .color = Group, 
            .size = 1.2,
            .alpha = 1,
            ellipse = TRUE,
            show.legend = FALSE # don't display the legend of stat_ellipse
          ) +
          scale_fill_manual(values=c("#00A087FF", "#3C5488FF")) +
          scale_color_manual(values=c("#00A087FF", "#3C5488FF")) 
        
        # The size of point also can be mapped to other variables such as Observe, or Shannon 
        # Then the alpha diversity and beta diversity will be displayed simultaneously.
        
        p2 <- mpse %>% 
          mp_plot_ord(
            .ord = pcoa, 
            .group = Group, 
            .color = Group, 
            .size = Observe, 
            .alpha = Shannon,
            ellipse = TRUE,
            show.legend = FALSE # don't display the legend of stat_ellipse 
          ) +
          scale_fill_manual(
            values = c("#00A087FF", "#3C5488FF"), 
            guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
          ) +
          scale_color_manual(
            values=c("#00A087FF", "#3C5488FF"),
            guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
          ) +
          scale_size_continuous(
            range=c(0.5, 3),
            guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
          )
        PCOA <- p1 + p2
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","PCOA1/","PCOA.png"), PCOA, width = 30, height = 15, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","PCOA1/","PCOA.svg"), PCOA, width = 30, height = 15, dpi = 300, units = "cm", device='svg')
        
        
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
          geom_tippoint(aes(color=Group)) +
          geom_tiplab(aes(label = sample.name), as_ylab = TRUE) +
          ggplot2::scale_x_continuous(expand=c(0, 0.01))
        
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Hierarchical_cluster/","hierarchical_cluster_samples.png"), p, width = 30, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Hierarchical_cluster/","hierarchical_cluster_samples.svg"), p, width = 30, dpi = 300, units = "cm", device='svg')
        
        #### Barplots #####
        
        RareAbund <- mpse %<>%
          mp_cal_abundance( 
            .abundance = RareAbundance
          ) %>%
          mp_cal_abundance( 
            .abundance=RareAbundance,
            .group=Group
          )
        
        ### Phyla ###
        
        phyla.tb <- RareAbund %>% 
          mp_extract_abundance(taxa.class=Phylum, topn=50)
        phyla.tb$label <- gsub(x = phyla.tb$label, pattern = ".*([a-z]__)", "\\1")
        # The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
        bar <- phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phylum="label")
        bar
        p1 <- p + 
          geom_fruit(
            data=bar,
            geom=geom_col,
            mapping = aes(x = RelRareAbundanceBySample, 
                          y = Sample, 
                          fill = Phylum
            ),
            orientation = "y",
            #offset = 0.4,
            pwidth = 3, 
            axis.params = list(axis = "x", 
                               title = "The relative abundance of Phyla (%)",
                               title.size = 4,
                               text.size = 2,
                               limits = c(0,120),
                               vjust = 1),
            grid.params = list()
          )
        p1 <- p1 + theme(legend.position = "bottom")
        ggsave(file=paste0(projectName,"/GROUPWISE/", Comp,"/figures/","Bar_plot/","Phylum_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/", Comp,"/figures/","Bar_plot/","Phylum_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        print (bar)
        BAR <- bar[c("Phylum","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Phylum_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Phylum Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        ### Class ###
        
        class.tb <- RareAbund %>% 
          mp_extract_abundance(taxa.class=Class, topn=50)
        class.tb$label <- gsub(x = class.tb$label, pattern = "^([a-z]__).*__(\\S+)", "\\1\\2")
        # The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
        bar <- class.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Class="label")
        p1 <- p + 
          geom_fruit(
            data=bar,
            geom=geom_col,
            mapping = aes(x = RelRareAbundanceBySample, 
                          y = Sample, 
                          fill = Class
            ),
            orientation = "y",
            #offset = 0.4,
            pwidth = 3, 
            axis.params = list(axis = "x", 
                               title = "The relative abundance of Class (%)",
                               title.size = 4,
                               text.size = 2,
                               limits = c(0,120),
                               vjust = 1),
            grid.params = list()
          )
        p1 <- p1 + theme(legend.position = "bottom")
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Class_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Class_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        BAR <- bar[c("Class","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Class_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Class Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        ## Order ##
        
        order.tb <- RareAbund %>% 
          mp_extract_abundance(taxa.class=Order, topn=50)
        order.tb$label <- gsub(x = order.tb$label, pattern = "^([a-z]__).*__(\\S+)", "\\1\\2")
        # The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
        bar <- order.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Order="label")
        p1 <- p + 
          geom_fruit(
            data=bar,
            geom=geom_col,
            mapping = aes(x = RelRareAbundanceBySample, 
                          y = Sample, 
                          fill = Order
            ),
            orientation = "y",
            #offset = 0.4,
            pwidth = 3, 
            axis.params = list(axis = "x", 
                               title = "The relative abundance of Order (%)",
                               title.size = 4,
                               text.size = 2,
                               limits = c(0,120),
                               vjust = 1),
            grid.params = list()
          )
        p1 <- p1 + theme(legend.position = "bottom")
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Order_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Order_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        BAR <- bar[c("Order","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Order_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Order Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        #write.table(BAR, file=paste0(projectName, "/GROUPWISE/", Comp,"/tables/","Order_abundance_table.xlsx"), sep = "\t", row.names = FALSE, quote = FALSE)
        
        ### Family ###
        
        family.tb <- RareAbund %>% 
          mp_extract_abundance(taxa.class=Family, topn=50)
        family.tb$label <- gsub(x = family.tb$label, pattern = "^([a-z]__).*__(\\S+)", "\\1\\2")
        # The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
        bar <- family.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Family="label")
        p1 <- p + 
          geom_fruit(
            data=bar,
            geom=geom_col,
            mapping = aes(x = RelRareAbundanceBySample, 
                          y = Sample, 
                          fill = Family
            ),
            orientation = "y",
            #offset = 0.4,
            pwidth = 3, 
            axis.params = list(axis = "x", 
                               title = "The relative abundance of Family (%)",
                               title.size = 4,
                               text.size = 2,
                               limits = c(0,120),
                               vjust = 1),
            grid.params = list()
          )
        
        p1 <- p1 + theme(legend.position = "bottom")
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Family_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Family_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        BAR <- bar[c("Family","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Family_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Family Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)

        
        ### Genus ###
        
        genus.tb <- RareAbund %>% 
          mp_extract_abundance(taxa.class=Genus, topn=50)
        genus.tb$label <- gsub(x = genus.tb$label, pattern = "^([a-z]__).*__(\\S+)", "\\1\\2")
        # The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
        bar <- genus.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Genus="label")
        p1 <- p + 
          geom_fruit(
            data=bar,
            geom=geom_col,
            mapping = aes(x = RelRareAbundanceBySample, 
                          y = Sample, 
                          fill = Genus
            ),
            orientation = "y",
            #offset = 0.4,
            pwidth = 3, 
            axis.params = list(axis = "x", 
                               title = "The relative abundance of Genus (%)",
                               title.size = 4,
                               text.size = 2,
                               limits = c(0,120),
                               vjust = 1),
            grid.params = list()
          )
        p1 <- p1 + theme(legend.position = "bottom")
        
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Genus_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Genus_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        BAR <- bar[c("Genus","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Genus_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Genus Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        #write.table(BAR, file=paste0(projectName, "/GROUPWISE/", Comp,"/tables/","Genus_abundance_table.xlsx"), sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Species:
        
        species.tb <- RareAbund %>% 
          mp_extract_abundance(taxa.class=Species, topn=50)
        species.tb$label <- gsub(x = species.tb$label, pattern = ".*([a-z]__)", "\\1")
        # The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
        bar <- species.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Species="label")
        p1 <- p + 
          geom_fruit(
            data=bar,
            geom=geom_col,
            mapping = aes(x = RelRareAbundanceBySample, 
                          y = Sample, 
                          fill = Species
            ),
            orientation = "y",
            #offset = 0.4,
            pwidth = 3, 
            axis.params = list(axis = "x", 
                               title = "The relative abundance of Species (%)",
                               title.size = 4,
                               text.size = 2,
                               limits = c(0,120),
                               vjust = 1),
            grid.params = list()
          )
        p1 <- p1 + theme(legend.position = "bottom")
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Species_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Species_barplots.svg"), p1, width = 50, height = 35,dpi = 300, units = "cm", device='svg')
        BAR <- bar[c("Species","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Species_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Species Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)

        # Cladogram <- plot_cladogram(mm_aldex, color = c("darkgreen","red")) + theme(plot.margin = margin(0, 0, 0, 0))
        # ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Cladogram/CLADOGRAM_PLOT.png"), Cladogram, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        # ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Cladogram/CLADOGRAM_PLOT.svg"), Cladogram, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        
        ## Statistical Analysis ##
        
        TAXA <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
        for (t in seq_along(TAXA))
        {
          ph1_norm <- norm_tss(ph1) # Normalize with TSS
          mm_aldex <- run_aldex(ps = ph1_norm, group = "Group", pvalue_cutoff = 1, p_adjust = "fdr", method = "wilcox.test", norm_para = "tss", taxa_rank = TAXA[t])
          TABLE <- mm_aldex@marker_table@.Data
          names(TABLE) <- mm_aldex@marker_table@names
          sheet <- createSheet(wb2, sheetName = paste0(TAXA[t], "Statistical_Table"))
          xlsx.addTitle(sheet, rowIndex=1, title=paste0(TAXA[t]," Statistical Table"), titleStyle = TITLE_STYLE2)
          addDataFrame(TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
          #write.table(TABLE, file = paste0(projectName,"/GROUPWISE/", Comp,"/tables/", TAXA[t], "_STATISTICAL_TABLE.xlsx"), sep = "\t", quote = FALSE)
          
          # Plot the Box plots ##
	        mm_aldex2 <- mm_aldex
          otu_table(mm_aldex2) <- (otu_table(mm_aldex2)*10^6)
          
          if (TAXA[t] == "Species")
          {
            markers <- marker_table(mm_aldex2)[!str_detect(string = tax_table(mm_aldex2), pattern = "_s"),"feature"]
            p_abd <- plot_abundance(mm_aldex2, group = "Group", markers = markers$feature[1:25]) + scale_fill_brewer(palette = "Set1") + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='svg')
          }
          else if (TAXA[t] == "Genus")
          {
            markers <- marker_table(mm_aldex2)[!str_detect(string = tax_table(mm_aldex2), pattern = "_g__"),"feature"]
            p_abd <- plot_abundance(mm_aldex2, group = "Group", markers = markers$feature[1:25]) + scale_fill_brewer(palette = "Set1") + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='svg')
          }
          else if (TAXA[t] == "Family")
          {
            markers <- marker_table(mm_aldex2)[!str_detect(string = tax_table(mm_aldex2), pattern = "_f__"),"feature"]
            p_abd <- plot_abundance(mm_aldex2, group = "Group", markers = markers$feature[1:25]) + scale_fill_brewer(palette = "Set1") + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='svg')
          }
          else {
            p_abd <- plot_abundance(mm_aldex2, group = "Group") + scale_fill_brewer(palette = "Set1") + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 35, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
          }
          print ("Plotted the Box plots")
          
          # Plot Heatmap:
          HEATMAP <- plot_heatmap(mm_aldex,  group = "Group")
          print ("Got the heatmap")
          png(file = paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Heatmap/",TAXA[t], "_HEATMAP_PLOT.png"), height = 10, width = 16, units = 'in', res = 300)
          plot(HEATMAP)
          dev.off()
          print ("Plotted the Heatmap")
        }
        
        print ("Generate the Excel Tables")
        saveWorkbook(wb, paste0(projectName, "/GROUPWISE/", Comp, "/tables/TAXA_ABUNDANCE_TABLE.xlsx"))
        saveWorkbook(wb2, paste0(projectName, "/GROUPWISE/", Comp, "/tables/STATISTICAL_ANALYSIS_TABLE.xlsx"))
        saveWorkbook(wb3, paste0(projectName, "/GROUPWISE/", Comp, "/tables/DIVERSITY_ANALYSIS_TABLE.xlsx"))
      }
    }
  }
}

Plot_16S_rRNA(biom = biom, tree = tree, metadata = metadata, projectName = projectName)
