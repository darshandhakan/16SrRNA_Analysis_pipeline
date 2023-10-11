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
suppressMessages(library(viridis))
suppressMessages(library(viridisLite))

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
  suppressMessages(library(microshades))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(viridis))
  suppressMessages(library(viridisLite))
  
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
  
  cat ("###########\n\nGenerated the All Sample Rarefaction Plot\n\n##########")
  
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
        
        ## Generate the Taxonomy Tables ##
        
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
        
        ########################################################
        
        cat ("########## ---------- \n\nSTART GROUPWISE ANALYSIS ----------##########")
        
        print (paste0("Group1 is ",GP[i]))
        print (paste0("Group2 is ",GP2[j]))
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
        plot_types <- c("Bar_plot", "Alpha_rarefaction_curve", "Hierarchical_cluster", "Alpha_Diversity", "Beta_Diversity", "Box_plot", "Heatmap")
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
        
        cat ("##########\n\nCompleted Rarefaction Analysis of these Groups\n\n##########")
        
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
        
        ## Generate the Tables of Alpha Diversity and Beta Diversity
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
          
          ## Perform pairwise tests:
          
          GP1_TABLE <- TAB2[TAB2$Group == GP[i], c("sample.name", ALPHA_VAL[al])]
          GP2_TABLE <- TAB2[TAB2$Group == GP2[j], c("sample.name", ALPHA_VAL[al])]
          assign("GP1_TABLE", value = GP1_TABLE, envir = .GlobalEnv)
          GP1_TABLE$sample.name <- gsub(x=GP1_TABLE$sample.name, pattern = "TL.*_", replacement = "")
          GP2_TABLE$sample.name <- gsub(x=GP2_TABLE$sample.name, pattern = "TL.*_", replacement = "")
          GP1_TABLE$sample.name <- gsub(x=GP1_TABLE$sample.name, pattern = "Control_", replacement = "")
          GP2_TABLE$sample.name <- gsub(x=GP2_TABLE$sample.name, pattern = "Control_", replacement = "")
          NAMES <- intersect(GP1_TABLE$sample.name, GP2_TABLE$sample.name)
          GP1_TABLE <- GP1_TABLE[GP1_TABLE$sample.name %in% NAMES,]
          GP2_TABLE <- GP2_TABLE[GP2_TABLE$sample.name %in% NAMES,]
          PAIRED_W_TEST <- wilcox.test(x = GP1_TABLE[,ALPHA_VAL[al]], y = GP2_TABLE[,ALPHA_VAL[al]], paired = TRUE, alternative = "two.sided")
          PAIRED_T_TEST <- t.test(x = GP1_TABLE[,ALPHA_VAL[al]], y = GP2_TABLE[,ALPHA_VAL[al]], paired = TRUE, alternative = "two.sided")
          
          
          RES <- c(paste0(ALPHA_VAL[al]), MEAN_GP1, MEAN_GP2, SD_GP1, SD_GP2, SEM_GP1, SEM_GP2, TP_VAL, P_VAL, PAIRED_T_TEST$p.value, PAIRED_W_TEST$p.value)
          ALPHA_STATS_TABLE <- rbind(ALPHA_STATS_TABLE, RES)
        }
        colnames(ALPHA_STATS_TABLE) <- c("Alpha Metric", paste0(GP[i], " [Mean]"), paste0(GP2[j], " [Mean]"), paste0(GP[i], " [Standard Deviation]"), paste0(GP2[j], " [Standard Deviation]"), paste0(GP[i], " [SEM]"), paste0(GP2[j], " [SEM]"), "FDR corrected p-value [Student's t-test]" ,"FDR corrected p-value [Wilcoxon's rank sum test]", "Paired t-test", "Paired Wilcoxon test")
        sheet <- createSheet(wb3, sheetName = paste0("ALPHA_DIVERSITY_STATISTICS"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("ALPHA DIVERSITY STATISTICAL ANALYSIS"), titleStyle = TITLE_STYLE3)
        addDataFrame(ALPHA_STATS_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE3, rownamesStyle = TABLE_ROWNAMES_STYLE3)
        
        cat("#########------------\n\nCOMPLETED ALPHA DIVERSITY ANALYSIS\n\n----------##########")
        
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
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("BETA DIVERSITY (BRAY CURTIS) DISTANCES ACROSS SAMPLES"), titleStyle = TITLE_STYLE3)
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
        
        betadiv <- B1 %>% mp_plot_dist(.distmethod = bray, .group = Group, group.test=TRUE, textsize=2) + theme(axis.text = element_text(size = 14, face = "bold"))
        
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Beta_Diversity/","beta_diversity.svg"), betadiv, width = 20, dpi = 300, units = "cm", device='svg')
        ggsave(file=paste0(projectName,"/GROUPWISE/", Comp,"/figures/","Beta_Diversity/","beta_diversity.png"), betadiv, width = 20, dpi = 300, units = "cm", device='png')
        
        #### PCOA Analysis ######
        
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
          scale_fill_manual(values=c("orangered", "dodgerblue")) +
          scale_color_manual(values=c("orangered", "#3C5488FF")) 
        
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
            values = c("orangered", "dodgerblue"), 
            guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
          ) +
          scale_color_manual(
            values=c("orangered", "dodgerblue"),
            guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
          ) +
          scale_size_continuous(
            range=c(0.5, 3),
            guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
          )
        PCOA <- p1 + p2
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Beta_Diversity/","PCOA.png"), PCOA, width = 30, height = 15, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Beta_Diversity/","PCOA.svg"), PCOA, width = 30, height = 15, dpi = 300, units = "cm", device='svg')
        
        cat ("##########----------\n\nCOMPLETED ALPHA AND BETA DIVERSITY ANALYSIS\n\n----------##########")
        
        ### The hierarchical cluster result of samples #####
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
        
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Hierarchical_cluster/","hierarchical_cluster_samples.png"), p, width = 30, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Hierarchical_cluster/","hierarchical_cluster_samples.svg"), p, width = 30, dpi = 300, units = "cm", device='svg')
        
        cat ("##########----------\n\nCompleted Clustering of the Samples\n\n----------##########")
        
        
        #### Barplots #####
        
        RareAbund <- mpse %<>%
          mp_cal_abundance( 
            .abundance = RareAbundance
          ) %>%
          mp_cal_abundance( 
            .abundance=RareAbundance,
            .group=Group
          )
        
        ################################### 
        ## Now plot the barplots of Taxa ##
        ###################################
        
        my_colors <- c(microshades_palette("micro_cvd_purple", 5, lightest = TRUE), microshades_palette("micro_cvd_blue", 5, lightest = TRUE), microshades_palette("micro_cvd_green", 5, lightest = TRUE), microshades_palette("micro_cvd_gray", 5, lightest = TRUE), microshades_palette("micro_cvd_orange", 5, lightest = TRUE))
        COL <- sample(x = my_colors, size = 25, replace = FALSE)
        
        ## Get the Phylum  
        phyla.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Phylum)
        bar <- phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phylum="label")
        BAR <- bar[c("Phylum","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Phylum_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Phylum Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        phyla.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Phylum, topn=20)
        phyla.tb$label <- gsub(x = phyla.tb$label, pattern = ".*([a-z]__)", "\\1")
        bar <- phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phylum="label")
        
        # Color the text using Group
        
        myPalette <- c("dodgerblue", "purple")
        names(myPalette) <- levels(factor(bar$Group))
        bar$Phylum <- gsub(x = bar$Phylum, pattern = "p__", replacement = "")
        p1 <- ggplot(bar, aes(x=paste0(Sample, " (", Group, ")"), y=RelRareAbundanceBySample, fill = Phylum)) + geom_bar(stat = "identity", position = "stack", width = 0.8) + coord_flip() + scale_fill_manual(values = COL) +  ylab(label = "Relative Abundance of Phylum") + xlab(label = "Samples (GROUP)") + 
              theme(axis.text.y = element_text(size = 8, face = "bold", color = myPalette[bar$Group]), axis.text.x = element_text(size = 12, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Phylum_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Phylum_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        
        cat ("##########\n\nPhylum Bar plots completed\n\n##########")
        
        ####### Class ########
        
        class.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Class)
        bar <- class.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Class="label")
        ## Add taxa Table
        BAR <- bar[c("Class","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Class_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Class Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        # Plot the Class Taxa.
        class.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Class, topn = 20)
        class.tb$label <- gsub(x = class.tb$label, pattern = ".*([a-z]__)", "\\1")
        bar <- class.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Class="label")
        myPalette <- c("dodgerblue", "purple")
        names(myPalette) <- levels(factor(bar$Group))
        bar$Class <- gsub(x = bar$Class, pattern = "p__", replacement = "")
        p1 <- ggplot(bar, aes(x=paste0(Sample, " (", Group, ")"), y=RelRareAbundanceBySample, fill = Class)) + geom_bar(stat = "identity", position = "stack", width = 0.8) + coord_flip() + scale_fill_manual(values = COL) +  ylab(label = "Relative Abundance of Class") + xlab(label = "Samples (GROUP)") + 
            theme(axis.text.y = element_text(size = 8, face = "bold", color = myPalette[bar$Group]), axis.text.x = element_text(size = 12, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Class_barplots.png"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Class_barplots.svg"), p1, width = 50, height = 35, dpi = 300, units = "cm", device='svg')
        
        cat ("##########\n\nClass Bar Plots Completed\n\n##########")
        
        ########## Order ############
        
        order.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Order)
        bar <- order.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Order="label")
        ## Add taxa Table
        BAR <- bar[c("Order","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Order_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Order Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        # Plot the Order Taxa.
        
        order.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Order, topn = 20)
        order.tb$label <- gsub(x = order.tb$label, pattern = ".*([a-z]__)", "\\1")
        bar <- order.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Order="label")
        myPalette <- c("dodgerblue", "purple")
        names(myPalette) <- levels(factor(bar$Group))
        bar$Order <- gsub(x = bar$Order, pattern = "p__", replacement = "")
        p1 <- ggplot(bar, aes(x=paste0(Sample, " (", Group, ")"), y=RelRareAbundanceBySample, fill = Order)) + geom_bar(stat = "identity", position = "stack", width = 0.8) + coord_flip() + scale_fill_manual(values = COL) +  ylab(label = "Relative Abundance of Order") + xlab(label = "Samples (GROUP)") + 
          theme(axis.text.y = element_text(size = 8, face = "bold", color = myPalette[bar$Group]), axis.text.x = element_text(size = 12, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
       
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Order_barplots.png"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Order_barplots.svg"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='svg')
        
        cat ("##########\n\nOrder Bar Plots Completed\n\n##########")
        
        ######### Family ############
        
        family.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Family)
        bar <- family.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Family="label")
        ## Add taxa Table
        BAR <- bar[c("Family","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Family_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Family Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        # Plot the Family Taxa.
        family.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Family)
        strings <- c("o__", "c__", "p__")
        family.tb <- family.tb[!str_detect(string = family.tb$label, pattern = paste(strings, collapse = '|')), ]
        family.tb$label <- gsub(x = family.tb$label, pattern = ".*([a-z]__)", "\\1")
        bar <- family.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Family="label")
        bar <- bar[rev(order(bar$RareAbundance)),]
        Family <- unique(bar$Family)[1:20]
        bar <- bar[bar$Family %in% Family,]
        myPalette <- c("dodgerblue", "purple")
        names(myPalette) <- levels(factor(bar$Group))
        bar$Family <- gsub(x = bar$Family, pattern = "p__", replacement = "")
        p1 <- ggplot(bar, aes(x=paste0(Sample, " (", Group, ")"), y=RelRareAbundanceBySample, fill = Family)) + geom_bar(stat = "identity", position = "stack", width = 0.8) + coord_flip() + scale_fill_manual(values = COL) +  ylab(label = "Relative Abundance of Family") + xlab(label = "Samples (GROUP)") + 
              theme(axis.text.y = element_text(size = 8, face = "bold", color = myPalette[bar$Group]), axis.text.x = element_text(size = 12, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Family_barplots.png"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Family_barplots.svg"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='svg')
        
        cat ("##########\n\nFamily Bar Plots completed\n\n##########")
        
        ############# Genus #################
        
        genus.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Genus)
        bar <- genus.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Genus="label")
        ## Add taxa Table
        BAR <- bar[c("Genus","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Genus_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Genus Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        # Plot the Genus Taxa.
        genus.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Genus)
        strings <- c("f__", "o__", "c__", "p__")
        genus.tb <- genus.tb[!str_detect(string = genus.tb$label, pattern = paste(strings, collapse = '|')), ]
        genus.tb$label <- gsub(x = genus.tb$label, pattern = ".*([a-z]__)", "\\1")
        bar <- genus.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Genus="label")
        bar <- bar[rev(order(bar$RareAbundance)),]
        Genus <- unique(bar$Genus)[1:20]
        bar <- bar[bar$Genus %in% Genus,]
        myPalette <- c("dodgerblue", "purple")
        names(myPalette) <- levels(factor(bar$Group))
        bar$Genus <- gsub(x = bar$Genus, pattern = "p__", replacement = "")
        p1 <- ggplot(bar, aes(x=paste0(Sample, " (", Group, ")"), y=RelRareAbundanceBySample, fill = Genus)) + geom_bar(stat = "identity", position = "stack", width = 0.8) + coord_flip() + scale_fill_manual(values = COL) +  ylab(label = "Relative Abundance of Genus") + xlab(label = "Samples (GROUP)") + 
          theme(axis.text.y = element_text(size = 8, face = "bold", color = myPalette[bar$Group]), axis.text.x = element_text(size = 12, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Genus_barplots.png"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Genus_barplots.svg"), p1, width = 60, height = 35, dpi = 300, units = "cm", device='svg')
        
        cat ("##########\n\nGenus Bar Plots completed\n\n##########")
        
        ############# Species #################
        
        species.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Species)
        bar <- species.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Species="label")
        ## Add taxa Table
        BAR <- bar[c("Species","Sample","RareAbundance")]
        BAR <- spread(BAR, key = "Sample", value = "RareAbundance")
        sheet <- createSheet(wb, sheetName = paste0("Species_Abundance_Table"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("Species Abundance Table"), titleStyle = TITLE_STYLE)
        addDataFrame(BAR, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        # Plot the Species Taxa #
        species.tb <- RareAbund %>% mp_extract_abundance(taxa.class=Species)
        strings <- c("f__", "o__", "c__", "p__", "g__")
        species.tb <- species.tb[!str_detect(string = species.tb$label, pattern = paste(strings, collapse = '|')), ]
        species.tb$label <- gsub(x = species.tb$label, pattern = ".*([a-z]__)", "\\1")
        bar <- species.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Species="label")
        bar <- bar[rev(order(bar$RareAbundance)),]
        Species <- unique(bar$Species)[1:20]
        bar <- bar[bar$Species %in% Species,]
        myPalette <- c("dodgerblue", "purple")
        names(myPalette) <- levels(factor(bar$Group))
        bar$Species <- gsub(x = bar$Species, pattern = "p__", replacement = "")
        p1 <- ggplot(bar, aes(x=paste0(Sample, " (", Group, ")"), y=RelRareAbundanceBySample, fill = Species)) + geom_bar(stat = "identity", position = "stack", width = 0.8) + coord_flip() + scale_fill_manual(values = COL) +  ylab(label = "Relative Abundance of Class") + xlab(label = "Samples (GROUP)") + 
          theme(axis.text.y = element_text(size = 8, face = "bold", color = myPalette[bar$Group]), axis.text.x = element_text(size = 12, face = "bold"),panel.background = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Species_barplots.png"), p1, width = 80, height = 50, dpi = 300, units = "cm", device='png')
        ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Bar_plot/","Species_barplots.svg"), p1, width = 80, height = 50, dpi = 300, units = "cm", device='svg')
        
        cat ("##########\n\nSpecies Bar Plots Completed\n\n##########")
        
        ## Statistical Analysis ##
        
        TAXA <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
        for (t in seq_along(TAXA))
        {
          ph1_norm <- norm_tss(ph1) # Normalize with TSS
          tax_table(ph1_norm) <- gsub(x=tax_table(ph1_norm), pattern = "^.__", replacement = "")
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
            markers <- markers[1:20]
            p_abd <- plot_abundance(mm_aldex2, group = "Group", markers = markers$feature) + scale_fill_viridis(discrete = TRUE, alpha = 0.5) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='svg')
          }
          else if (TAXA[t] == "Genus")
          {
            markers <- marker_table(mm_aldex2)[!str_detect(string = tax_table(mm_aldex2), pattern = "_g__"),"feature"]
            markers <- markers[1:20]
            p_abd <- plot_abundance(mm_aldex2, group = "Group", markers = markers$feature) + scale_fill_viridis(discrete = TRUE, alpha = 0.5) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='svg')
          }
          else if (TAXA[t] == "Family")
          {
            markers <- marker_table(mm_aldex2)[!str_detect(string = tax_table(mm_aldex2), pattern = "_f__"),"feature"]
            markers <- markers[1:20]
            p_abd <- plot_abundance(mm_aldex2, group = "Group", markers = markers$feature) + scale_fill_viridis(discrete = TRUE, alpha = 0.5) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='svg')
          }
          else if (TAXA[t] == "Order")
          {
            markers <- marker_table(mm_aldex2)[!str_detect(string = tax_table(mm_aldex2), pattern = "_o__"),"feature"]
            p_abd <- plot_abundance(mm_aldex2, group = "Group", markers = markers$feature) + scale_fill_viridis(discrete = TRUE, alpha = 0.5)  + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.png"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='png')
            ggsave(file=paste0(projectName,"/GROUPWISE/",Comp,"/figures/","Box_plot/",TAXA[t], "_BOX_PLOT.svg"), p_abd, width = 50, height = 60, dpi = 300, units = "cm", device='svg')
          }
          else {
            p_abd <- plot_abundance(mm_aldex2, group = "Group") + scale_fill_viridis(discrete = TRUE, alpha = 0.5)+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + xlab(label = "Relative Normalized Abundance * 10^6 in Log10 scale") + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
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
