

#-----------------------------------------------------------------------
# Summary: Histogram plots with Gene names
#-----------------------------------------------------------------------
# Script to make histogram plots for 16p-22q FC Correlation Per Gene with 
# gene names shown for 16p11.2 (or 22q11.2) region genes


#-----------------------------------------------------------------------
# set wd
#-----------------------------------------------------------------------

# SET the Working DIRECTORY
in_wd <- 'D:/SSDN_KD/SSDN_Postdoc/Temp_Analysis/AHBA_analysis_MIST_apr2020/GeneExp_and_CNV_FCsignatures/code'
setwd(in_wd)

# specify the data and plots directories
data_dir <- '../data'
plots_dir <- '../plots'

#-----------------------------------------------------------------------
# libraries
#-----------------------------------------------------------------------

library(readxl)
library(ggplot2)
library(ggrepel)
library(tidyverse)

#-----------------------------------------------------------------------
# functions
#-----------------------------------------------------------------------


fHistPlot_save_png <- function(plot_file_name,in_df_hist,in_xlabel,in_ylabel,key.gns,in_nbins){
  
  # define plot width, height and resolution
  in_width <- 0.8*8
  in_height <- 0.8*6
  in_png_res <- 1000
  
  # Define xlim and ylim (to make all plots within same Correlation value range)
  xlim1 <- -0.8
  xlim2 <- 0.8
  ylim1 <- 0
  ylim2 <- 500
  
  # Fontsize for Xticks and Yticks
  in_fontsize <- 24
  
  # Make a new df: plot_df with column for key.gns (select gene names to be displayed)
  selected_obs <- key.gns
  values = in_df_hist$Corr
  id = in_df_hist$gene
  pval = in_df_hist$pval
  
  plot_df <- tibble(id = id,
                    values = values, pval = pval) %>%
    mutate(obs_labels = ifelse(id %in% selected_obs, id, NA),
           bins = as.factor( as.numeric( cut(values, in_nbins)))) # cutting into in_nbins
  
  # Make a Subset df: label_df for key.gns (select gene names to be displayed)
  label_df<- plot_df %>% filter(id %in% selected_obs) %>% left_join(plot_df, by = 'bins') %>% 
    group_by(values = values.x, obs_labels = obs_labels.x) %>% count
  
  # Add p-value column (To be used for scaling Gene Names (label) font-size)
  temp_subset = in_df_hist[in_df_hist$gene %in% selected_obs,]
  temp_subset = temp_subset[match(label_df$obs_labels, temp_subset$gene),]
  label_df$pval = temp_subset$pval
  
  
  # Make Histogram and save as png
  png(plot_file_name,units="in", width=in_width, height=in_height, res=in_png_res)
  print({ ggplot(plot_df, aes(values)) + 
      coord_cartesian(xlim=c(xlim1,xlim2),ylim=c(ylim1,ylim2))+
      geom_histogram(bins = in_nbins,linetype="dashed", fill = "lightblue", 
                     alpha = 0.5, color = "gray") +
      xlab(in_xlabel) +  ylab(in_ylabel) +
      # add select gene names (text using geom_text_repel) 
      geom_text_repel(data =label_df, aes(label = obs_labels, y = n, size= pval),   
                      fontface = "bold", min.segment.length = 0,  # Draw all lines  
                      colour = "royalblue", segment.color = "royalblue", 
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.30, "lines")) + 
      # Remove legend, and adjust the minimum label font
      scale_size(guide = 'none',range = c(2,6), limits = c(0,4)) +  
      # adjust yTick labels: display select values only
      scale_y_continuous(breaks=c(0,200,400)) +
      theme_light() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(axis.text=element_text(size=in_fontsize,colour="black"),axis.title=element_text(size=in_fontsize)) +
      # Add two vertical lines: 5th and 95th percentile of Correlation values
      geom_vline(xintercept = quantile(plot_df$values,0.95),lty=2) +
      geom_vline(xintercept = quantile(plot_df$values,0.05),lty=2)
  })
  
  dev.off()

}



#-----------------------------------------------------------------------
# 1. read data: All GeneList (15633 used in our analysis) + CNV region gene lists
#-----------------------------------------------------------------------

# All gene name list (used in our analysis)
gene_set_all <- read.csv2(paste0(data_dir,"/gene_set_all.csv"), header = TRUE, stringsAsFactors = FALSE , sep = ",")

# CNV region gene name list
gene_set16p <- read.csv2 (paste0(data_dir,"/cnv_genes_cnv16p112.csv"),header = FALSE, stringsAsFactors = FALSE , sep = ",")
gene_set22q <- read.csv2 (paste0(data_dir,"/cnv_genes_cnv22q112.csv"),header = FALSE, stringsAsFactors = FALSE ,  sep = ",")


#--------------------------------------------------------------------------------
# 2. Read Correlation and pval per Gene
#--------------------------------------------------------------------------------

# read Correlation Per Gene and Pvalues from .xlsx file
in_xsl_filename <- paste0(data_dir,"/tab_CorrPerGene_Pval_16p22q_MIST64.xlsx")
in_data_corr <- as.data.frame(read_excel(in_xsl_filename,sheet = "Corr"))
in_data_pval <- as.data.frame(read_excel(in_xsl_filename,sheet ="Pval"))

# set first column (/variable) as rownames and drop the column
rownames(in_data_corr) <- in_data_corr[, 1] ## set rownames
in_data_corr <- in_data_corr[, -1]          ## remove the first column
in_data_pval <- in_data_pval[, -1]          ## remove the first column


#--------------------------------------------------------------------------------
# 3. Define the n_bins and CorrPerGene and CNVgene lists
#--------------------------------------------------------------------------------

# Define the number of histogram bins
in_nbins <- 100

# define the CorrPerGene set and CNVgene sets
list_CorrPerGene_set <- c("FC16pdel","FC22qdel")
list_CNVgene_set <- c("Genes16p","Genes22q")

# xlabel and ylabel (fixed) to be used
list_xlabel <- c("Correlation with 16p11.2 FC profile","Correlation with 22q11.2 FC profile")
in_ylabel <- "Number of genes"

#--------------------------------------------------------------------------------
# 4. Make histograms
#--------------------------------------------------------------------------------
# loop over CorrPerGene sets: FC16pdel and FC22qdel
for(loop_c in c(1:length(list_CorrPerGene_set)))
{
  # make a df with Corr, gene names, and pvalues (-log10)
  in_df_hist <- as.data.frame(as.numeric(in_data_corr[,loop_c]))
  in_df_hist$gene <- gene_set_all[,1]
  in_df_hist$pval <- -1*log10(as.numeric(in_data_pval[,loop_c])) # p-values in -log10 scale
  names(in_df_hist) <- c("Corr","gene","pval")
  
  in_xlabel <- list_xlabel[loop_c]

  # loop over CNVgene sets: genes names to be highlighted in the histogram
  for( loop_g in c(1:length(list_CNVgene_set)))
  {
    
    if( loop_g == 1){
      # 16p11.2 gene names to be displayed
      key.gns <- gene_set16p[,1]
    }
    else if( loop_g == 2){
      # 22q11.2 gene names to be displayed
      key.gns <- gene_set22q[,1]
    }
  
    # png file name (saved in plots_dir)
    plot_file_name <- paste0(plots_dir,"/hist_CorrPerGene_",list_CorrPerGene_set[loop_c],"_v_",list_CNVgene_set[loop_g],".png")
    
    # call the function to make histogram + show gene names
    fHistPlot_save_png(plot_file_name,in_df_hist,in_xlabel,in_ylabel,key.gns,in_nbins)
    
  }
  
}


#--------------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------------


