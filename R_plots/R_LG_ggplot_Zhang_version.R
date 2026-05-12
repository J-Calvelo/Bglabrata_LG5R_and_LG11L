setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/New_Snail_Oregon/FST_data_iM_Zhang_2024")

library(ggplot2)
library(cowplot)
dir.create("Explore_Graphs", showWarnings = FALSE)
min_snp_count <- 5

data <- read.table("Good_Mean_Fsts_Windows_10000_Freqs_Zhang_iM_2024_masked_GWAS_RS.pileup.txt", sep = '\t', header = TRUE)
data <- data[data$Sites >= min_snp_count, ]

for (Scaffold in unique(data$Scaffold)) {
 out_PDF1 <- paste("Explore_Graphs/FST_",Scaffold,"_explore.pdf", sep = "" )

 title=paste("FST for Scaffold:",Scaffold)

 lg_data <- data[data$Scaffold == Scaffold,]
 Max_pos <- max(lg_data$Window) 
 grafica <- ggplot(lg_data, aes(x = Window, y = Fst) )  +
   theme_classic() +
   ggtitle(title) +
   scale_x_continuous(breaks = seq(1, Max_pos, by = 200000) ) +
   geom_point(size=0.5) +
   xlab("Window") + 
   ylab("Fst") +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 3))
 ggsave(out_PDF1, plot=grafica, width=25, height=8) 
}

data <- data[data$Sites >= min_snp_count,]

# LG10_Left plot figure
chr="CM074212.1"
start_region <-4648498
end_region <-6946461

additiona_start_region <- 15116011
additiona_end_region <-16477347


main_peak_start <- 5600000
main_peak_end <- 6200000

out_PDF <- "LG10_Left_Fst_Peak_Lines.pdf"
title <- "Average Fst window value: CM074212.1"

subset_data <- data[data$Scaffold == chr,]
Max_pos <- max(subset_data$Window) 

grafica <- ggplot(subset_data, aes(x = Window, y = Fst) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(0, Max_pos, by = 500000) ) +
  geom_vline(xintercept=start_region, linetype="dashed", color="blue2", size=0.5) +
  geom_vline(xintercept=main_peak_start, color="red2", linewidth=0.5) +
  geom_vline(xintercept=main_peak_end, color="red2", linewidth=0.5) +
  geom_vline(xintercept=end_region, linetype="dashed", color="blue2", size=0.5) +
  geom_point(size=0.2) +
  xlab("Window") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 3))
ggsave(out_PDF, plot=grafica, width=25, height=8) 


out_PDF <- "LG10_Left_Fst_Peak.pdf"
title <- "Average Fst window value: CM074212.1"
chr="CM074212.1"
subset_data <- data[data$Scaffold == chr,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"
start_region <- 4648498/1000000
end_region <- 6946461/1000000
subset_data <- cbind(subset_data,megabases)

additiona_start_region <- 15116011/1000000
additiona_end_region <-16477347/1000000

Max_pos <- max(subset_data$Window)/1000000 + 1
grafica <- ggplot(subset_data, aes(x = MB, y = Fst) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(0, Max_pos, by = 1) ) +
  geom_vline(xintercept=start_region, linetype="dashed", color="blue2", size=0.5) +
  geom_vline(xintercept=end_region, linetype="dashed", color="blue2", size=0.5) +
  
  geom_vline(xintercept=additiona_start_region, linetype="dashed", color="red2", size=0.5) +
  geom_vline(xintercept=additiona_end_region, linetype="dashed", color="red2", size=0.5) +
  
  geom_point(size=0.2) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=12, height=6) 


chr="CM074206.1"
start_region <- 43022826/1000000
end_region <- 47628293/1000000

title <- "Average Fst window value: CM074206.1"
out_PDF <- "LG5_Right_Fst_Peak.pdf"
subset_data <- data[data$Scaffold == chr,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"
subset_data <- cbind(subset_data,megabases)

Min_Fst_pos <- min(subset_data$Fst)
Max_pos <- max(subset_data$Window)/1000000 + 1
grafica <- ggplot(subset_data, aes(x = MB, y = Fst) )  +
  theme_classic() +
  ggtitle(title) +
  ylim(Min_Fst_pos, 0.08) +
  #scale_y_continuous(breaks = seq(0, 0.08, by = 0.02) ) +
  scale_x_continuous(breaks = seq(0, Max_pos, by = 1) ) +
  geom_vline(xintercept=start_region, linetype="dashed", color="red2", size=0.5) +
  geom_vline(xintercept=end_region, linetype="dashed", color="red2", size=0.5) +
  geom_point(size=0.2) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=12, height=6) 

############################################################################################

chr="CM074212.1"
start_region <- 4600000
end_region <- 7000000
out_PDF <- "LG10_Left_Fst_Peak_close_in.pdf"
subset_data <- data[data$Scaffold == chr,]
subset_data <- subset_data[subset_data$Window >= start_region,]
subset_data <- subset_data[subset_data$Window <= end_region,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"

start_region <- 4600000/1000000
end_region <- 7000000/1000000

subset_data <- cbind(subset_data,megabases)

title <- "Fst over Gene Models: LG10_Left "

grafica <- ggplot(subset_data, aes(x = MB, y = Fst) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(start_region, end_region, by = 0.5) ) +
  geom_point(size=1) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=12.1, height=6) 

############################################################################################

chr="CM074212.1"
start_region <- 4600000
end_region <- 7000000
out_PDF <- "LG10_Left_Het_Peak_close_in.pdf"
subset_data <- data[data$Scaffold == chr,]
subset_data <- subset_data[subset_data$Window >= start_region,]
subset_data <- subset_data[subset_data$Window <= end_region,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"

start_region <- 4600000/1000000
end_region <- 7000000/1000000

subset_data <- cbind(subset_data,megabases)

title <- "Het over Gene Models: LG10_Left "

grafica <- ggplot(subset_data, aes(x = MB, y = HetE) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(start_region, end_region, by = 0.5) ) +
  geom_point(size=1) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=12.1, height=6) 

############################################################################################
############################################################################################

# CM074206.1:43022826-47628293+

chr="CM074206.1"
start_region <- 43022826
end_region <- 47628293
out_PDF <- "LG5_Right_FST_Peak_close_in.pdf"
subset_data <- data[data$Scaffold == chr,]
subset_data <- subset_data[subset_data$Window >= start_region,]
subset_data <- subset_data[subset_data$Window <= end_region,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"

start_region <- 43000000/1000000
end_region <- 48000000/1000000


subset_data <- cbind(subset_data,megabases)
title <- "FST over Gene Models: LG5_Right "

grafica <- ggplot(subset_data, aes(x = MB, y = Fst) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(start_region, end_region, by = 0.5) ) +
  geom_point(size=1) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=12.1, height=6)   


############################################################################################
############################################################################################


############################################################################################

chr="CM074212.1"
start_region <- 15000000
end_region <- 16500000

out_PDF <- "LG10_Middle_Fst_Peak_close_in.pdf"
subset_data <- data[data$Scaffold == chr,]
subset_data <- subset_data[subset_data$Window >= start_region,]
subset_data <- subset_data[subset_data$Window <= end_region,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"

start_region <- start_region/1000000
end_region <- end_region/1000000

subset_data <- cbind(subset_data,megabases)

title <- "Fst over Gene Models: LG10_Middle "

grafica <- ggplot(subset_data, aes(x = MB, y = Fst) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(start_region, end_region, by = 0.1) ) +
  geom_point(size=1) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=18, height=6) 

out_PDF <- "LG10_Middle_HetE_Peak_close_in.pdf"
title <- "HetE over Gene Models: LG10_Middle "

grafica <- ggplot(subset_data, aes(x = MB, y = HetE) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(start_region, end_region, by = 0.1) ) +
  geom_point(size=1) +
  xlab("Position (MB)") + 
  ylab("HetE") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=18, height=6) 

chr="CM074212.1"
start_region <- 15000000
end_region <- 16500000

out_PDF <- "LG10_coverage_group_A_Fst_Peak_close_in.pdf"
subset_data <- data[data$Scaffold == chr,]
subset_data <- subset_data[subset_data$Window >= start_region,]
subset_data <- subset_data[subset_data$Window <= end_region,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"

start_region <- start_region/1000000
end_region <- end_region/1000000

subset_data <- cbind(subset_data,megabases)

title <- "MeanDepth: LG10_Middle "

grafica <- ggplot(subset_data, aes(x = MB, y = MeanPercentGroupA) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(start_region, end_region, by = 0.1) ) +
  geom_point(size=1) +
  xlab("Position (MB)") + 
  ylab("MeanPercentGroupA %") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=18, height=6) 



########################################################################################

	
chr="CM074212.1"
start_region <- 15100000
end_region <- 16500000

out_PDF <- "LG10_Middle_Fst_Peak_for_plot.pdf"
subset_data <- data[data$Scaffold == chr,]
subset_data <- subset_data[subset_data$Window >= start_region,]
subset_data <- subset_data[subset_data$Window <= end_region,]
megabases <- as.data.frame(subset_data$Window/1000000)
colnames(megabases) <- "MB"

start_region <- start_region/1000000
end_region <- end_region/1000000

subset_data <- cbind(subset_data,megabases)

title <- "Fst over Gene Models: LG10_Middle "

grafica <- ggplot(subset_data, aes(x = MB, y = Fst) )  +
  theme_classic() +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(start_region, end_region, by = 0.1) ) +
  geom_point(size=1) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  theme(axis.text.x = element_text(vjust = 1, size = 8))
ggsave(out_PDF, plot=grafica, width=18, height=6) 

###################################################################################
###################################################################################
###################################################################################

# Full Figures
data_Fst <- read.table("Good_Mean_Fsts_Windows_10000_Freqs_Zhang_iM_2024_masked_GWAS_RS.pileup.txt", sep = '\t', header = TRUE)
data_Fst <- data_Fst[data_Fst$Sites >= 150, ]
chromosome_gruide <- read.table("Chromosome_Guide.txt", sep = '\t', header = TRUE)

plot_dataset <- data.frame()
add_chr_nums <- data.frame()

postion_correction <- 0
chr_num=1
color_for_plot <- c() 

for (Scaffold in chromosome_gruide$Scaffold) {
  print(Scaffold)
  chromosome_number <- chromosome_gruide[chromosome_gruide$Scaffold == Scaffold, 2]
  pre_subdata <- data_Fst[data_Fst$Scaffold == Scaffold,]
  megabases <- as.data.frame(c(pre_subdata$Window)/1000000 + postion_correction)
  temp_num_note <-  c(as.character(chr_num), (postion_correction + (max(megabases) - min(megabases))/2) )
  add_chr_nums <- rbind(add_chr_nums, temp_num_note)
  colnames(add_chr_nums) <-c("chr", "Position") 

  postion_correction <- max(megabases)
  postion_correction + (max(megabases) - min(megabases))/2
  
  colnames(megabases) <- "MB"
  add_chr_group <- rep(chromosome_number, length(pre_subdata$Scaffold))
  
  pre_subdata <- cbind(pre_subdata, megabases, Chr = add_chr_group)
  plot_dataset <- rbind(plot_dataset,pre_subdata)
  
  if((chr_num %% 2) == 0) {
    color_for_plot <- c(color_for_plot, "red2")
  } else {
    color_for_plot <- c(color_for_plot, "blue2")
  }
  
  chr_num <- chr_num + 1
}

title <- "Fst Plot on iM_Zhong"
out_PDF <- "Genome_Wide_Fst_Plot.pdf"

grafica_fst <- ggplot(plot_dataset, aes(x = MB, y = Fst, color=factor(Chr)))  +
  theme_classic() +
  ggtitle(title ) +
  geom_point(size=0.1) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  scale_color_manual(values=color_for_plot) + 
  theme(axis.text.x = element_text(vjust = 1, size = 18), legend.position="none", axis.text.y = element_text(size = 18), plot.title = element_text(size = 24), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) 

counter = 1
for (i in 1:nrow(add_chr_nums)) {
  
  extract_info <- add_chr_nums[i,]
  chr_note <- extract_info$chr
  chr_x_pos <- as.numeric(extract_info$Position)
  
  grafica_fst <- grafica_fst + 
    annotate("text", x=chr_x_pos, y= -0.003 , label=chr_note, size = 4)
}

ggsave(out_PDF, plot=grafica_fst, width=12, height=6)

full_fst_plot <- plot_dataset

###########################################################################################################

data_Fst <- read.table("Good_Mean_Fsts_Windows_10000_Freqs_Zhang_iM_2024_masked_GWAS_RS.pileup.txt", sep = '\t', header = TRUE)
data_Fst <- data_Fst[data_Fst$Sites >= 5, ]

data_TMHMM_Full <- read.table("TMHMM_on_window_anot.txt", sep = '\t', header = TRUE)
data_TMHMM <- data_TMHMM_Full[data_TMHMM_Full$Annot_source == "Bg2024",]

plot_dataset <- data.frame()
postion_correction <- 0
chr_num=1
color_for_plot <- c() 

for (Scaffold in chromosome_gruide$Scaffold) {
  print(Scaffold)
  print(postion_correction)
  pre_subdata_fst <- data_Fst[data_Fst$Scaffold == Scaffold, ] 
  
  chromosome_number <- chromosome_gruide[chromosome_gruide$Scaffold == Scaffold, 2]
  pre_subdata <- data_TMHMM[data_TMHMM$Chr == Scaffold,]
  megabases_main <- as.data.frame(c(((pre_subdata$Start + ((pre_subdata$End - (pre_subdata$Start +1) )/2)))/1000000))
  megabases_main <- as.data.frame(megabases_main + postion_correction)

  postion_correction <- max(megabases_main)
  colnames(megabases_main) <- "MB"

  add_chr_group <- as.data.frame(rep(chromosome_number, length(pre_subdata$Chr)))
  colnames(add_chr_group) <- "Chr_num"
  
  pre_subdata <- cbind(pre_subdata, megabases_main, Chr = add_chr_group)
  plot_dataset <- rbind(plot_dataset,pre_subdata)
  
  if((chr_num %% 2) == 0) {
    color_for_plot <- c(color_for_plot, "red2")
  } else {
    color_for_plot <- c(color_for_plot, "blue2")
  }
  
  #############################################################
  ### Generate FST local peak
  # Prepare subplots: FST
  pre_subdata_fst <- data_Fst[data_Fst$Scaffold == Scaffold, ]
  megabases <- as.data.frame(c(pre_subdata_fst$Window)/1000000)
  colnames(megabases) <- "MB"
  pre_subdata_fst <- cbind(pre_subdata_fst,megabases)
  
  # Prepare subplots: TMHMM
  pre_subdata_TMHMM <- data_TMHMM_Full[data_TMHMM_Full$Chr == Scaffold, ]
  megabases <- as.data.frame(c(((pre_subdata_TMHMM$Start + ((pre_subdata_TMHMM$End - (pre_subdata_TMHMM$Start +1) )/2)))/1000000))
  colnames(megabases) <- "MB"
  pre_subdata_TMHMM <- cbind(pre_subdata_TMHMM, megabases)
  title <- paste("Nº TM1 Genes in iM_Zhong: ",Scaffold , sep="" )
  out_PDF <-paste("Fst_plus_TMHMM_",Scaffold,".pdf", sep="" )
  
  x_max <- max(c(pre_subdata_fst$MB, pre_subdata_TMHMM$MB ))
  
  grafica_fst_local <- ggplot(pre_subdata_fst, aes(x = MB, y = Fst))  +
    theme_classic() +
    ggtitle(title) +
    geom_point(size=0.1) +
    ylab("Fst") +
    xlab(NULL) +
    xlim(0, x_max) +
    scale_color_manual(values=color_for_plot) + 
    theme(axis.text.x=element_blank(), axis.ticks.x = element_line(color = "white"), axis.line.x =element_line(color = "white"), axis.text.y = element_text(size = 14), plot.title = element_text(size = 18), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))

  grafica_tmhmm_local <- ggplot(pre_subdata_TMHMM , aes(x = MB, y = N_Genes_Single_TMHMM, color=Annot_source, group=Annot_source))  +
    theme_classic() +
    geom_line(linewidth = 1) +
    xlab("Position (MB)") + 
    ylab("Nº TM1 genes") +
    xlim(0, x_max) +
    theme(axis.text.x = element_text(vjust = 1, size = 14), axis.text.y = element_text(size = 14), legend.position="bottom", plot.title = element_text(size = 18), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))

  print_final_plot <- plot_grid(grafica_fst_local, grafica_tmhmm_local, align= "v", ncol=1,rel_heights = c(1,2))  
  ggsave(out_PDF, plot=print_final_plot, width=12, height=6) 

  
  chr_num <- chr_num + 1
}

title <- "TMHMM genes in iM_Zhong genome (Annotation: Bg_2024) "
out_PDF <- "Genome_Wide_TMHMM_Plot.pdf"

grafica_TMHMM <- ggplot(plot_dataset, aes(x = MB, y = N_Genes_Single_TMHMM, color=factor(Chr_num)))  +
  theme_classic() +
  ggtitle(title) +
  geom_line(linewidth = 0.8) +
  xlab("Position (MB)") + 
  ylab("Nº TM1 genes") +
  scale_color_manual(values=color_for_plot) + 
  theme(axis.text.x = element_text(vjust = 1, size = 14), axis.text.y = element_text(size = 14), legend.position="none", plot.title = element_text(size = 18), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
ggsave(out_PDF, plot=grafica_TMHMM, width=12, height=6) 

print_final_plot <- plot_grid(grafica_fst, grafica_TMHMM, align= "v", ncol=1, rel_heights = c(1,1.5))
ggsave("Genome_Fst_plus_TMHMM_Anot_Bg_2024.pdf", plot=print_final_plot, width=12, height=8) 

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# Close up images of FST

title <- "(A) Fst across iM_Zhong genome"
grafica_fst <- ggplot(full_fst_plot, aes(x = MB, y = Fst, color=factor(Chr)))  +
  theme_classic() +
  ggtitle(title) +
  geom_point(size=0.1) +
  geom_hline(yintercept=0.04, linetype="dashed", color="cyan3", linewidth=0.5) +
  xlab("Position (MB)") + 
  ylab("Fst") +
  scale_color_manual(values=color_for_plot) + 
  theme(axis.text.x = element_text(vjust = 1, size = 18), plot.title=element_text(size = 22, face="bold") , legend.position="none", axis.text.y = element_text(size = 18), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))

counter = 1
for (i in 1:nrow(add_chr_nums)) {
  
  extract_info <- add_chr_nums[i,]
  chr_note <- extract_info$chr
  chr_x_pos <- as.numeric(extract_info$Position)
  
  grafica_fst <- grafica_fst + 
    annotate("text", x=chr_x_pos, y= -0.003 , label=chr_note, size = 4)
}


data_Fst <- read.table("Good_Mean_Fsts_Windows_10000_Freqs_Zhang_iM_2024_masked_GWAS_RS.pileup.txt", sep = '\t', header = TRUE)
data_Fst <- data_Fst[data_Fst$Sites >= 5, ]

Scaffold <- "CM074206.1"
pre_subdata_fst <- data_Fst[data_Fst$Scaffold == Scaffold, ] 
megabases <- as.data.frame(c(pre_subdata_fst$Window)/1000000)
colnames(megabases) <- "MB"
pre_subdata_fst <- cbind(pre_subdata_fst,megabases)

title <- "(B) Fst across Chr 5 "
grafica_fst_LG5 <- ggplot(pre_subdata_fst, aes(x = MB, y = Fst))  +
  theme_classic() +
  annotate("text",x=41700000/1000000, y= 0.042, label="BgSRR1", size = 6) +
  annotate("segment",x=39635000/1000000, y= 0.040, xend=43445001/1000000, yend=0.040, color="blue2", linewidth=0.8) +
  geom_vline(xintercept=43022826/1000000,linetype="dashed", color="blue2", linewidth=0.5) +
  geom_vline(xintercept=47628293/1000000,linetype="dashed", color="blue2", linewidth=0.5) +

  annotate("text",x=45485000/1000000, y= 0.052, label="LG5R", size = 6) +  
  #geom_vline(xintercept=44551840/1000000,linetype="dashed", color="red2", linewidth=0.5) +
  #geom_vline(xintercept=45787293/1000000,linetype="dashed", color="red2", linewidth=0.5) +
  ggtitle(title) +
  geom_point(size=0.1) +
  ylab("Fst") +
  xlab("Position (MB)") + 
  theme(axis.text.x = element_text(vjust = 1, size = 18), plot.title=element_text(size = 22, face="bold") , legend.position="none", axis.text.y = element_text(size = 18), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))

Scaffold <- "CM074212.1"
pre_subdata_fst <- data_Fst[data_Fst$Scaffold == Scaffold, ] 
megabases <- as.data.frame(c(pre_subdata_fst$Window)/1000000)
colnames(megabases) <- "MB"
pre_subdata_fst <- cbind(pre_subdata_fst,megabases)

title <- "(C) Fst across Chr 11"
grafica_fst_LG11 <- ggplot(pre_subdata_fst, aes(x = MB, y = Fst))  +
  theme_classic() +
  ggtitle(title) +
  geom_point(size=0.1) +
  ylab("Fst") +
  xlab("Position (MB)") + 
  annotate("text",x=5797479/1000000, y= 0.052, label="LG11L", size = 6 ) +  
  annotate("text",x=31500000/1000000, y= 0.045, label="RADres1", size = 6 ) +
  geom_vline(xintercept=4648498/1000000,linetype="dashed", color="blue2", linewidth=0.5) +
  geom_vline(xintercept=6946461/1000000,linetype="dashed", color="blue2", linewidth=0.5) +
  theme(axis.text.x = element_text(vjust = 1, size = 18), plot.title=element_text(size = 22, face="bold") , legend.position="none", axis.text.y = element_text(size = 18), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))

print_final_plot <- plot_grid(grafica_fst, grafica_fst_LG5, grafica_fst_LG11 , align= "v", ncol=1, rel_heights = c(2,1.5,1.5) )   
ggsave("Fst_Main_Figure_JC.pdf", plot=print_final_plot, width=12, height=16) 
ggsave("Fst_Main_Figure_JC.png", plot=print_final_plot, width=12, height=16) 


data_Fst_Zhang <- read.table("iM_Zhang_2024_data.txt", sep = '\t', header = TRUE)

Scaffold <- "BgChr_5"
pre_subdata_fst <- data_Fst_Zhang[data_Fst_Zhang$CHROM == Scaffold, ] 
megabases <- as.data.frame(c(pre_subdata_fst$BIN_START)/1000000)
colnames(megabases) <- "MB"
pre_subdata_fst <- cbind(pre_subdata_fst,megabases)

grafica_fst_zhang <- ggplot(pre_subdata_fst, aes(x = MB, y = MEAN_FST))  +
  theme_classic() +
  ggtitle("Chr 5: Zhang et al. 2024") +
  annotate("text",x=41700000/1000000, y= 0.55, label="BgSRR1", size = 6) +
  geom_vline(xintercept=39635000/1000000,linetype="dashed", color="blue2", linewidth=0.5) +
  geom_vline(xintercept=43445001/1000000,linetype="dashed", color="blue2", linewidth=0.5) +

  annotate("text",x=45787293/1000000, y= 0.5, label="LG5R", size = 6 ) +  
  geom_vline(xintercept=44551840/1000000,linetype="dashed", color="red2", linewidth=0.5) +
  geom_vline(xintercept=45787293/1000000,linetype="dashed", color="red2", linewidth=0.5) +

  geom_point(size=0.1) +
  ylab("Fst") +
  xlab("Position (MB)") + 
  theme(axis.text.x = element_text(vjust = 1, size = 14), plot.title=element_text(size = 18, face="bold"), legend.position="none", axis.text.y = element_text(size = 16))

recombination_data_LG5 <- read.table("iM_Zhang_2024_recombination_bins.txt", sep = '\t', header = TRUE)
recombination_data_LG5 <- recombination_data_LG5[recombination_data_LG5$SNPs >= 5, ] 

counter <- 1
usable_heights <- c(-0.05, -0.07, -0.09, -0.11, -0.13)
max_levels <- length(usable_heights)

for (i in 1:nrow(recombination_data_LG5)) {
  start_bin <- recombination_data_LG5[i,2]/1000000
  len_bin <- recombination_data_LG5[i,3]/1000000
  end_bin <- start_bin + len_bin
  
  plot_height <- usable_heights[counter]
  plot_height
  grafica_fst_zhang <- grafica_fst_zhang + 
    annotate("segment",x=start_bin, y= plot_height, xend=end_bin, yend=plot_height, color="orange2", linewidth=1.2)
    
  counter <- counter + 1
  if(counter > max_levels) {
    counter <- 1
  }
}

print_final_plot <- plot_grid(grafica_fst_LG5, grafica_fst_zhang, align= "v", ncol=1)
ggsave("Sup_File_Fst_Zhang_vs_Us_Chr5.pdf", plot=print_final_plot, width=10, height=14) 
ggsave("Sup_File_Fst_Zhang_vs_Us_Chr5.png", plot=print_final_plot, width=10, height=14) 



#############################################################################
#############################################################################
#############################################################################

data_Fst <- read.table("Good_Mean_Fsts_Windows_10000_Freqs_Zhang_iM_2024_masked_GWAS_RS.pileup.txt", sep = '\t', header = TRUE)
data_Fst <- data_Fst[data_Fst$Sites >= 5, ]
data_TMHMM_Full <- read.table("TMHMM_on_window_anot.txt", sep = '\t', header = TRUE)


Scaffold <- "CM074207.1"
pre_subdata_fst <- data_Fst[data_Fst$Scaffold == Scaffold, ]
megabases <- as.data.frame(c(pre_subdata_fst$Window)/1000000)
colnames(megabases) <- "MB"
pre_subdata_fst <- cbind(pre_subdata_fst,megabases)

pre_subdata_TMHMM <- data_TMHMM_Full[data_TMHMM_Full$Chr == Scaffold, ]
megabases <- as.data.frame(c(((pre_subdata_TMHMM$Start + ((pre_subdata_TMHMM$End - (pre_subdata_TMHMM$Start +1) )/2)))/1000000))
colnames(megabases) <- "MB"
pre_subdata_TMHMM <- cbind(pre_subdata_TMHMM, megabases)
x_max <- max(c(pre_subdata_fst$MB, pre_subdata_TMHMM$MB ))
start_region <- 8827438/1000000
end_region <- 9252765/1000000

title <- paste("Fst across ", Scaffold, sep= "")

grafica_fst_local <- ggplot(pre_subdata_fst, aes(x = MB, y = Fst))  +
  theme_classic() +
  ggtitle(title) +
  geom_point(size=0.1) +
  geom_vline(xintercept=start_region, linetype="dashed", color="red2", linewidth=0.5) +
  geom_vline(xintercept=end_region, linetype="dashed", color="red2", linewidth=0.5) +
  ylab("Fst") +
  xlab(NULL) +
  xlim(0, x_max) +
  scale_color_manual(values=color_for_plot) + 
  theme(axis.text.x=element_blank(), axis.ticks.x = element_line(color = "white"), axis.line.x =element_line(color = "white"), axis.text.y = element_text(size = 14), plot.title = element_text(size = 18), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))

title <- paste("Nº TM1 Genes in iM_Zhong: ",Scaffold , sep="" )
grafica_tmhmm_local <- ggplot(pre_subdata_TMHMM , aes(x = MB, y = N_Genes_Single_TMHMM, color=Annot_source, group=Annot_source))  +
  theme_classic() +
  geom_line(linewidth = 1) +
  xlab("Position (MB)") + 
  ylab("Nº TM1 genes") +
  geom_vline(xintercept=start_region, linetype="dashed", color="red2", linewidth=0.5) +
  geom_vline(xintercept=end_region, linetype="dashed", color="red2", linewidth=0.5) +
  xlim(0, x_max) +
  theme(axis.text.x = element_text(vjust = 1, size = 14), axis.text.y = element_text(size = 14), legend.position="bottom", plot.title = element_text(size = 18), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))

print_final_plot <- plot_grid(grafica_fst_local, grafica_tmhmm_local, align= "v", ncol=1,rel_heights = c(1,2))  
ggsave("GRC_TMHMM_plots.pdf", plot=print_final_plot, width=12, height=6) 

###################################################################################

Scaffold <- "CM074215.1"
pre_subdata_fst <- data_Fst[data_Fst$Scaffold == Scaffold, ]
megabases <- as.data.frame(c(pre_subdata_fst$Window)/1000000)
colnames(megabases) <- "MB"
pre_subdata_fst <- cbind(pre_subdata_fst,megabases)

pre_subdata_TMHMM <- data_TMHMM_Full[data_TMHMM_Full$Chr == Scaffold, ]
megabases <- as.data.frame(c(((pre_subdata_TMHMM$Start + ((pre_subdata_TMHMM$End - (pre_subdata_TMHMM$Start +1) )/2)))/1000000))
colnames(megabases) <- "MB"
pre_subdata_TMHMM <- cbind(pre_subdata_TMHMM, megabases)
x_max <- max(c(pre_subdata_fst$MB, pre_subdata_TMHMM$MB ))
start_region <- 28686194/1000000
end_region <- 30287776/1000000

title <- paste("Fst across ", Scaffold, sep= "")
grafica_fst_local <- ggplot(pre_subdata_fst, aes(x = MB, y = Fst))  +
  theme_classic() +
  ggtitle(title) +
  geom_point(size=0.1) +
  geom_vline(xintercept=start_region, linetype="dashed", color="red2", linewidth=0.5) +
  geom_vline(xintercept=end_region, linetype="dashed", color="red2", linewidth=0.5) +
  ylab("Fst") +
  xlab(NULL) +
  xlim(0, x_max) +
  scale_color_manual(values=color_for_plot) + 
  theme(axis.text.x=element_blank(), axis.ticks.x = element_line(color = "white"), axis.line.x =element_line(color = "white"), axis.text.y = element_text(size = 14), plot.title = element_text(size = 18), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))

title <- paste("Nº TM1 Genes in iM_Zhong: ",Scaffold , sep="" )
grafica_tmhmm_local <- ggplot(pre_subdata_TMHMM , aes(x = MB, y = N_Genes_Single_TMHMM, color=Annot_source, group=Annot_source))  +
  theme_classic() +
  geom_line(linewidth = 1) +
  xlab("Position (MB)") + 
  ylab("Nº TM1 genes") +
  geom_vline(xintercept=start_region, linetype="dashed", color="red2", linewidth=0.5) +
  geom_vline(xintercept=end_region, linetype="dashed", color="red2", linewidth=0.5) +
  xlim(0, x_max) +
  theme(axis.text.x = element_text(vjust = 1, size = 14), axis.text.y = element_text(size = 14), legend.position="bottom", plot.title = element_text(size = 18), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))

print_final_plot <- plot_grid(grafica_fst_local, grafica_tmhmm_local, align= "v", ncol=1,rel_heights = c(1,2))  
ggsave("PTC2_TMHMM_plots.pdf", plot=print_final_plot, width=12, height=6) 
