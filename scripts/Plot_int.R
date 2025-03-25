#!/usr/bin/env Rscript

# To accept arguments from command line

args = commandArgs(trailingOnly=TRUE)



# Loading required libraries

library(RIdeogram)
library(stringr)
library(ggplot2)



# Loading file with information of second mapping to host genome

int <- read.csv(args[1],header = TRUE)



# Setting output directory

outdir <- args[2]



# Loading human genome and human karyotype

data(human_karyotype, package="RIdeogram")

human_karyotype <- human_karyotype[-c(24),]



# Loading cytogenetic bands file

gene_density <- read.csv("scripts/bands.txt")



# Density plot for lengths of contigs mapped to second mapping to host genome

cols = c("black","#0000ff","#ffc0cb","#ffa500")



# Normal plot

g <- ggplot(int, aes(x = len,fill=Type)) +
  geom_density(alpha = 0.5, color = "black") +
  scale_fill_manual(values = cols)+theme_minimal()+
  xlab("Contig length (bp)")+ylab("Density")+
  ggtitle("Length distribution of mapped contigs")

ggsave(paste(outdir,"contig_length.jpg",sep=""),g,width = 1200,height = 800,units = "px",dpi=200, create.dir = TRUE)



# Plot of log lengths

g_log <- ggplot(int, aes(x = log10(len), fill = Type)) +
  geom_density(alpha = 0.5, color = "black") + 
  scale_fill_manual(values = cols)+theme_minimal()+
  xlab("Contig length log10(bp)")+ylab("Density")+
  ggtitle("Length distribution of mapped contigs","Log 10 transformation")

ggsave(paste(outdir,"contig_log_length.jpg",sep=""),g_log,width = 1200,height = 800,units = "px",dpi=200, create.dir = TRUE)



# Plotting the ideograms



# Choosing 300 random contigs

for_plot <- int[int[,2]!="triangle",]

if(nrow(for_plot)>300){
    rows <- sample(nrow(for_plot), 300)
  }else{
    rows <- seq(nrow(for_plot))
  }



# Including all viral contigs mapped to host genome too

for_plot <- rbind(for_plot[rows,],int[int[,2]=="triangle",])



# Creating ideogram of all selected contigs

ideogram(karyotype = human_karyotype, overlaid = gene_density, label = for_plot,
         label_type = "marker",colorset1 = c("#ffffff", "#000000"))

convertSVG("chromosome.svg", device = "png")

file.rename(from="chromosome.png",to=paste(outdir,"all.png"))

file.rename(from="chromosome.svg",to=paste(outdir,"all.svg"))

if("triangle"%in%levels(int[,2])){



# Plotting ideogram of only viral contigs

ideogram(karyotype = human_karyotype, overlaid = gene_density, label = int[int[,2]=="triangle",],
         label_type = "marker",colorset1 = c("#ffffff", "#000000"))

convertSVG("chromosome.svg", device = "png")

file.rename(from="chromosome.png",to=paste(outdir,"viral.png"))

file.rename(from="chromosome.svg",to=paste(outdir,"viral.svg"))

}

  

# Ideograms for size category

for(i in 1:length(for_plot$len)){
  if(for_plot$len[i]<=200){
    for_plot$Type[i] <- "<200"
    for_plot$color[i] <- "021893"
  }else if(for_plot$len[i]>200&for_plot$len[i]<=500){
    for_plot$Type[i] <- "200-500"
    for_plot$color[i] <- "3537ae"
  }else if(for_plot$len[i]>500&for_plot$len[i]<=1000){
    for_plot$Type[i] <- "500-1000"
    for_plot$color[i] <- "740699"
  }else if(for_plot$len[i]>1000&for_plot$len[i]<=1500){
    for_plot$Type[i] <- "1000-1500"
    for_plot$color[i] <- "b70b0b"
  }else if(for_plot$len[i]>1500){
    for_plot$Type[i] <- ">1500"
    for_plot$color[i] <- "990623"
  }
}



# Plotting the ideogram

ideogram(karyotype = human_karyotype, overlaid = gene_density, label = for_plot,
         label_type = "marker",colorset1 = c("#ffffff", "#000000"))

convertSVG("chromosome.svg", device = "png")

file.rename(from="chromosome.png",to=paste(outdir,"size.png"))

file.rename(from="chromosome.svg",to=paste(outdir,"size.svg"))

