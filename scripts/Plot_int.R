#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(RIdeogram)
library(stringr)

int <- read.csv(args[1])

outdir <- args[2]

data(human_karyotype, package="RIdeogram")

gene_density <- read.csv("scripts/bands.txt")

for_plot <- int[int[,2]!="triangle",]

rows <- sample(nrow(for_plot), 300)

for_plot <- rbind(for_plot[rows,],int[int[,2]=="triangle",])

ideogram(karyotype = human_karyotype, overlaid = gene_density, label = for_plot,
         label_type = "marker",colorset1 = c("#ffffff", "#000000"))

convertSVG("chromosome.svg", device = "png")

file.rename(from="chromosome.png",to=paste(outdir,"all.png"))

file.rename(from="chromosome.svg",to=paste(outdir,"all.svg"))

ideogram(karyotype = human_karyotype, overlaid = gene_density, label = int[int[,2]=="triangle",],
         label_type = "marker",colorset1 = c("#ffffff", "#000000"))

convertSVG("chromosome.svg", device = "png")

file.rename(from="chromosome.png",to=paste(outdir,"viral.png"))

file.rename(from="chromosome.svg",to=paste(outdir,"viral.svg"))
