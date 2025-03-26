#!/usr/bin/env Rscript

# To accept arguments from command line

args = commandArgs(trailingOnly=TRUE)



# Loading required libraries

library(ggplot2)
library(stringi)



# Loading csv file with data to plot

data <- read.csv(args[1], header=FALSE)



# Setting the output directory

outdir <- args[2]



# Function to plot the reference genome (as a red rectangle)

plot_ref <- function(x){

  g <- ggplot(data=x)+geom_rect(data=x, mapping=aes(xmin=0, xmax=V4, ymin=-0.15, ymax=0.15),color="black",fill="red")+
    theme_minimal()+ggtitle(x$V1)+scale_x_continuous(position = "top")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length = unit(.25, "cm"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line = element_line(color="black"),
          axis.line.y = element_blank(),
          legend.position = "none", 
          axis.text = element_text(color="black", size=14) )
  
  return(g)
}



# Function to plot unmapped regions of contigs

plot_not <- function(g,x1,x2,ypos){

  g <- g+geom_segment(x=x1,xend=x2,y=ypos,yend=ypos)

  print(g)

  return(g)
}



# Function to plot mapped regions of contigs

plot_mapped <- function(g,x1,x2,ypos,fill,alpha){

  g <- g+geom_rect(aes(xmin=x1,xmax=x2,ymin=ypos-0.15,ymax=ypos+0.15),
                   color="black",fill=fill,alpha=alpha)

  print(g)

  return(g)
}



# Function to plot whole contigs (x is the row of the sample_aln.csv file)

plot_contig <- function(g,x,ypos){
  
  cigar <- as.numeric(strsplit(x[[6]],split=",")[[1]])
  
  if(x[[3]]=="primary"){
    color = "blue"
    alpha = 1
  }else if(x[[3]]=="not primary"){
    color = "orange"
    alpha = 0.5
  }else if(x[[3]]=="supplementary"){
    color = "pink"
    alpha = 0.75
  }
  
  actual = cigar[1]
  last = x[[5]]
  
  for(i in 2:length(cigar)){

    if(i==2 & actual == 0){
      if(last-cigar[i]<0){
        g <- plot_not(g,last,0,ypos)
        actual=1
      }else{
        g <- plot_not(g,last,last-cigar[i],ypos)
        actual=1
      }
      next
    }

    next_pos = last+(cigar[i])
    
    if(actual==0){
      if(next_pos>x[[4]]){
        if(i==length(cigar)){
          g <- plot_not(g,last,x[[4]],ypos)
        }else if(i<length(cigar)){
          g <- plot_not(g,last,x[[4]],ypos)
          g <- plot_not(g,0,next_pos-x[[4]],ypos)
          last = next_pos-x[[4]]
          actual=1
        }
      }else{
        g <- plot_not(g,last,next_pos,ypos)
        last = next_pos
        actual = 1
      }
    }else if(actual==1){
      if(next_pos>x[[4]]){
        g <- plot_mapped(g,last,x[[4]],ypos,color,alpha)
        g <- plot_mapped(g,0,next_pos-x[[4]],ypos,color,alpha)
        last = next_pos-x[[4]]
        actual = 0
      }else{
        g <- plot_mapped(g,last,next_pos,ypos,color,alpha)
        last = next_pos
        actual = 0
      }
    }
  }

  return(g)
  
}



# Defining function for plotting all the contigs mapped to viral genome

plot_aln <- function(df){

  g <- plot_ref(df[1,])

  for(i in 1:dim(df)[1]){

    g <- plot_contig(g,df[i,],i*(-0.4))

  }
  
  return(g)
  
}



# Defining function for adding annotations of contigs

plot_annot <- function(g,x,ypos){

  if(!is.na(x[[9]])){
    annot <- strsplit(x[[9]], split=",")[[1]]
  } else {
    annot <- character(0)  # Empty vector
  }

  for(i in 1:(length(annot)/3)){

    x1 <- as.integer(annot[i*3-1])
    x2 <- as.integer(annot[i*3])

    if(x[[7]]!="True"){
      x1 <- x[[8]]-x1
      x2 <- x[[8]]-x2
    }

    cigar <- as.numeric(strsplit(x[[6]],split=",")[[1]])

    if(cigar[1]==0){
      x1 <- x1-cigar[2]
      x2 <- x2-cigar[2]
    }

    if(x[[5]]+x1>0 & x[[5]]+x2>0){
      g <- g+annotate("segment", x = x[[5]]+x1, xend = x[[5]]+x2,
                  y = ypos, yend = ypos,
                  arrow = arrow(ends = "both", angle = 90, length = unit(.1,"cm")))+
      annotate("text", x = x[[5]]+((x1+x2)/2), y = ypos-0.125, label = annot[i*3-2])
    }else if(x[[5]]+x1<0 & x[[5]]+x2>0){
      g <- g+annotate("segment", x = 0, xend = x[[5]]+x2,
                  y = ypos, yend = ypos,
                  arrow = arrow(ends = "both", angle = 90, length = unit(.1,"cm")))+
      annotate("text", x = (x[[5]]+x2)/2, y = ypos-0.125, label = annot[i*3-2])
    }
 

    
  }
  
  return(g)
}



# Defining function for plotting all contigs mapped to a viral genome and their annotations

plot_aln_annot <- function(df){

  g <- plot_ref(df[1,])
  
  for(i in 1:dim(df)[1]){
    g <- plot_contig(g,df[i,],i*(-0.6))
    if(stri_length(df[i,9])>1){
      g <- plot_annot(g,df[i,],i*(-0.6)-0.25)
    }
  }
  
  return(g)
  
}



# Number of different viral genomes with contigs mapped to them

lev <- levels(as.factor(data$V1))



# Iterate through each viral genome and plot with and without annotations

if(!"0"%in%lev){

  for(i in 1:length(lev)){
    
    temp <- data[which(data[,1]==lev[i]),]
    
    g_temp <- plot_aln(temp)

    ggsave(paste(outdir,lev[i],"_mapped.jpg",sep=""),g_temp,width = 1200,height = 800,units = "px",dpi=200)


    if(any(sapply(temp[,9],stri_length)>1)){

      g_temp_annot <- plot_aln_annot(temp)

      ggsave(paste(outdir,lev[i],"_mapped_annot.jpg",sep=""),g_temp_annot,width = 1200,height = 800,units = "px",dpi=200)
    }
  }
}else{
  g <- ggplot(data, aes(x=0,y=0))+theme_void()+
  annotate("text",x=0,y=0,label="No viral sequence found",size=10)

  ggsave(paste(outdir,"NO_VIRUS.jpg",sep=""),g,width = 1200,height = 800,units = "px",dpi=200)
}

