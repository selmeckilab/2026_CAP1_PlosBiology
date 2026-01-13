
library(ggplot2)
library(tidyverse)
library(readxl)
library(writexl)

cap1_freqs <- data.frame(readxl::read_xlsx("/Volumes/Xin/CAP1/allele frequency/recap1predictions/cap1_frequencies_Mutect2_AMS5178_240725.xlsx"))

spacer=1

cap1_freqs_plot <- rbind(cap1_freqs[1:4,],cap1_freqs[rep(5,spacer+1),])

cap1_freqs_plot[c(4+1:spacer),c("AMS.ID","Strain.Number")] <- c(paste0("spacer",1:spacer),rep("",spacer))

cap1_freqs_plot[5,"AF"] <- mean(c(as.numeric(cap1_freqs_plot[4,"AF"]),as.numeric(cap1_freqs_plot[6,"AF"])))

allele_frequencies_plot <- (cap1_freqs_plot) %>%
  mutate(AMS.ID=factor(AMS.ID, levels=AMS.ID)) %>% data.frame() 

allele_frequencies_plot %>%
  ggplot(aes(x=AMS.ID,y=as.numeric(AF),group=1)) +
  geom_rect(aes(xmin=0,xmax="AMS5726/AMS5727/AMS5728",ymin=0,ymax=1,fill="Progenitor")) +
  geom_area(aes(fill="Mut")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(breaks=cap1_freqs_plot$AMS.ID,labels=cap1_freqs_plot$Strain.Number) +
  scale_fill_manual(values=c("Mut"=rgb(0,191,196,maxColorValue = 255),"Progenitor"=rgb(128,130,133,maxColorValue = 255)), breaks=c("Progenitor","Mut"),name="Cap1 Genotype") +
  geom_line(colour="white",linewidth=1) +
  geom_point(data=subset(allele_frequencies_plot,!(Strain.Number=="")),colour="white",size=2) +
  labs(x="Passage",y="Allele Frequency") +
  theme(panel.background = element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=60,hjust=1,vjust=1),
        axis.ticks.x=element_line(colour=sapply(allele_frequencies_plot$Strain.Number,function(x)ifelse(x=="","transparent","black"))),
        plot.margin=unit(c(0,0,0,0),"null")) +
  coord_cartesian(xlim=c(1.6,5.5), ylim=c(0.04,1))
