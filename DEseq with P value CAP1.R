library(ggplot2)
library(tximport)
library(DESeq2)
library(ggpubr)
library(EnhancedVolcano)
library(ggalt)

setwd("/Volumes/Extreme SSD/RNA_Seq/CAP1 RNA seq/quants")
samples <- read.csv("CAP1_manifest.csv", header = TRUE)
head(samples)
files <- file.path(samples$SampleID, "quant.sf")
all(file.exists(files))

## head(tx2gene)
txi <- tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE)
head(txi$counts)
#write.csv(as.data.frame(txi$counts), file = "readcounts.csv")

sampleTable <- data.frame(Name = paste0("S", 1:12), Strains = rep(c("AMS2401", "AMS6430"), each=6),
                          Condition = rep(c("FLC2","YPAD"), each=3,2))
sampleTable$Strains <-factor(sampleTable$Strains, levels = c("AMS2401","AMS6430"))
sampleTable$Condition <-factor(sampleTable$Condition, levels = c("YPAD","FLC2"))
#rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~Strains + Condition + Strains:Condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

resultsNames(dds)
# The effect of CAP1 mutation in YPAD  
res <- results(dds, contrast = c("Strains","AMS6430","AMS2401"))

write.csv(as.data.frame(res), file = "AMS6430_vs_2401_YPAD_no thereshold_newmodel.csv")
summary(res,alpha=0.01)
res2<-as.data.frame(res)
ves <- res2[res2$padj<=0.05 & abs(res2$log2FoldChange) >=1 & !is.na(res2$padj),]
write.csv(as.data.frame(ves), file = "AMS6430_vs_2401_YPAD_no thereshold_newmodel_outputlogfold1.csv")

View(ves)

# The difference in response to fluconazole between AMS6430 and AMS2401 

res <- results(dds, name = "StrainsAMS6430.ConditionFLC2")

write.csv(as.data.frame(res), file = "AMS6430_vs_2401_FLC2vsYPAD_no thereshold_newmodel.csv")
summary(res,alpha=0.01)
res2<-as.data.frame(res)
ves <- res2[res2$padj<=0.05 & abs(res2$log2FoldChange) >=1.5 & !is.na(res2$padj),]
write.csv(as.data.frame(ves), file = "AMS6430_vs_2401_FLC2vsYPAD_1.5_newmodel_output_newmodel.csv")

View(ves)

# The effect of fluconazole in CAP1 mutation 
res <- results(dds, alpha =0.05,list(c("Condition_FLC2_vs_YPAD", "StrainsAMS6430.ConditionFLC2")))
View (res)
write.csv(as.data.frame(res), file = "AMS6430_vs_AMS6430_FLC2vsYPAD_no thereshold_newmodel.csv")
summary(res,alpha=0.01)
res2<-as.data.frame(res)
ves <- res2[res2$padj<=0.05 & abs(res2$log2FoldChange) >=1.5 & !is.na(res2$padj),]
write.csv(as.data.frame(ves), file = "AMS6430_vs_AMS6430_FLC2vsYPAD_1.5_newmodel_output_newmodel.csv")

View(ves)
# The effect of fluconazole in WT mutation 
res <- results(dds, contrast= c("Condition","FLC2", "YPAD"))
write.csv(as.data.frame(res), file = "AMS2401_vs_AMS2401_FLC2vsYPAD_no thereshold_newmodel.csv")
summary(res,alpha=0.01)
res2<-as.data.frame(res)
ves <- res2[res2$padj<=0.05 & abs(res2$log2FoldChange) >=1.5 & !is.na(res2$padj),]
write.csv(as.data.frame(ves), file = "AMS2401_vs_AMS2401_FLC2vsYPAD_1.5_newmodel_output_newmodel.csv")



vc <- read.csv("/Volumes/Extreme SSD/RNA_Seq/CAP1 RNA seq/quants/AMS6430_vs_2401_FLC2vsYPAD_no thereshold_newmodel.csv", header = TRUE, row.names = 1)

EnhancedVolcano(vc,
                lab = vc$GeneNameorSymbol,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('CAP1','MDR1','SOD6','AOX2','CIP1','DBP8','OYE32','EBP1','OYE23','SOD3','CAT1','IFR1','GST1','GST1','ERG3','ERG1'),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'AMS6430 vs AMS2401_FLC2vsYPAD',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,)

library(scales)   # for trans_new

# ---- read your data ----
vc <- read.csv(
  "/Volumes/Extreme SSD/RNA_Seq/CAP1 RNA seq/quants/AMS6430_vs_2401_YPAD_no thereshold_newmodel.csv",
  header = TRUE,
  row.names = 1
)

# ---- define custom piecewise transform and its inverse ----
# Mapping (original x):
#  -10..0  -> slope 1.5  (length 15)  = 30% of total
#   0..10  -> slope 3    (length 30)  = 60% of total
#  10..20  -> slope 0.5  (length 5)   = 10% of total
stretch_fun <- function(x) {
  ifelse(x < 0,
         1.5 * x,
         ifelse(x <= 10,
                3 * x,
                30 + 0.5 * (x - 10)))
}

inv_stretch_fun <- function(y) {
  ifelse(y < 0,
         y / 1.5,             # inverse of 1.5*x
         ifelse(y <= 30,
                y / 3,         # inverse of 3*x
                10 + 2 * (y - 30)))  # inverse of 30 + 0.5*(x-10)
}

stretch_trans <- trans_new(
  name = "stretch_log2FC",
  transform = stretch_fun,
  inverse   = inv_stretch_fun
)

# ---- make volcano with ORIGINAL x so FCcutoff works ----
p <- EnhancedVolcano(
  vc,
  lab = vc$GeneNameorSymbol,
  x = 'log2FoldChange',         # keep original → FCcutoff logic stays correct
  y = 'padj',
  selectLab = c('CAP1','OYE23','CIP1','OYE32', 'MRF1', 'EBP1','IFR1',
                'MDR1','SOD1','AOX2','DBP8',
                'GST1','NRG1','BRG1', 'UPC2','ERG3','ERG11','ERG1','SLD1','ARE2', 'FLU1', 'PDR16'),
  xlab = bquote(~Log[2]~ 'fold change'),
  title = 'AMS6430 vs AMS2401_YPAD',
  pCutoff = 0.05,
  FCcutoff = 1.0,               # still cuts at log2 = ±1
  pointSize = 3.0,
  labSize = 6.0,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black',
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0
)

# ---- apply custom x scaling (space redistribution) ----
p +
  scale_x_continuous(
    trans = stretch_trans,
    limits = c(-10, 20),
    breaks = c(-10, 0, 10, 20),
    labels = c(-10, 0, 10, 20),
    expand = expansion(mult = 0.02))
 

#plotMA(res)
ma <- read.csv("/Volumes/Extreme SSD 2/RNA_Seq/2023_08 erg251/quants/AMS5853_FLC1_vs_AMS5853_YPAD_1.0.csv", header = TRUE, row.names = 1)
ma_plot <- ggmaplot(vc, main = expression("AMS5853_FLC1_vs_AMS5853_YPAD_1.0"),
                    fdr = 0.05, fc = 2, size = 3,
                    palette = c("#B31B21", "#1465AC", "darkgray"),
                    legend = "top", top = 0,
                    alpha = 0.7,
                    font.label = c("italic", 8),
                    font.legend = "bold",
                    font.main = "bold") +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
  )
ma_plot

vst <- vst(dds, blind = FALSE)
Sc_pca <- plotPCA(vst, intgroup = "Condition") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) 
Sc_pca
