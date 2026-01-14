# run in command line:
# bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%AF\t%ANN\n]' Calbicans_Mutect2_5178.vcf.gz > variant_table.txt

library(vcfR)
library(tidyverse)

# Read in and format variant table from bcftools

table <- read.table("variant_table.txt")
samples_table <- table[!grepl("AMS5178",table[,1]),]
colnames(samples_table) <- c("SAMPLE","CHROM","POS","REF","ALT","GT","AF","ANN")

samples_filtered_1 <- samples_table %>%
  subset(!grepl(",\\.",AF) & !grepl("\\.,",AF) & AF!=".")

samples_filtered <- samples_filtered_1[which(sapply(samples_filtered_1[,"AF"],function(x)all(unlist(lapply(str_split(x,","),function(y)as.numeric(y)>=0.2))))),]

ann <- lapply(samples_filtered$ANN,function(x){
  split_info <- str_split(x,"\\,")[[1]]
  split_split <- lapply(split_info,function(i)str_split(i,"\\|")[[1]])
  refined_info <- split_split[which(sapply(split_split,function(z)z[3]!="MODIFIER"))]
  refined_info
})
ann_table <- do.call(rbind,unlist(ann,recursive=F))


# Expand the SNP table to account for multiple alleles at the same position

expanded_samples <- apply(samples_filtered,1,function(x)data.frame(c(str_split(x[1:7],"\\,")),x[8]))

expanded.dataframe <- do.call(rbind, lapply(expanded_samples,function(x){  export <- t(matrix(unlist(x),nrow=8,byrow=T)) 
cbind(export,c(1:nrow(export)))
}))

# Fix annotation for each entry

ann_split <- apply(expanded.dataframe,1,function(x){
  ann_sample <- str_split(x[8],"\\,")
  ann_split <- lapply(unlist(ann_sample),function(i)str_split(i,"\\|")[[1]])
  refined_split <- rev(ann_split[which(sapply(ann_split,function(z)z[3]!="MODIFIER"))])
  query <- which(unlist(lapply(refined_split,function(z)z[1]==x[5])))
  if(length(query)!=0){
    ann_export <- refined_split[[query]]
  } else {ann_export <- rep("remove",16)}
  
  ann_export
})

ann_matrix <- do.call(rbind,ann_split)
expanded.variants <- cbind(expanded.dataframe[,1:7],ann_matrix)

ann_header <- "Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos / cDNA.length|CDS.pos / CDS.length|AA.pos / AA.length|Distance|ERRORS / WARNINGS / INFO"
ann_header <- unlist(str_split(ann_header,"\\|"))
colnames(expanded.variants) <- c("SAMPLE","CHROM","POS","REF","ALT","GT","AF",ann_header)

annotated_table <- data.frame(expanded.variants)
annotated_table$Gene_Name <- gffSplit$gene_name[match(annotated_table$Gene_ID,gffSplit$gene_id)]

snp_table_export <- data.frame(cbind(annotated_table[,c(1,1)],"AMS5178",annotated_table[,c("CHROM","POS","Gene_Name","Gene_ID","REF","ALT","GT","AF","Annotation","Annotation_Impact","HGVS.c","HGVS.p")]))
colnames(snp_table_export) <- c("AMS ID","Strain Number","Progenitor","CHROM","POS","GENE","ORF","REF","ALT","GT","AF","ANNOTATION","IMPACT","CODING CHANGE","AA CHANGE")

snp_table <- snp_table_export

conversion <- "AMS6442\tChr3_lineage1_P0
AMS6443\tChr3_FLC1_L1_P1
AMS6444	Chr3_FLC1_L1_P2
AMS6445	Chr3_FLC1_L1_P3
AMS6446	Chr3_FLC1_L1_P1_S4
AMS6447	Chr3_FLC1_L1_P1_S12
AMS6448	Chr3_FLC1_L1_P1_S16
AMS6449	Chr3_FLC1_L1_P1_S18
AMS6450	Chr3_FLC1_L1_P1_S20
AMS6451	Chr3_FLC1_L1_P1_S21
AMS6452	Chr3_FLC1_L1_P1_S22"

conversion_dict <- matrix(unlist(str_split(unlist(str_split(conversion, "\n")),"\t"),recursive=F),ncol=2,byrow=T)

snp_table$`Strain Number`<- conversion_dict[match(data.frame(snp_table)[,"AMS.ID"],conversion_dict[,1]),2]
snp_table$`Strain Number`[which(is.na(snp_table$`Strain Number`))] <- snp_table$`AMS ID`[which(is.na(snp_table$`Strain Number`))]

writexl::write_xlsx(snp_table,"Mutect2_AMS5178_evolved_strains_240725.xlsx")
writexl::write_xlsx(snp_table %>% subset(GENE == "CAP1"),"cap1_frequencies_Mutect2_AMS5178_240725.xlsx")
