Illumina whole genome sequencing (array_bbduk_align_sort.sh)

Genomic DNA was isolated using a phenol-chloroform extraction as described previously. Libraries were prepared using the Illumina DNA Prep kit and IDT 10bp UDI indices, and sequenced on an Illumina NextSeq 2000, producing 2x151bp reads. Demultiplexing, quality control, and adapter trimming were performed with bcl-convert (https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html)(v3.9.3). Adapter and quality trimming were performed with BBDuk (BBTools v38.94. Trimmed reads were aligned to the C. albicans reference genome (SC5314_version_A21-s02-m09-r08) using BWA-MEM (v0.7.17) with default parameters. Aligned reads were sorted, duplicate reads were marked and the resulting BAM file was indexed with Samtools (v1.10). Quality of trimmed FASTQ and BAM files was assessed for all strains with FastQC (v0.11.7), Qualimap (v2.2.2-dev) and MultiQC (v1.16)

Variant Calling and SNP Analysis
(i) In vitro evolved isolates (Calbicans_Mutect2.sh)

De novo variant calling was performed using GATK (v.4.1.2) as previously described. Briefly, variant calling was performed with Mutect2, assigning each progenitor as ‘normal’ and the evolved isolates as ‘tumour’. After filtering, merged VCF files were created for each progenitor group followed by additional hard filtering. Variants were annotated with SnpEff (v.5.0e; database built from SC5314 v.A21-s02-m09-r08, with alternate yeast nuclear codon table). 

(ii) Clinical and environmental isolates (Calbicans_bcftools_filtering.sh)

WGS data for 101 bloodstream isolates and 199 publicly available data sets (BioProjects PRJNA193498 and PRJNA432884) were processed as described above. SNPs and small indels were called using HaplotypeCaller (GATK v4.4.0). SNPs were filtered on the parameters QD < 2.00, QUAL < 30.0, SOR > 3.0, FS > 60.0, MQ < 40.0, MQRankSum < -12.5, and ReadPosRankSum < -8.0. Indels were filtered on the parameters QD < 2.0, QUAL < 30.0, FS > 200.0, and ReadPosRankSum < -20.0. Bcftools (v1.17) was used to calculate variant allele frequency (VAF) per sample, filter for heterozygous variants with VAF between 0.15 and 0.85 and homozygous variants with VAF > 0.98, and to exclude known repetitive regions as annotated in the SC5314 A21-s02-m09-r08 GFF (rRNA, repeat_region, retrotransposon) and telomere-proximal regions, defined here as extending from each chromosome end to the first non-repetitive genome feature. Filtered variants were annotated using SnpEff. CAP1 variants were compiled into Table S2 using GATK VariantsToTable and bcftools query.

Data visulization figure 1B (allele_frequencies.R)
