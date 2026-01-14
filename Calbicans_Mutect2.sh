#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --time=2:00:00
#SBATCH -p 
#SBATCH -o %x_%u_%A_%a.out
#SBATCH -e %x_%u_%A_%a.err
#SBATCH --array=1-15


line=${SLURM_ARRAY_TASK_ID} # If working on a cluster, change to the appropriate variable

reference_fasta= # Reference genome .fasta file

normal_bam= # Alignment .bam file to use as the "normal" sample
normal_strain= # Name to give the "normal" sample

sample_file= # File with sample information, with one bam per line. Bams are titled "${strain_name}"_*.bam
tumor_bam=$(awk -v val=$line 'NR == val { print $1}' $sample_file)
strain=$(basename "$tumor_bam" | cut -d "_" -f 1)


module load gatk/4.1.2

mkdir raw_vcfs_${normal_strain}

gatk Mutect2 --input ${tumor_bam} --input ${normal_bam} --normal ${normal_strain} --output ${strain}_Mutect2.vcf --reference ${reference_fasta}

gatk FilterMutectCalls --variant ${strain}_Mutect2.vcf --reference ${reference_fasta} --output raw_vcfs_${normal_strain}/${strain}_filtered_Mutect2.vcf

module load bcftools
module load htslib

bgzip raw_vcfs_${normal_strain}/${strain}_filtered_Mutect2.vcf
bcftools index raw_vcfs_${normal_strain}/${strain}_filtered_Mutect2.vcf.gz
