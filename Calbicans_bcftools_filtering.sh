#!/bin/bash

#SBATCH --job-name=variantCalling
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=2
#SBATCH --time=03:00:00
#SBATCH -p 
#SBATCH -o variant_%j.output
#SBATCH -e variant_%j.error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --array=1

set -ue
set -o pipefail

array=${SLURM_ARRAY_TASK_ID} # If working on a cluster, change to the appropriate variable

species="Calbicans"
split_vcf_raw= # File with one Mutect2 filtered vcf per line
output= # Name to give the output .vcf file

snpEffDB=SC5314_s02m09r08 # Requires preconfigured snpEff database; refer to snpEff documentation!
filterArgs="(ANN[*].IMPACT has 'HIGH' | ANN[*].IMPACT has 'MODERATE')"

vcftools=/home/selmecki/shared/software/VCFtools/bin/vcftools

#load bcftools
module load bcftools

# merge if different samples
bcftools merge --force-samples -l ${split_vcf_raw} -o ${output}_concat.vcf

###Filter with bcftools

# For Mutect2

bcftools view -O u "${output}_concat.vcf" \
| bcftools view \
-i "INFO/MMQ>=40 & FORMAT/AF[0]>=0.2 & FORMAT/AD[*]>4 & FORMAT/F1R2[*]>=1 & FORMAT/F2R1[*]>0" \
| bcftools filter -G 10 -o "${output}_filtered.vcf"

bcftools view ${output}_filtered.vcf > ${output}_vcfTools.recode.vcf

###Annotate with SnpEff & SnpSift

module load java/openjdk-21.0.2 

java -Xmx4g -jar "/home/selmecki/shared/software/snpEff/snpEff.jar" \
-c "/home/selmecki/shared/software/snpEff/snpEff.config" \
"${snpEffDB}" \
"${output}_vcfTools.recode.vcf" > "${output}_snpeff.vcf"
 

java -Xmx4g -jar "/home/selmecki/shared/software/snpEff/SnpSift.jar" filter \
"${filterArgs}" "${output}_snpeff.vcf" > "${output}.vcf"

#Index with bcftools
module load htslib
bgzip "${output}.vcf"
bcftools index -t "${output}.vcf.gz"
bcftools stats "${output}.vcf.gz" > "${output}.vcf.stats"
