#!/bin/bash
#SBATCH -M teach
#SBATCH -A hugen2072-2025s
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=LAA196@pitt.edu
#SBATCH -c 24
#SBATCH -t 120:00:00
set -ve

#module to load
module load fastqc/0.11.9
module load gatk/4.5.0.0
module load cutadapt/2.10
module load gcc/8.2.0 bwa/0.7.17 samtools/1.21

cd $SLURM_SCRATCH

# Quality check
fastqc $SLURM_SUBMIT_DIR/p5/p5_1.fastq.gz -t 48 --outdir=.
fastqc $SLURM_SUBMIT_DIR/p5/p5_2.fastq.gz -t 48 --outdir=.

# Adapter trimming
cutadapt -j 0 -m 10 -q 20 $SLURM_SUBMIT_DIR/p5/p5_1.fastq.gz $SLURM_SUBMIT_DIR/p5/p5_2.fastq.gz \
-a AGATCGGAAGAG -A AGATCGGAAGAG \
-o $SLURM_SUBMIT_DIR/p5_1_trimmed.fastq.gz -p $SLURM_SUBMIT_DIR/p5_2_trimmed.fastq.gz

fastqc $SLURM_SUBMIT_DIR/p5_1_trimmed.fastq.gz -t 48 --outdir=$SLURM_SUBMIT_DIR/
fastqc $SLURM_SUBMIT_DIR/p5_2_trimmed.fastq.gz -t 48 --outdir=$SLURM_SUBMIT_DIR/

# Alignment
bwa mem -t 48 \
$SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
$SLURM_SUBMIT_DIR/p5_1_trimmed.fastq.gz \
$SLURM_SUBMIT_DIR/p5_2_trimmed.fastq.gz \
-R "@RG\tID:P5\tLB:P5\tSM:P5\tPL:ILLUMINA" | \
samtools view -bh | samtools sort > $SLURM_SUBMIT_DIR/project.bam
samtools index $SLURM_SUBMIT_DIR/project.bam
rm  $SLURM_SUBMIT_DIR/p5_1_trimmed.fastq.gz $SLURM_SUBMIT_DIR/p5_2_trimmed.fastq.gz


run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/project.bam" ]; then
        echo "Script failed before completion — restoring project.bam"
        cp project.bam* $SLURM_SUBMIT_DIR/
    fi
}
trap run_on_exit EXIT

cp $SLURM_SUBMIT_DIR/project.bam* . && rm $SLURM_SUBMIT_DIR/project.bam*


# Alignment Quality Control
#  MarkDuplicatesSpark
gatk MarkDuplicatesSpark -I $SLURM_SCRATCH/project.bam \
-O $SLURM_SUBMIT_DIR/project_dupsmarked.bam

if [ $? -eq 0 ]; then
    echo "MarkDuplicatesSpark finished successfully. Copying result to submit dir..."
else
    echo "MarkDuplicatesSpark failed."
fi


run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/project_dupsmarked.bam" ]; then
        echo "Script failed before completion — restoring project_dupsmarked.bam"
        cp project_dupsmarked.bam* $SLURM_SUBMIT_DIR/
    fi
}
trap run_on_exit EXIT
cp $SLURM_SUBMIT_DIR/project_dupsmarked.bam* . && rm $SLURM_SUBMIT_DIR/project_dupsmarked.bam*

#  Base quality recalibration
gatk BaseRecalibrator -I $SLURM_SCRATCH/project_dupsmarked.bam \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
--known-sites $SLURM_SUBMIT_DIR/p5/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known-sites $SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
-O $SLURM_SUBMIT_DIR/project_BQSR.table

if [ $? -eq 0 ]; then
    echo "Base quality recalibration finished successfully. Copying result to submit dir..."
    cp $SLURM_SCRATCH/project_dupsmarked.bam* $SLURM_SUBMIT_DIR/
else
    echo "Base quality recalibration failed."
fi

run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/project_dupsmarked.bam" ]; then
        echo "Script failed before completion — restoring project_dupsmarked.bam"
        cp project_dupsmarked.bam* $SLURM_SUBMIT_DIR/
    fi
}
trap run_on_exit EXIT
cp $SLURM_SUBMIT_DIR/project_dupsmarked.bam* . && rm $SLURM_SUBMIT_DIR/project_dupsmarked.bam*
#  ApplyBQSR
gatk ApplyBQSR -R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-I $SLURM_SCRATCH/project_dupsmarked.bam \
--bqsr-recal-file $SLURM_SUBMIT_DIR/project_BQSR.table \
-O $SLURM_SCRATCH/project_dupsmarked_cleaned.bam

if [ $? -eq 0 ]; then
    echo "ApplyBQSR finished successfully. Copying result to submit dir..."
    cp $SLURM_SCRATCH/project_dupsmarked_cleaned.bam* $SLURM_SUBMIT_DIR/
else
    echo "ApplyBQSR failed."
fi

rm $SLURM_SUBMIT_DIR/project_BQSR.table

samtools flagstat -@ 48 $SLURM_SUBMIT_DIR/project_dupsmarked_cleaned.bam \
                  > $SLURM_SUBMIT_DIR/project_alignment_statistics.out


samtools depth $SLURM_SUBMIT_DIR/project_dupsmarked_cleaned.bam | gzip $SLURM_SUBMIT_DIR/project_depth_statistics.out.gz

$SLURM_SUBMIT_DIR/project_depth_statistics.out.gz awk '{if ($3 >= 0) {total_bases += $3}} END {print total_bases / NR}' \
> $SLURM_SUBMIT_DIR/project_avg_read_depth.out

$SLURM_SUBMIT_DIR/project_depth_statistics.out.gz | awk '{sum += $3} END {print sum}' > $SLURM_SUBMIT_DIR/project_total_depth.out

rm $SLURM_SUBMIT_DIR/project_depth_statistics.out.gz

run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/project_dupsmarked_cleaned.bam" ]; then
        echo "Script failed before completion — restoring project_dupsmarked.bam"
        cp project_dupsmarked_cleaned.bam* $SLURM_SUBMIT_DIR/
    fi
}
trap run_on_exit EXIT
cp $SLURM_SUBMIT_DIR/project_dupsmarked_cleaned.bam* . && rm $SLURM_SUBMIT_DIR/project_dupsmarked_cleaned.bam*

# Genotyping
#  HaplotypeCaller
gatk HaplotypeCaller -I $SLURM_SCRATCH/project_dupsmarked_cleaned.bam \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-O $SLURM_SUBMIT_DIR/project_genotypes.g.vcf.gz \
-ERC GVCF \
-OVI \
--native-pair-hmm-threads 48

if [ $? -eq 0 ]; then
    echo "HaplotypeCaller finished successfully. Continuing with GenotypeGVCFs"
else
    echo "HaplotypeCaller failed."
    cp $SLURM_SCRATCH/project_dupsmarked_cleaned.bam* $SLURM_SUBMIT_DIR/
fi

run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/project_genotypes.g.vcf.gz*" ]; then
        echo "Script failed before completion — restoring project_genotypes.g.vcf.gz"
        cp project_genotypes.g.vcf.gz* $SLURM_SUBMIT_DIR/
    fi
}
trap run_on_exit EXIT

cp $SLURM_SUBMIT_DIR/project_genotypes.g.vcf.gz* . && rm $SLURM_SUBMIT_DIR/project_genotypes.g.vcf.gz*

#  GenotypeGVCFs
gatk GenotypeGVCFs -R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V $SLURM_SCRATCH/project_genotypes.g.vcf.gz \
-O $SLURM_SCRATCH/project_genotypes.vcf.gz

if [ $? -eq 0 ]; then
    echo "GenotypeGVCFs finished successfully. Continuing with GenotypeGVCFs"
    cp $SLURM_SCRATCH/project_genotypes.vcf* $SLURM_SUBMIT_DIR/
else
    echo "GenotypeGVCFs failed."
fi

# Genotype Quality Control
gatk MakeSitesOnlyVcf -I $SLURM_SUBMIT_DIR/project_genotypes.vcf.gz \
-O $SLURM_SUBMIT_DIR/project_sites_only.vcf.gz

if [ $? -eq 0 ]; then
    echo "Genotype Quality Control finished successfully." 
else
    echo "Genotype Quality Control failed."
fi

# VariantRecalibrator - INDEL
gatk VariantRecalibrator -mode INDEL \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V $SLURM_SUBMIT_DIR/project_sites_only.vcf.gz \
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
-resource:axiomPoly,known=false,training=true,truth=false,prior=10 \
	$SLURM_SUBMIT_DIR/p5/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=2 \
	$SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
-resource:mills,known=false,training=true,truth=true,prior=12 \
	$SLURM_SUBMIT_DIR/p5/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--java-options "-XX:ParallelGCThreads=48" \
-O $SLURM_SUBMIT_DIR/project_indels.recal \
--tranches-file $SLURM_SUBMIT_DIR/project_indels.tranches 

# VariantRecalibrator - SNP
gatk VariantRecalibrator -mode SNP \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V $SLURM_SUBMIT_DIR/project_sites_only.vcf.gz \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
-resource:hapmap,known=false,training=true,truth=true,prior=15 \
	$SLURM_SUBMIT_DIR/p5/hapmap_3.3.hg38.vcf.gz \
-resource:omni,known=false,training=true,truth=true,prior=12 \
	$SLURM_SUBMIT_DIR/p5/1000G_omni2.5.hg38.vcf.gz \
-resource:1000G,known=false,training=true,truth=false,prior=10 \
	$SLURM_SUBMIT_DIR/p5/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=7 \
	$SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
--java-options "-XX:ParallelGCThreads=48" \
-O $SLURM_SUBMIT_DIR/project_snps.recal \
--tranches-file $SLURM_SUBMIT_DIR/project_snps.tranches


# ApplyVQSR - INDEL
gatk ApplyVQSR -mode INDEL \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V $SLURM_SUBMIT_DIR/project_genotypes.vcf.gz \
--recal-file $SLURM_SUBMIT_DIR/project_indels.recal \
--tranches-file $SLURM_SUBMIT_DIR/project_indels.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
--java-options "-XX:ParallelGCThreads=48" \
-O $SLURM_SUBMIT_DIR/project_genotypes_indelqc.vcf.gz

# ApplyVQSR - SNP
gatk ApplyVQSR -mode SNP \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V $SLURM_SUBMIT_DIR/project_genotypes_indelqc.vcf.gz \
--recal-file $SLURM_SUBMIT_DIR/project_snps.recal \
--tranches-file $SLURM_SUBMIT_DIR/project_snps.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
--java-options "-XX:ParallelGCThreads=48" \
-O $SLURM_SUBMIT_DIR/project.vcf.gz

# Metrics
gatk CollectVariantCallingMetrics \
-I $SLURM_SUBMIT_DIR/project.vcf.gz \
--DBSNP $SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
-O $SLURM_SUBMIT_DIR/project_genotype_metrics
