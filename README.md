## End-to-End Variant Calling Pipeline (NGS, GATK, HPC)

This project implements a high-throughput, end-to-end next-generation sequencing (NGS) analysis pipeline for variant discovery using industry-standard tools. The workflow processes raw sequencing data through quality control, alignment, post-processing, variant calling, and variant quality recalibration.

---

## Project Overview

Accurate variant discovery from sequencing data requires multiple processing and quality control steps. This pipeline follows GATK Best Practices to transform raw FASTQ files into high-confidence variant calls, ensuring data quality and reproducibility in a high-performance computing (HPC) environment.

---

## Objectives

- Perform quality control on raw sequencing data  
- Trim adapters and low-quality bases  
- Align reads to the reference genome  
- Process alignments for downstream analysis  
- Perform variant calling and genotyping  
- Apply variant quality recalibration  
- Generate alignment and variant quality metrics  

---

## Workflow Summary

### Quality Control

- Performed initial QC using FastQC  
- Evaluated read quality, GC content, and sequencing artifacts  

---

### Adapter Trimming

- Removed adapter sequences and low-quality bases using Cutadapt  
- Applied minimum length and quality thresholds  

---

### Alignment

- Aligned reads to the human reference genome (GRCh38) using BWA-MEM  
- Converted SAM to sorted BAM format  
- Indexed BAM files for efficient access  

---

### Post-Alignment Processing

- Marked duplicate reads using GATK MarkDuplicatesSpark  
- Performed base quality score recalibration (BQSR) using known variant sites  
- Generated recalibrated BAM files  

---

### Alignment Quality Metrics

- Calculated alignment statistics using samtools  
- Computed:
  - Read depth  
  - Coverage statistics  
  - Alignment summary metrics  

---

### Variant Calling

- Generated intermediate GVCF using GATK HaplotypeCaller  
- Performed joint genotyping using GenotypeGVCFs  

---

### Variant Quality Control

- Created sites-only VCF  
- Applied Variant Quality Score Recalibration (VQSR) for:
  - INDELs  
  - SNPs  

- Used multiple training datasets for recalibration  

---

### Final Variant Set

- Produced filtered, high-confidence variant calls  
- Generated variant-level quality metrics  

---

## Technologies Used

- Bash (Unix scripting)  
- SLURM (HPC job scheduling)  
- FastQC  
- Cutadapt  
- BWA  
- Samtools  
- GATK  

---

## Key Features

- End-to-end NGS processing pipeline  
- Implementation of GATK Best Practices  
- Parallelized execution using HPC resources  
- Automated file handling and checkpoint recovery  
- Multi-step quality control and validation  

---

## Key Skills Demonstrated

- NGS data processing and variant calling  
- High-performance computing (SLURM workflows)  
- Alignment and post-processing of sequencing data  
- Variant quality control and recalibration  
- Pipeline automation and error handling  
- Handling large-scale genomic datasets  

---

## Notes

- Pipeline is optimized for HPC environments using SLURM  
- Intermediate files are managed to reduce storage usage  
- Variant recalibration uses established reference datasets  
- Error handling ensures recovery of intermediate outputs  

---

## Author

Lakshita Arunkumar
