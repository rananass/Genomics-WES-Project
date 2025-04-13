# Whole Genome Variant Analysis Pipeline

## Introduction
This workflow documents the steps followed to process whole-genome sequencing (WGS) data for a trio (father, mother, and proband). It includes quality control, trimming, alignment, variant calling, annotation, and filtering. Each step includes the rationale behind using specific tools and parameters.

## Step 1: Checking Software Versions
Before running the analysis, software versions were checked to ensure compatibility.

### Commands:
fastqc --version

trimmomatic --version

sudo apt update && sudo apt install samtools -y

samtools --version

bwa index --version

apt update && apt install -y openjdk-17-jdk

java -version

gatk --version

### Rationale:
Ensuring that all tools are installed and compatible with the workflow.

## Step 2: Quality Control (FastQC)
FastQC was used to assess the quality of raw sequencing reads.

### Commands:
cd /path/to/raw_data
ls -l *.fq
fastqc *.fq
ls -l *.html *.zip
### Rationale:
FastQC provides metrics such as base quality, GC content, adapter contamination, and sequence duplication levels.

## Step 3: Read Trimming (Trimmomatic)
Trimmomatic was used to remove low-quality reads and adapters.

### Commands:
trimmomatic PE -phred33 father_R1.fq father_R2.fq -baseout father MINLEN:5
trimmomatic PE -phred33 mother_R1.fq mother_R2.fq -baseout mother MINLEN:5
trimmomatic PE -phred33 proband_R1.fq proband_R2.fq -baseout proband MINLEN:5
ls
### Rationale:
Removing low-quality bases improves downstream alignment and variant calling. MINLEN:5 ensures that very short reads are discarded.

## Step 4: Reference Genome Preparation
Downloading and indexing the reference genome (GRCh38).

### Commands:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict
### Rationale:
Indexing the genome allows efficient alignment and variant calling.

## Step 5: Read Alignment (BWA-MEM)
Reads were aligned to the reference genome using BWA-MEM.

### Commands:
bwa mem -t 4 GRCh38.primary_assembly.genome.fa father_R1.fq father_R2.fq > father_aligned.sam
bwa mem -t 4 GRCh38.primary_assembly.genome.fa mother_R1.fq mother_R2.fq > mother_aligned.sam
bwa mem -t 4 GRCh38.primary_assembly.genome.fa proband_R1.fq proband_R2.fq > proband_aligned.sam
### Rationale:
BWA-MEM efficiently maps reads to the reference genome, allowing downstream analysis.

## Step 6: Convert SAM to BAM (Samtools View)
SAM files were converted to BAM format for better storage and processing efficiency.

### Commands:
samtools view -S -b father_aligned.sam > father_aligned.bam
samtools view -S -b mother_aligned.sam > mother_aligned.bam
samtools view -S -b proband_aligned.sam > proband_aligned.bam
### Rationale:
BAM files are more compact and efficient for variant calling.

## Step 7: Sorting and Indexing BAM Files
Sorting ensures alignments are ordered by genomic coordinates.

### Commands:
samtools sort father_aligned.bam -o father_sorted.bam
samtools sort mother_aligned.bam -o mother_sorted.bam
samtools sort proband_aligned.bam -o proband_sorted.bam
samtools index father_sorted.bam
samtools index mother_sorted.bam
samtools index proband_sorted.bam
bwa index GRCh38.primary_assembly.genome.fa
ls -lh GRCh38.primary_assembly.genome.fa.*
### Rationale:
Sorted BAM files are required for downstream variant calling and analysis.

## Step 8: Variant Calling (GATK HaplotypeCaller)
GATK HaplotypeCaller was used to detect variants.
