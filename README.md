# hgt-detection
It is a software takes use of existing shotgun NGS datasets to detect HGT breakpoints, identify the transferred genome segments, and reconstructs the inserted local haplotype.
## Table of Contents
1. [Installation](#installation)
2. [HGT-detection usage](#hgt-detection-usage)
3. [Example workflow](#example-workflow)
## Installation
### Requirements
- Softwares
  - BWA(0.7.12+)
  - Samtools(1.3.1+)
  - Bedtools(2.20.1+)
  - LUMPY
  - Python(3.5+)
- Python packages(3.5+)
  - pysam
  - sklearn
  - parasail
  - numpy
  - scipy
  - lmfit
### Install
git clone --recursive https://github.com/lichen2018/hgt-detection.git
## HGT-detection usage
### 1. Detect raw HGT breakpoints.
```
usage: python hgt-detection/scripts/get_raw_bkp.py [options]
```
#### Required arguments  
```
  -r FILE  Metagenomic Reference 
  -u FILE  unique reads bam file
  -o FILE  raw breakpoints file
  ```
#### Option arguments
  ```
  -t INT  number of threads [4]
  ```
### 2. Detect accurate HGT breakpoints.
```
usage: python hgt-detection/scripts/get_accurate_bkp.py [options]
```
#### Required arguments
  ```
  -r        FILE  Metagenomic Reference
  -s        FILE  split reads bam file
  --raw_bkp FILE  raw breakpoints file
  -o        FILE  accurate reakpoints file
  ```
#### Option arguments
  ```
  -t INT  number of threads [10]
  ```
### 3. Get HGT references.
```
usage: python hgt-detection/scripts/get_reference.py [options]
```
#### Required arguments
  ```
  -r        FILE  Metagenomic Reference
  -c        FILE  coverage file
  -id       STR   Sample name or id
  --acc_bkp FILE  accurate reakpoints file
  --out_dir STR   path to the directory where results should be stored
  ```
## Example workflow
### Preprocessing
```
# Align the data
bwa mem -M -t 8 -R "@RG\tID:id\tSM:sample\tLB:lib" Metagenomic_reference.fasta sample.1.fq sample.2.fq \
  | samtools view -bhS -> sample.unsort.bam

# Sort bam file
samtools sort -o sample.bam sample.unsort.bam

# Extract split reads
samtools view -h sample.bam \
  | lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
  | samtools view -Sb > sample.unsort.splitters.bam

# Sort split reads bam file
samtools sort -o sample.splitters.bam sample.unsort.splitters.bam

# Extract unique reads bam file
samtools view -q 20 -b sample.bam > sample.unique.bam

# Calculate coverage
bedtools genomecov -ibam sample.bam -bg > sample.coverage.txt
```
### Running HGT-detection
```
# 1. Detect raw HGT breakpoints.
python hgt-detection/scripts/get_raw_bkp.py -r meta_ref.fasta -u sample.unique.bam -o sample.raw.txt

# 2. Detect accurate HGT breakpoints.
python hgt-detection/scripts/get_accurate_bkp.py -r meta_ref.fasta -s sample.splitters.bam --raw_bkp sample.raw.txt -o sample.acc.txt

# 3. Get HGT references.
python hgt-detection/scripts/get_reference.py -r meta_ref.fasta -c sample.coverage.txt -id sample --acc_bkp sample.acc.txt --out_dir ./
```
