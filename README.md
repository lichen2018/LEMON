
# hgt-detection
It is a software takes use of existing shotgun NGS datasets to detect HGT breakpoints, identify the transferred genome segments, and reconstructs the inserted local haplotype.
## Table of Contents
1. [Installation](#readme)
2. [hgt-detection usage](#readme)
3. [Example workflow](#readme)
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
git clone --recursive https://github.com/lichen2018/hgt-detection.git
