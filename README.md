# Bulk_RNA_Seq_Analysis

Reproducible repository for bulk RNA-seq analyses (colon organoids and keratinocyte cell lines) – from my undergrad 2nd-3rd year (2022-2023)

## Project summary
This repository contains scripts, metadata, raw and processed count tables, and results for differential expression analyses and downstream visualizations. The recent reorganization added a recommended folder layout, helper files for reproducible environments, and a `.gitignore` tuned to ignore intermediate data files.

## Repository Structure

```
.
├── colon organoids/      # Analysis scripts and data for colon organoid experiments
├── keratinocyte/         # Analysis scripts and data for keratinocyte cell lines
├── pipeline/             # Reusable pipeline scripts for RNA-seq analysis
├── data/                 # Raw and processed data files
├── results/              # Output files from analyses
├── scripts/             # Utility scripts and shared functions
├── requirements.txt     # Python package dependencies
└── R_requirements.R    # R package dependencies
```

## From SRA / GEO FASTQ to gene-level counts (example)

The following shell snippet shows a common end-to-end workflow: fetch FASTQ from SRA/GEO, optionally trim, align to a reference, and generate gene-level counts with `featureCounts`. There is an alternative Salmon (quasi-mapping) example included.


### Prerequisites
- **Bash tools for FASTQ to counts:**
	- SRA Toolkit (`prefetch`, `fasterq-dump`): download FASTQ from SRA/GEO
	- FASTQC: quality control of raw FASTQ files
	- fastp: trimming and filtering reads
	- HISAT2 or STAR: alignment to reference genome
	- samtools: BAM file manipulation
	- subread (featureCounts): gene-level counting
	- Salmon: fast transcript quantification (alignment-free)

**Install commands (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install fastqc hisat2 star samtools
# SRA Toolkit
conda install -c bioconda sra-tools
# fastp
conda install -c bioconda fastp
# subread (featureCounts)
conda install -c bioconda subread
# Salmon
conda install -c bioconda salmon
```
Or use [mamba](https://github.com/mamba-org/mamba) for faster conda installs.

**Role in pipeline:**
- FASTQC: run after download/trimming for QC reports
- fastp: trim/filter before alignment
- HISAT2/STAR: align reads to genome
- samtools: sort/index BAM files
- featureCounts: count reads per gene
- Salmon: alternative quantification (transcript/gene)

### Example: end-to-end (HISAT2 + featureCounts)

```bash
# variables
SRR=SRR22450503   # replace with SRA run accession 
THREADS=8
REF_FA=/path/to/genome.fa
GTF=/path/to/annotations.gtf
HISAT2_INDEX=/path/to/hisat2_index_prefix
OUTDIR=data/raw/${SRR}
mkdir -p ${OUTDIR}

# 1) fetch FASTQ
prefetch ${SRR}
fasterq-dump --split-files --outdir ${OUTDIR} ${SRR}

# 2) trim (optional)
fastp -i ${OUTDIR}/${SRR}_1.fastq -I ${OUTDIR}/${SRR}_2.fastq \
			-o ${OUTDIR}/${SRR}_1.trim.fastq -O ${OUTDIR}/${SRR}_2.trim.fastq \
			-w ${THREADS} -h ${OUTDIR}/${SRR}.fastp.html -j ${OUTDIR}/${SRR}.fastp.json

# 3) align with HISAT2
hisat2 -p ${THREADS} -x ${HISAT2_INDEX} \
	-1 ${OUTDIR}/${SRR}_1.trim.fastq -2 ${OUTDIR}/${SRR}_2.trim.fastq \
	| samtools sort -@ ${THREADS} -o ${OUTDIR}/${SRR}.sorted.bam
samtools index ${OUTDIR}/${SRR}.sorted.bam

# 4) count reads at gene level
featureCounts -T ${THREADS} -a ${GTF} -o ${OUTDIR}/${SRR}.featureCounts.txt ${OUTDIR}/${SRR}.sorted.bam

# The counts file contains gene-level raw counts suitable for DESeq2/edgeR downstream.
```

### Quick alternative: Salmon (fast, alignment-free)

```bash
# create salmon index (once)
salmon index -t ${REF_FA} -i salmon_index --type quasi -k 31

# quantify
salmon quant -i salmon_index -l A \
	-1 ${OUTDIR}/${SRR}_1.trim.fastq -2 ${OUTDIR}/${SRR}_2.trim.fastq \
	-p ${THREADS} -o ${OUTDIR}/salmon_${SRR}

# export gene-level counts (tximport in R or use salmon2gene tools)
```

### Notes
- Replace variables (paths, indexes) with your environment values.
- For GEO datasets that provide raw FASTQ links, use `wget` or `curl` instead of `prefetch`.
- If your reads are single-end adjust `fasterq-dump` and alignment commands accordingly.
- Save these pipeline scripts in `pipeline/` (created below) and run from project root so outputs go to `data/raw/` or `data/processed/`.
## Contents
- `colon organoids/`, `keratinocyte/` — existing analysis folders with R scripts, count matrices and results.
- `data/` — recommended place for data; contains `raw/` and `processed/` subfolders.
- `scripts/` — recommended place for analysis and utility scripts.
- `results/` — figures, tables and exported outputs.
- `R_requirements.R` — R and Bioconductor package installer to reproduce the R environment.
- `requirements.txt` — Python dependencies for any auxiliary scripts.

## Getting started

Prerequisites
- R (>= 4.0) and Rtools on Windows if compiling packages.
- Optional: Python 3.8+ if you use Python scripts.

Quick start
1. Clone the repository and change into it:

```powershell
git clone <repo-url>
cd "RNA Seq Analysis"
```

2. Install R packages (run inside an R session):

```r
source("R_requirements.R")
```

3. (Optional) Create and activate a Python virtual environment and install Python dependencies:

```powershell
python -m venv .venv; .\.venv\Scripts\Activate.ps1; pip install -r requirements.txt
```

4. Place raw immutable inputs (FASTQ, original count matrices) in `data/raw/`. Place derived outputs in `data/processed/`.

## Recommended repository layout
- `data/raw/` — immutable original inputs.
- `data/processed/` — derived tables (normalized counts, filtered results).
- `scripts/` — modular, documented scripts; place helper functions in `scripts/utils/` as needed.
- `results/` — final figures, tables and CSVs used for reporting.

## Notes on `.gitignore` and tracking CSVs
The repository `.gitignore` was updated to ignore `*.csv` and common R artifacts to avoid committing large or intermediate tables. Important: files already tracked by git will remain tracked. To stop tracking a file already in the repo, run:

```powershell
git rm --cached path\\to\\file.csv
git commit -m "Stop tracking large CSV"
```

If you want specific CSVs to remain tracked, I can add explicit negation rules (e.g., `!colon organoids/meta_colon.csv`) to the `.gitignore`.

## How to run analyses
- Inspect the R scripts in `colon organoids/` and `keratinocyte/` to identify pipeline steps.
- Prefer running analyses from the repository root so that relative paths (e.g., `data/processed/`) resolve correctly.
- Example: run an R script from PowerShell:

```powershell
Rscript scripts/my_analysis.R
```

## Reproducibility recommendations
- Record package versions using `sessionInfo()` in R and include this output in result folders.
- Use `data/raw` vs `data/processed/` split to prevent accidental modification of raw inputs.
- Use branches and PRs for large reorganizations.

## Contributing
- Open an issue to propose large structural changes.
- Use branches for work and submit PRs with small, reviewable changes.

## License & contact
This project uses the license in `LICENSE`. For questions, open an issue or contact me at [joy21.dev.pd@gmail.com](mailto:joy21.dev.pd@gmail.com).

---

