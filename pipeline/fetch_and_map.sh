#!/usr/bin/env bash
# Simple pipeline to fetch SRA runs (via prefetch/fasterq-dump), trim, align (HISAT2) and count (featureCounts)
# Usage: ./fetch_and_map.sh SRR22450503 /path/to/genome.fa /path/to/annotations.gtf /path/to/hisat2_index_prefix

#!/usr/bin/env bash
# Real analysis pipeline: batch SRA/FASTQ processing, auto-detect paired/single-end, log summary

# Usage (real analysis):
#   ./pipeline/fetch_and_map.sh pipeline/manifest.tsv /path/genome.fa /path/annotations.gtf /path/hisat2_index_prefix [threads]
#
#   manifest.tsv (tab-separated, generated from SRR_Acc_List.txt):
#   sample_id	type	paired
#   SRR22450503	SRA	1
#   SRR22450504	SRA	1
#   ...
#
# Columns:
#   sample_id: SRA run accession or FASTQ basename
#   type: SRA (download) or FASTQ (local)
#   paired: 1 for paired-end, 0 for single-end

set -euo pipefail

MANIFEST=${1:-}
REF_FA=${2:-}
GTF=${3:-}
HISAT2_INDEX=${4:-}
THREADS=${5:-8}

if [ -z "${MANIFEST}" ] || [ -z "${REF_FA}" ] || [ -z "${GTF}" ] || [ -z "${HISAT2_INDEX}" ]; then
  echo "Usage: $0 manifest.tsv /path/genome.fa /path/annotations.gtf /path/hisat2_index_prefix [threads]"
  exit 1
fi

if [ ! -f "${MANIFEST}" ]; then
  echo "Manifest file not found: ${MANIFEST}"
  exit 1
fi

LOG=results/pipeline_summary_$(date +%Y%m%d_%H%M%S).log
mkdir -p results/
echo "Sample\tType\tPaired\tCountsFile\tTotalReads" > "$LOG"

while IFS=$'\t' read -r SAMPLE TYPE PAIRED; do
  OUTDIR=data/raw/${SAMPLE}
  mkdir -p "${OUTDIR}"
  if [ "$TYPE" = "SRA" ]; then
    prefetch "$SAMPLE"
    fasterq-dump --split-files --outdir "$OUTDIR" "$SAMPLE"
    if [ "$PAIRED" = "1" ]; then
      FQ1="$OUTDIR/${SAMPLE}_1.fastq"
      FQ2="$OUTDIR/${SAMPLE}_2.fastq"
    else
      FQ1="$OUTDIR/${SAMPLE}.fastq"
      FQ2=""
    fi
  else
    # FASTQ mode: expects files in OUTDIR
    FQ1="$OUTDIR/${SAMPLE}_1.fastq"
    FQ2="$OUTDIR/${SAMPLE}_2.fastq"
    if [ ! -f "$FQ1" ]; then FQ1="$OUTDIR/${SAMPLE}.fastq"; FQ2=""; fi
  fi

  # Trim (optional)
  if command -v fastp >/dev/null 2>&1; then
    if [ -n "$FQ2" ]; then
      fastp -i "$FQ1" -I "$FQ2" -o "$OUTDIR/${SAMPLE}_1.trim.fastq" -O "$OUTDIR/${SAMPLE}_2.trim.fastq" -w "$THREADS" -h "$OUTDIR/${SAMPLE}.fastp.html" -j "$OUTDIR/${SAMPLE}.fastp.json"
      FQ1="$OUTDIR/${SAMPLE}_1.trim.fastq"
      FQ2="$OUTDIR/${SAMPLE}_2.trim.fastq"
    else
      fastp -i "$FQ1" -o "$OUTDIR/${SAMPLE}.trim.fastq" -w "$THREADS" -h "$OUTDIR/${SAMPLE}.fastp.html" -j "$OUTDIR/${SAMPLE}.fastp.json"
      FQ1="$OUTDIR/${SAMPLE}.trim.fastq"
      FQ2=""
    fi
  fi

  # Align
  if command -v hisat2 >/dev/null 2>&1; then
    if [ -n "$FQ2" ]; then
      hisat2 -p "$THREADS" -x "$HISAT2_INDEX" -1 "$FQ1" -2 "$FQ2" | samtools sort -@ "$THREADS" -o "$OUTDIR/${SAMPLE}.sorted.bam"
    else
      hisat2 -p "$THREADS" -x "$HISAT2_INDEX" -U "$FQ1" | samtools sort -@ "$THREADS" -o "$OUTDIR/${SAMPLE}.sorted.bam"
    fi
    samtools index "$OUTDIR/${SAMPLE}.sorted.bam"
  else
    echo "hisat2 not found. Please install HISAT2 or adjust script to use STAR/Salmon."
    continue
  fi

  # Count
  if command -v featureCounts >/dev/null 2>&1; then
    featureCounts -T "$THREADS" -a "$GTF" -o "$OUTDIR/${SAMPLE}.featureCounts.txt" "$OUTDIR/${SAMPLE}.sorted.bam"
    COUNTS="$OUTDIR/${SAMPLE}.featureCounts.txt"
    TOTAL=$(samtools view -c "$OUTDIR/${SAMPLE}.sorted.bam")
    echo -e "$SAMPLE\t$TYPE\t$PAIRED\t$COUNTS\t$TOTAL" >> "$LOG"
  else
    echo "featureCounts not found. Install subread or use an alternative counting method."
    continue
  fi
done < <(tail -n +2 "$MANIFEST")
