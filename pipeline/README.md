Pipeline folder

This folder contains a simple bash pipeline `fetch_and_map.sh` that demonstrates fetching SRA FASTQ files, trimming, aligning (HISAT2) and counting (featureCounts).

Usage example (from repository root):

```bash
./pipeline/fetch_and_map.sh SRRxxxxxx /path/to/genome.fa /path/to/annotations.gtf /path/to/hisat2_index_prefix 8
```

Prerequisites
- SRA Toolkit (prefetch, fasterq-dump)
- fastp (optional)
- HISAT2 and samtools (or STAR)
- subread (featureCounts)

Notes
- The script is intentionally simple for clarity. For production use consider Nextflow/Snakemake and containerized environments (Docker/Singularity).
