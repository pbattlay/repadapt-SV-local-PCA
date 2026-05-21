# RepAdapt SV Local PCA

## Requirements
These scripts are designed for HPC clusters running SLURM. Update the `module load` 
commands in `script1.sh` and `script2.sh` to match your cluster's module names.

### Cluster modules
- R (4.4.0 or later, ideally built with MKL for vectorised LD computation)
- plink2
- bcftools
- java

### R packages
```r
install.packages(c("data.table", "BEDMatrix", "ggplot2", "diptest"))

# lostruct is not on CRAN — install from GitHub:
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("petrelharp/local_pca/lostruct")
```

### BEAGLE
Download `beagle.jar` from the [BEAGLE website](https://faculty.washington.edu/browning/beagle/beagle.html) 
and place it in the `scripts/` directory, or update `BEAGLE_JAR` in your config file to point to its location.

## Setup
Set the base directory
```bash
base_dir="/home/pbattlay/om62_scratch2/rasvlpca_auto_2"
```

Set the repadapt ID or name for the species
```bash
raid=rawg0053
```

Ensure the config file is in the base directory and named correctly. Update the config file (normally `$raid` and `$vcf` will change per species)
```bash
source config_${raid}.sh
```

The `scripts/` directory should contain `script1.sh`, `script2.sh`, `script3.sh` and `beagle.jar`
```bash
ls ${base_dir}/scripts
```

Set up directories in the species directory
```bash
mkdir "${base_dir}/${raid}"

cd "${base_dir}/${raid}"

mkdir -p errout ref vcf \
         mds/plots mds/candidates \
         local_pca/raw local_pca/tables local_pca/plots local_pca/axis_scores \
         local_pca/collapsed \
         breakpoints/plots \
         final/plots
```

Make a contig; size table from the VCF header
```bash
module load bcftools

bcftools view -h "${vcf}" | awk '
BEGIN { OFS="\t" }
/^##contig=<.*ID=/ {
  id = len = ""
  if (match($0, /ID=([^,>]+)/, a)) id = a[1]
  if (match($0, /length=([0-9]+)/, b)) len = b[1]
  if (id != "" && len != "") print id, len
}' > ${base_dir}/${raid}/ref/contigs.tsv
```

Not all VCF headers will contain this information (but RepAdapt VCFs will). If you can't get contig length information from the VCF header, use a FASTA index
```bash
fai=""
awk '{print $1"\t"$2}' $fai > ${base_dir}/${raid}/ref/contigs.tsv
```

Make a list of 'chromosomes' (scaffolds >= `min_chr_length`)
```bash
awk -v min="${min_chr_length}" '$2 >= min { print $1 }' \
    ${base_dir}/${raid}/ref/contigs.tsv > ${base_dir}/${raid}/ref/chrs.txt
```

## Script 1 — filter and impute
Runs as a SLURM array (one task per chromosome). For each chromosome, the script filters the input VCF (designed for RepAdapt 'raw' VCFs) for genotype quality and depth, missingness, MAF, and biallelic SNPs, imputes missing data with BEAGLE, and outputs a PLINK BED fileset. Intermediate VCFs are removed on completion. Note that if a VCF doesn't contain the correct INFO fields, filtering will remove all varaints!
```bash
cd ${base_dir}

sbatch --account="${account}" \
    --job-name="s1_${raid}" \
    --array="1-$(wc -l < ${base_dir}/${raid}/ref/chrs.txt)" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script1.sh" "${base_dir}/${base_dir}/config_${raid}.sh"
```

## Script 1b - diagnostic plots
Runs a single SLURM job, produces diagnostic plots only to `${base_dir}/vcf/plots` (not required for script2)
- `scaffold_sizes.pdf` histogram of reference scaffold size distribution with `min_chr_length` filter marked
- `snps_per_window.pdf` histogram of SNP counts per window distribution with `MIN_SNPS_PER_WINDOW` marked
- `pca.pdf` scatterplot of PC1 and PC2 (100k random SNPs post-filtering)
```bash
cd ${base_dir}

sbatch --account="${account}" \
    --job-name="s1b_${raid}" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script1b.sh" "${base_dir}/config_${raid}.sh"
```

## Script 2 - MDS scan and local PCA
Runs as a SLURM array (one task per chromosome). Each task runs a single R session with four sequential sections.
- Section A performs windowed local PCA across the chromosome, computes MDS on pairwise window distances, and identifies candidate inversion regions as MAD Z-score outliers on each MDS axis using a permutation test (adapted from Huang et al. 2020 Mol Ecol).
- Section B validates each candidate with local PCA clustering: samples are assigned to three genotype classes (LL/LR/RR) and filtered by cluster separation, dip test, R² clustering quality, and heterozygote excess.
- Section C collapses overlapping candidates with correlated genotypes into representative inversions.
- Section D computes chromosome-wide LD with each candidate's inversion genotype vector and scans the binned LD profile to refine breakpoint coordinates.
- Section D2 performs a second collapse round using LD-refined coordinates, writes the final per-chromosome breakpoint and genotype output files, and produces diagnostic plots for each candidate.
```bash
sbatch --account="${account}" \
    --job-name="s2_${raid}" \
    --array="1-$(wc -l < ${base_dir}/${raid}/ref/chrs.txt)" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script2.sh" "${base_dir}/config_${raid}.sh"
```

## Script 3 - merge results
Runs a single SLURM job to merge outputs across chromosomes, collect plots for surviving SV candidates and asses genome-wide LD between candidates.

```bash
sbatch --account="${account}" \
    --job-name="s3_${raid}" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script3.sh" "${base_dir}/config_${raid}.sh"
```

### WARNING! This erases results from scripts 2 and 3. Useful for testing before a rerun
```bash
rm -f ${base_dir}/${raid}/mds/*
rm -f ${base_dir}/${raid}/mds/*/*
rm -f ${base_dir}/${raid}/local_pca/*
rm -f ${base_dir}/${raid}/local_pca/*/*
rm -f ${base_dir}/${raid}/breakpoints/*
rm -f ${base_dir}/${raid}/breakpoints/*/*
rm -f ${base_dir}/${raid}/final/*
rm -f ${base_dir}/${raid}/final/*/*
```
