# RepAdapt SV Local PCA

## Requirements
These scripts are designed for HPC clusters running SLURM. Update the `module load` 
commands in the `script1*.sh` files and `script2.sh` to match your cluster's module names.

### Cluster modules
- R (4.4.0 or later, ideally built with MKL for vectorised LD computation)
- plink2
- bcftools
- bedtools (used by script0 for the mappability mask)
- java

### genmap (mappability)
genmap is **not** a cluster module — it lives in a conda env. Point `GENMAP_BIN` at the top of `script0.sh` (e.g. `/home/<user>/.conda/envs/genmap/bin/genmap`); `GENMAP_K`, `GENMAP_E`, and `MAPPABILITY_THRESHOLD` are set there too. Recreate the env if it's missing:
```bash
conda create -n genmap -c bioconda genmap
```

### Reference FASTAs
Mappability is decoupled from the per-species pipeline. Collect the reference FASTA(s) the VCFs were mapped to into one directory (e.g. `${base_dir}/FASTA/`) and run script0 over it (see Script 0); it writes a `<prefix>.mappable.bed` next to each FASTA. The mask is a property of the reference, so species sharing a reference share one `.bed`. The per-species config no longer holds `FASTA` — it just points `MAPPABILITY_MASK` at the relevant `.bed`.

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
base_dir="/home/pbattlay/om62_scratch2/REPADAPT_LPCA"
```

Set the repadapt ID or name for the species
```bash
raid=rawg0001
```

Ensure the config file is in the base directory and named correctly. Update the config file (normally `$raid` and `$vcf` will change per species)
```bash
source config_${raid}.sh
```

The `scripts/` directory should contain `script0.sh`, `script1a.sh`, `script1b.sh`, `script1c.sh`, `script2.sh`, `script3.sh`, `nuke_script2_and_3_results.sh` and `beagle.jar`
```bash
ls ${base_dir}/scripts
```

Set up directories in the species directory
```bash
mkdir -p ${base_dir}/${raid}/errout \
         ${base_dir}/${raid}/ref \
         ${base_dir}/${raid}/vcf/plots \
         ${base_dir}/${raid}/mds/plots \
         ${base_dir}/${raid}/mds/candidates \
         ${base_dir}/${raid}/local_pca/raw \
         ${base_dir}/${raid}/local_pca/tables \
         ${base_dir}/${raid}/local_pca/plots \
         ${base_dir}/${raid}/local_pca/axis_scores \
         ${base_dir}/${raid}/local_pca/collapsed \
         ${base_dir}/${raid}/breakpoints/plots \
         ${base_dir}/${raid}/final/plots \
         ${base_dir}/${raid}/inputs
```

### Chromosome names
PLINK uses the VCF chromosome names directly, and script1a/1c write per-chromosome files named after them. Some RepAdapt VCFs carry names containing a pipe, e.g. `lcl|Bpe_Chr1`, which breaks `--chr` handling and file naming. If so, rename them in the VCF *before* building the contig table below, then point `vcf=` in your config at the renamed file:
```bash
module load bcftools
zcat "${vcf}" | sed 's/lcl|/lcl_/g' | bgzip > "${vcf%.vcf.gz}_renamed.vcf.gz"
tabix -p vcf "${vcf%.vcf.gz}_renamed.vcf.gz"
```
Adjust the `sed` pattern to your prefix. If your chromosome names are already pipe-free, skip this.

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
#fai=""
#awk '{print $1"\t"$2}' $fai > ${base_dir}/${raid}/ref/contigs.tsv
```

Make a list of 'chromosomes' (scaffolds >= `min_chr_length`)
```bash
awk -v min="${min_chr_length}" '$2 >= min { print $1 }' \
    ${base_dir}/${raid}/ref/contigs.tsv > ${base_dir}/${raid}/ref/chrs.txt
```

### Optional inputs
Four optional files live in `${base_dir}/${raid}/inputs/`, and the config points at them by default. **Presence = active**: each script ignores any file it can't find, so to disable one you remove or empty the file (not the config variable).

| file | format | used by |
|------|--------|---------|
| `mappability_mask.bed` | 0-based BED keep-list (from script0) | script1a |
| `remove.tsv` | sample IIDs, one per line, no header | script1b |
| `truth.tsv` | `chr  start  end`, 1-based, no header | script3 overview plot |
| `manual_candidates.tsv` | `chr  start  end`, 1-based, no header | script2 |

`mappability_mask.bed` is the only 0-based input (genmap/bedtools/plink convention); `truth.tsv` and `manual_candidates.tsv` are 1-based inclusive, matching the pipeline's candidate coordinates. Regions in `manual_candidates.tsv` are injected into the script2 candidate set and tested for genotype clustering exactly like MDS candidates — this bypasses MDS *detection*, not the cluster/het validation, so a region that doesn't cluster is still dropped. Chromosome names in all four files must match the (possibly renamed) VCF names.

## Script 0 — mappability mask
Decoupled from the per-species pipeline: runs **once per reference**, as a SLURM array over a directory of FASTAs. Each task builds a uniquely-mappable mask from one FASTA with genmap (`index` → `map -K${GENMAP_K} -E${GENMAP_E} -bg`), keeps per-base positions clearing `MAPPABILITY_THRESHOLD` (1 = unique), and writes a genome-wide `<prefix>.mappable.bed` next to the FASTA. genmap parameters and `GENMAP_BIN` are set at the top of `script0.sh`; the FASTA directory is the only argument. The job is guarded — reruns skip the genmap index/map if the bedgraph is cached.

script1a restricts SNPs to this mask **before** missingness/`--mind`/MAF and before imputation, which removes collapsed-paralog SNPs — high-`QUAL`/`GQ`/`DP`, balanced-heterozygote artefacts that the quality filters do *not* catch. Masking is **optional**: if `MAPPABILITY_MASK` is absent, script1a still runs but warns loudly, and you keep those artefacts. Point each species config's `MAPPABILITY_MASK` at the relevant `.bed` (in `inputs/`, or directly at the file next to the FASTA).

Submit one array task per FASTA — the count comes from script0 itself, so it can't drift from the manifest the worker indexes:
```bash
FASTA_DIR="${base_dir}/FASTA"
N=$(bash "${base_dir}/scripts/script0.sh" --count "${FASTA_DIR}")
sbatch --account="${account}" --job-name="s0" \
    --array="1-${N}" "${base_dir}/scripts/script0.sh" "${FASTA_DIR}"
```

## Script 1a — filter and diagnostic plots
Runs as a single SLURM job. Filters the input VCF (designed for RepAdapt 'raw' VCFs) across all analysis scaffolds in one pass — first restricting SNPs to the mappability mask (`MAPPABILITY_MASK`, from script0; optional — script1a warns loudly if absent), then genotype quality and depth, per-site missingness (`--geno`), per-individual missingness (`--mind`), MAF, and biallelic SNPs — and writes a single PLINK BED fileset (`vcf/${raid}.filtered`). If not using a verified RepAdapt VCF, check that your VCF contains the correct INFO fields. It then writes diagnostic plots to `${base_dir}/${raid}/vcf/plots`:
- `scaffold_sizes.pdf` histogram of reference scaffold size distribution with `min_chr_length` marked
- `window_sizes.pdf` histogram of the physical span (kb) of each `SNP_WINDOW_SIZE`-SNP window — a decision aid for the scan window size
- `pca_preimpute.pdf` scatterplot of PC1 vs PC2 (100k random SNPs, post-filtering / pre-imputation)
- `pca_preimpute.tsv` table of `sample`, `PC1`, `PC2` — use this to identify outlier samples to exclude

Sample exclusion is deliberately *not* applied here, so per-site missingness is computed on the full sample set.
```bash
sbatch --account="${account}" \
    --job-name="s1a_${raid}" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script1a.sh" "${base_dir}/config_${raid}.sh"
```

## Script 1b — sample exclusion
Runs as a single SLURM job. After reviewing the PCA, list the sample IIDs to drop (one per line — the `sample` column of `pca_preimpute.tsv`) in `inputs/remove.tsv` (the default `REMOVE_SAMPLES` path). This step drops those individuals from the filtered BED and writes `vcf/${raid}.kept`; it removes individuals only and does not recompute site missingness. If the file is absent or empty it is a no-op and script1c reads `vcf/${raid}.filtered` directly. It can be re-run with an updated list without re-running script1a.
```bash
sbatch --account="${account}" \
    --job-name="s1b_${raid}" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script1b.sh" "${base_dir}/config_${raid}.sh"
```

## Script 1c — MAF refilter and impute
Runs as a SLURM array (one task per chromosome). For each chromosome it reads the analysis BED (`vcf/${raid}.kept` if samples were excluded, otherwise `vcf/${raid}.filtered`), re-applies the MAF filter on the retained samples (allele frequencies change after exclusion; `--geno` is *not* re-applied), exports a per-chromosome VCF, imputes missing data with BEAGLE, and writes a PLINK BED fileset (`vcf/${chr}`). Intermediate VCFs are removed on completion.
```bash
sbatch --account="${account}" \
    --job-name="s1c_${raid}" \
    --array="1-$(wc -l < ${base_dir}/${raid}/ref/chrs.txt)" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script1c.sh" "${base_dir}/config_${raid}.sh"
```

## Script 2 - MDS scan and local PCA
Runs as a SLURM array (one task per chromosome). Each task runs a single R session with four sequential sections.
- Section A performs windowed local PCA across the chromosome, computes MDS on pairwise window distances, and identifies candidate inversion regions as MAD Z-score outliers on each MDS axis using a permutation test (adapted from Huang et al. 2020 Mol Ecol).
- Section B validates each candidate with local PCA clustering: samples are assigned to three genotype classes (LL/LR/RR) and filtered by cluster separation, dip test, R² clustering quality, and heterozygote excess.
- Section C performs the first collapse: candidates that spatially **overlap** and whose inversion genotypes are correlated (`|r| >= R_THRESH`, default 0.8) are merged, keeping the **widest** candidate (by physical span) as the representative — MDS-preferred on ties — so a manual candidate's full extent wins over a narrower correlated MDS fragment. The required overlap is its own safety, so a moderate correlation gate is enough.
- Section D computes chromosome-wide LD with each candidate's inversion genotype vector and scans the binned LD profile to refine breakpoint coordinates.
- Section D2 performs the **long-range collapse**: candidates need not overlap — any pair within a physical edge-to-edge gap of `COLLAPSE_BUFFER` bp (default `Inf`, i.e. chromosome-wide) and correlated at `|r| >= R_THRESH_LONG` (default 0.9) is merged. Because there is no overlap to lean on, `R_THRESH_LONG` is stricter than `R_THRESH`. Each merged union is then **re-clustered for validation**: if it re-clusters as one clean three-class inversion the merge is kept (whole-region genotype + stats, `genotype_source = merged`); if it does not, the merge is **reverted** to its constituent fragments with their original coordinates and genotypes restored (`genotype_source = fragment_reverted`). Singletons whose LD-refined region fails to re-cluster keep their fragment genotype (`genotype_source = fragment`). Section D2 then writes the final per-chromosome breakpoint and genotype files and a diagnostic plot per candidate.

Because the long-range gate has no overlap safety, `R_THRESH_LONG` is a per-species dial: 0.9 suits sunflower (within-inversion `r >= 0.95` separates cleanly from distinct/supergene `r <= 0.76`), but species with elevated background LD need it higher — e.g. capeweed has within-scaffold distinct pairs reaching `r ~ 0.91`, so set its `R_THRESH_LONG` to ~0.95 (or hold the long-range pass off that species) and lean on the re-cluster guard.
```bash
sbatch --account="${account}" \
    --job-name="s2_${raid}" \
    --array="1-$(wc -l < ${base_dir}/${raid}/ref/chrs.txt)" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script2.sh" "${base_dir}/config_${raid}.sh"
```

## Script 3 - merge results
Runs a single SLURM job to merge outputs across chromosomes, collect plots for surviving SV candidates, assess genome-wide LD between candidates, and draw a genome overview — one chromosome per row with each candidate's merged (green) / LD (red) / MDS (blue) spans, plus an optional gold truth track from `TRUTH_SET`.

```bash
sbatch --account="${account}" \
    --job-name="s3_${raid}" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script3.sh" "${base_dir}/config_${raid}.sh"
```

### Purging a previous run before a rerun
Files from different runs carry run-specific names (`chr:start-end`), so they accumulate rather than overwrite. `nuke_script2_and_3_results.sh` erases all script 2 and script 3 outputs (`mds/`, `local_pca/`, `breakpoints/`, `final/`) for a species so the next run starts clean.

**Run it once, before launching the script 2 array — never from inside script 2.** Script 2 is a SLURM array (one task per chromosome, running in parallel); a purge inside it would have each task delete the others' output mid-run. Chain it as a dependency so the array cannot start until the purge succeeds:
```bash
clean=$(sbatch --parsable \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/nuke_script2_and_3_results.sh" "${base_dir}/config_${raid}.sh")

sbatch --account="${account}" \
    --job-name="s2_${raid}" \
    --dependency=afterok:${clean} \
    --array="1-$(wc -l < ${base_dir}/${raid}/ref/chrs.txt)" \
    -D "${base_dir}/${raid}" \
    "${base_dir}/scripts/script2.sh" "${base_dir}/config_${raid}.sh"
```
Or run the purge by hand before the array — same effect. It deletes files only and leaves the directory tree intact, and refuses to run if `base_dir`/`raid` are unset.

## Results
A final list of candidate regions and their genotypes are written to `${base_dir}/${raid}/final/${raid}_merged.tsv` and `${base_dir}/${raid}/final/${raid}_merged_genotypes.tsv` respectively. Between SV LD is written and plotted to `${base_dir}/${raid}/final/${raid}_LD_heatmap.*`, and a genome overview of all candidates (with the optional truth track) to `${base_dir}/${raid}/final/${raid}_overview.pdf`.

# MDS plots
MDS plots for final candidates can be found in `${base_dir}/${raid}/final/plots/`. Blue bars show MDS scan candidate region, red bars show LD-refined regions, green bars (if present) show region merged across multiple candidate regions. MDS scan plots for all candidates before filtering can be found in `${base_dir}/${raid}/mds/plots/` and `${base_dir}/${raid}/breakpoints/plots/`
<img width="1500" height="420" alt="h1s7_h1s7_60600001_67920000_bp" src="https://github.com/user-attachments/assets/2c5bfa2d-ee47-49e5-b0e6-cf30320fa4d1" />

# Local PCA plots
Local PCA plots for final candidates can be found in `${base_dir}/${raid}/final/plots/`; initial plots for all MDS candidates can be found in `${base_dir}/${raid}/local_pca/plots/`
<img width="1300" height="1040" alt="h1s7_h1s7_60600001_67920000_lp" src="https://github.com/user-attachments/assets/703b4626-0726-4eec-881d-d45ffcf7b80d" />
