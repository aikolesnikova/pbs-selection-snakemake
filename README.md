# Population Branch Statistic (PBS) analysis workflow

*A reproducible Snakemake pipeline for population-specific genetic differentiation, LD-based candidate prioritisation, and functional annotation.*

The **Population Branch Statistic (PBS)** quantifies allele-frequency differentiation along one branch of a three-population comparison using pairwise *F*ST estimates. This workflow calculates PBS across autosomal variants, identifies variants in the extreme empirical tail of the focal population's distribution, reduces linkage between candidates using PLINK, and prepares the resulting variants for annotation with Ensembl VEP.

The pipeline is designed for population-genomic research on Siberian populations, with West, Central, and Northeast Siberian groups compared against Southwest European and mainland East Asian reference populations. **High PBS indicates unusual population-specific differentiation; it is not sufficient evidence of positive selection on its own.**

## Analysis workflow

```text
Whole-genome VCF + combined focal-cohort VCF + population metadata
                              │
                   Validate samples and inputs
                              │
               Recalculate AC/AN; filter SNPs
                 (both allele counts ≥ 5)
                              │
                 Extract regional triplet VCFs
                              │
               Hudson FST → PBS by chromosome
                              │
            Merge chromosomes; empirical PBS ranks
                              │
               PLINK LD clumping of candidates
                              │
              Recover REF/ALT; validate with FASTA
                              │
                 GRCh37 candidate VCFs
                              │
                Offline Ensembl VEP (optional)
```

Snakemake resolves dependencies and runs independent chromosome/region combinations in parallel. The default research setup covers 22 autosomes and three focal regions. The workflow produces per-chromosome PBS tables, genome-wide empirical ranks, LD-clumped candidate sets, and indexed, reference-checked VCFs; VEP can add transcript-level consequence annotations.

### Implementation

| Tool | Role |
|---|---|
| **Python** (NumPy, pandas, cyvcf2) | Allele counting, Hudson *F*ST, PBS, ranking and validation |
| **bcftools** | VCF filtering, cohort extraction, allele-count recalculation and reference checks |
| **PLINK 1.9** | LD-based candidate clumping |
| **Ensembl VEP** | Offline variant-consequence annotation with a GRCh37 cache |
| **Snakemake / Conda** | Dependencies, reproducible execution and software environments |
| **SLURM** | Distribution of jobs across HPC nodes |

## Repository structure

```text
Snakefile                        Workflow rules and resource requirements
config.yaml                      Configuration for the bundled test dataset
config.production.example.yaml   Template for a real-data analysis
scripts/
  pbs_calculations.py            Allele-frequency and PBS calculations
  core.py                        Validation and result-processing utilities
  stage.py                       Execution of workflow stages
envs/
  analysis.yaml                  Dependencies for workflow jobs
  portfolio.yaml                 Combined environment for local testing
profiles/slurm/config.yaml      Generic SLURM executor profile
tests/                          Unit tests, data generator, integration checks
VALIDATION.md                   Tests completed and remaining coverage
```

## Quick start

Install [Conda](https://docs.conda.io/) or Mamba on Linux, then run from the repository root:

```bash
conda env create -f envs/portfolio.yaml
conda activate pbs-portfolio
bash tests/smoke.sh
```

This command generates a small **synthetic test dataset** and runs the unit tests and end-to-end integration tests, including variant filtering, PBS calculation, empirical ranking, LD clumping, VCF construction and reference checks. It also exercises the empty-variant case. No individual-level human genotypes are distributed in this repository.

The default `config.yaml` uses the test dataset, two short contigs, VEP disabled, and an intentionally permissive empirical candidate threshold (`0.8`) to exercise downstream stages. **It is not a production configuration.** Generated files are written to `tests/data/` and `results/`, both excluded from Git.

The combined `envs/portfolio.yaml` recipe is provided for installation; successful creation of that exact environment from scratch has not yet been verified. See [Validation](VALIDATION.md) for tested environments and commands.

## Running on real data

Create a private configuration from the production template:

```bash
cp config.production.example.yaml config.local.yaml
```

Update the VCFs and their existing indexes, GRCh37 FASTA (`1`–`22` contigs with a `.fai` index), population metadata and resource settings. The `regions` section expects tab-separated sample-to-population maps and three-row cyclic triplet files for each focal region. The combined focal-cohort VCF must contain the union of all mapped focal samples. Input files are read-only; outputs must be placed in a new subdirectory of `results/`.

**Key analysis settings** in `config.production.example.yaml`:

| Parameter | Research setting |
|---|---|
| Input variants | Biallelic A/C/G/T SNPs |
| Allele-count criterion | `AC ≥ 5` and `AN − AC ≥ 5` in the combined focal cohort |
| PBS ranking | Descending per-region rank / number of non-missing values |
| Candidate threshold | Empirical rank fraction `< 0.001` |
| LD reference | Combined focal cohort (`all_siberians`) |
| PLINK clumping | Candidate-only, `r² = 0.5`, `250 kb` |
| *F*ST policy | `corrected_clipped` |
| Reference | GRCh37, numeric autosomes |

The empirical rank fraction is a **prioritisation score**, not a calibrated hypothesis-test *P* value. The `corrected_clipped` policy constrains Hudson *F*ST to `[0, 0.999999]` before PBS calculation. An alternative `archive_legacy` policy is available for comparison, but includes a known opposite-fixation edge case. Choice of *F*ST handling, LD reference and clumping strategy can affect candidate sets.

Review the planned jobs and then run input validation:

```bash
snakemake --configfile config.local.yaml --cores 2 --dry-run
snakemake --configfile config.local.yaml --cores 1 results/production/preflight.json
```

The preflight path above assumes `output_dir: results/production`. After checking the inputs and resource limits, the full workflow can be started with:

```bash
snakemake --configfile config.local.yaml --cores 2
```

To let Snakemake deploy each rule's `envs/analysis.yaml`, add `--software-deployment-method conda` to the corresponding command. The reference FASTA and input VCF indexes must exist before starting.

### Functional annotation

VEP is optional and disabled in the bundled configurations. Enabling `vep.enabled` requires an existing VEP installation, an offline **GRCh37** cache matching its version, and an accessible reference FASTA. The workflow writes annotated VCFs with transcript consequences in the `CSQ` field. Optional CADD annotation requires compatible local plugin and indexed score files. Neither annotation resource is bundled with the code.

### HPC execution

A generic SLURM executor profile is included under `profiles/slurm/`. Copy it to the ignored `profiles/local/` directory and add settings required by your HPC site:

```bash
cp -r profiles/slurm profiles/local
# Edit profiles/local/config.yaml with your account, partition or cluster
snakemake --profile profiles/local --configfile config.local.yaml --dry-run
```

After reviewing the job graph, run the desired target without `--dry-run`. Snakemake can execute chromosome-level tasks across compute nodes; the controller must run in a location approved by the HPC centre. Full-genome resource requirements depend on cohort size and should be assessed before production execution.

## Validation

Workflow validation was performed on the **VSC wICE cluster** (KU Leuven Tier-2 HPC). Automated unit tests and end-to-end integration testing covered VCF preprocessing, population-specific allele frequency estimation, PBS calculations, empirical significance assessment, LD clumping, and candidate variant generation.

Additional checks using real whole-genome sequencing data covered cohort configuration, sample selection, and SLURM job execution. Functional annotation was evaluated independently using offline Ensembl VEP with the GRCh37 reference assembly.

See [VALIDATION.md](VALIDATION.md) for the tested scope and software versions.

## References and tool documentation

[bcftools](https://samtools.github.io/bcftools/bcftools) · [PLINK 1.9 clumping](https://www.cog-genomics.org/plink/1.9/postproc) · [Snakemake SLURM executor](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) · [Ensembl VEP](https://grch37.ensembl.org/info/docs/tools/vep/script/vep_options.html)

## Example analyses

The repository includes two examples illustrating the main computational components of the workflow.

**[PBS analysis](examples/pbs/)** demonstrates population-specific allele frequencies, PBS calculations, empirical ranking, LD clumping, and candidate variant generation using simulated genotype data. The example includes input and output VCF files, PBS statistics, and prioritised variants.

**[Functional annotation](examples/vep/)** demonstrates offline Ensembl VEP annotation of GRCh37-compatible variants. The example includes input and annotated VCF files, together with a summary of predicted coding and noncoding consequences.

Both examples provide small, directly inspectable input and output files. The VEP demonstration uses GRCh37 coordinates, while the PBS integration test uses an artificial reference assembly.
