# Validation and testing

Testing covers the core computational workflow, input validation on real WGS data, and submission of selected Snakemake jobs through SLURM. Individual-level human genomic data are not included in this repository.

## Unit and integration tests

**15 October 2026 | VSC wICE (KU Leuven Tier-2 HPC)**

The public project archive was extracted into a fresh directory, and `bash tests/smoke.sh` completed successfully using an existing `popgen-analysis` Conda environment. The test suite generates 24 variants across two artificial contigs and 13 artificial samples, then evaluates:

- Hudson *F*ST, PBS calculations, allele counts, missing genotypes and empirical rank ties.
- AC/AN recalculation, allele-count filtering and three-population VCF extraction.
- Chromosome-wise PBS output, numerical chromosome merging and PLINK LD clumping.
- Source REF/ALT matching, reference-checked and indexed VCF output.
- Rejection of inconsistent input data, handling of an empty candidate set and repeat-execution checks.

The integration test completed successfully, including a **36/36-step Snakemake run** and verification of the empty-variant workflow. The empty case produced valid indexed VCFs without invoking PLINK clumping. The test script also invokes Snakemake dry runs to inspect the dependency graph and resumability; it does not enforce a zero-pending-jobs assertion.

The installed analysis environment was reused. Creation and execution of the combined `envs/portfolio.yaml` environment from scratch have not yet been tested. VEP and CADD command arguments are unit-tested but their tools are not invoked by this integration test.

## Selected checks with real WGS data

**15 October 2026 | VSC wICE**

The production input-validation rule (`preflight`) completed successfully using real WGS data. It checked VCF membership, the population maps and triplet definitions, GRCh37 reference compatibility and availability of the matching VEP cache.

| Focal population | Focal samples | Total samples in triplet |
|---|---:|---:|
| Northeast Siberia | 25 | 60 |
| West Siberia | 18 | 53 |
| Central Siberia | 31 | 66 |

There were **74 focal samples** across the three regions. Each triplet shared 17 Southwest European and 18 mainland East Asian reference samples.

A separate offline **Ensembl VEP 109.3** run annotated **84/84** SNPs from real WGS data using a GRCh37 cache. All variants were retained and carried nonempty `CSQ` fields. This confirms standalone annotation, not execution of the Snakemake-managed VEP rule.

## SLURM execution

Two jobs were submitted and completed successfully through Snakemake's SLURM executor on VSC wICE:

| Rule | Result |
|---|---|
| `preflight` | Completed successfully |
| `sample_lists` (Northeast) | Completed; generated 60 triplet and 25 focal sample IDs |

The run used **Snakemake 9.22.0**, **SLURM executor plugin 2.8.0**, and **Conda 26.9.3**. It used the site-specific execution profile. The generic `profiles/slurm/` template in this repository has not been independently submitted to an HPC scheduler.

## Scope

The **end-to-end PBS-to-VCF analysis was tested using generated data**. On real WGS data, validation covered inputs, sample extraction and job submission, along with a standalone VEP test. A full genome-wide PBS run with this packaged workflow, the Snakemake-managed VEP rule, and optional CADD annotation have not yet been executed end to end. The source layout in this release has passed static checks; the full integration test should be repeated after these packaging changes.
