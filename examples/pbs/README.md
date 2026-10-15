# Population Branch Statistic (PBS) analysis

This example demonstrates population-specific genetic differentiation and candidate variant prioritisation using simulated genotype data.

The input VCF contains 24 variants across two artificial chromosomes and 13 samples. The workflow calculates allele frequencies and PBS statistics, assigns empirical P-values, performs LD clumping with PLINK, and generates a VCF of prioritised candidates.

## Example files

| File | Description |
|---|---|
| `input.vcf` | Input VCF containing simulated genotypes. |
| `NE.pbs.tsv` | PBS statistics and empirical P-values. |
| `NE.leads.tsv` | Candidate variants selected after LD clumping. |
| `NE.candidates.vcf` | Final candidate VCF. |

## Results

The example contains 24 input variants and produces 14 Northeast candidate variants after the prioritisation and clumping steps.

The integration test can be reproduced using:

`bash tests/smoke.sh`

The test uses an artificial reference assembly and a deliberately permissive candidate threshold (0.8) to exercise the complete analytical workflow. These results demonstrate the computational implementation and have no biological interpretation.
