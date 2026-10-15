# Functional annotation with Ensembl VEP

This example illustrates functional annotation of single-nucleotide
variants using Ensembl VEP 109.3 and the GRCh37 reference assembly.

Twelve demonstration variants were selected to illustrate different
predicted consequences, including coding, intronic, regulatory-proximal
and intergenic annotations.

Reference alleles were obtained from GRCh37. Alternative alleles were
assigned for demonstration purposes.

## Files

- `input.vcf`: Variants before annotation.
- `annotated.vcf`: Variants with VEP transcript-level annotations.
- `consequences.tsv`: Representative functional consequences and genes.

The `CSQ` field in `annotated.vcf` retains all transcript annotations.
The TSV provides one representative consequence per variant.

These variants are not presented as observed population polymorphisms
or evidence of natural selection. The example demonstrates annotation
functionality independently of the PBS integration-test dataset.
