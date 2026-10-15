from pathlib import Path
import json
import pandas as pd
from cyvcf2 import VCF
p=Path('results/test_empty')
for r in ['NE','West','Central']:
    assert pd.read_csv(p/f'pbs/{r}/merged_pval.tsv',sep='\t').empty
    assert pd.read_csv(p/f'final/{r}.leads.tsv',sep='\t').empty
    assert json.loads((p/f'clumps/{r}.status.json').read_text())['status']=='no_candidates'
    assert json.loads((p/f'final/{r}.validated.json').read_text())['records']==0
    v=VCF(str(p/f'final/{r}.GRCh37.vcf.gz'))
    assert list(v)==[]
    v.close()
print('PASS empty-cohort path: header-only PBS, no PLINK clumping, valid indexed empty VCF')
