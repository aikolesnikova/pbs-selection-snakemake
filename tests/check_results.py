import json
from pathlib import Path
import pandas as pd
from cyvcf2 import VCF
from numpy import allclose

base=Path('results/test')
for c in ['1','2']:
    v=VCF(str(base/f'mac/chr{c}.vcf.gz'))
    recs=list(v)
    assert any(x.POS==80 for x in recs) # MAC exactly five must be retained
    assert not any(x.POS==20 for x in recs)
    assert all(x.ID==f'{c}:{x.POS}' for x in recs)
    for x in recs:
        assert x.INFO.get('AC')>=5 and x.INFO.get('AN')-x.INFO.get('AC')>=5
    v.close()
    ids=pd.read_csv(base/'ld/all_siberians.bim',sep=r'\s+',header=None)[1].astype(str)
    assert {x.ID for x in recs}<=set(ids)
for r,focal in [('NE','NE_Siberia'),('West','West_Siberia'),('Central','Central_Siberia')]:
    merged=pd.read_csv(base/f'pbs/{r}/merged_pval.tsv',sep='\t',dtype={'CHROM':str})
    expected=sum(int((base/f'mac/chr{c}.nrecords.txt').read_text()) for c in ['1','2'])
    assert len(merged)==expected
    assert merged.CHROM.tolist()==sorted(merged.CHROM.tolist(),key=int)
    columns=[c for c in merged if c.startswith('PBS_')]
    assert len(columns)==3 and f'PBS_{focal}_SW_Europe_Mainland_E_Asia' in columns
    if r=='NE':
        assert merged.loc[merged.POS==60,columns].isna().all().all()
    for col in columns:
        valid=merged[col].notna()
        assert allclose(merged.loc[valid,'Pval_'+col],merged.loc[valid,col].rank(method='min',ascending=False)/valid.sum(),rtol=1e-12,atol=1e-15)
    leads=pd.read_csv(base/f'final/{r}.leads.tsv',sep='\t')
    assert len(leads)>0 and leads.ID.is_unique
    assert (leads.REF=='A').all() and (leads.ALT=='G').all()
    assert json.loads((base/f'final/{r}.validated.json').read_text())['records']==len(leads)
    vcf=VCF(str(base/f'final/{r}.GRCh37.vcf.gz'))
    assert '#CHROM\tPOS' in vcf.raw_header and '##CHROM\tPOS' not in vcf.raw_header
    assert '##reference=synthetic_fixture' in vcf.raw_header
    assert [x.ID for x in vcf]==leads.ID.tolist()
    vcf.close()
print('PASS real-tool integration: MAC >= 5, cohort tags, strict extraction, PBS/merge, ranks, PLINK IDs/clumps, allele annotation, reference-checked VCF')
