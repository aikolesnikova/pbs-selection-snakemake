import sys
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
import yaml
sys.path.insert(0,str(Path('scripts').resolve()))
from core import *

def fails(fn):
    try:
        fn()
    except (ValueError,FileNotFoundError):
        return
    raise AssertionError('Expected rejection')

# Independent algebra: p1=0.5 (4 alleles), p2=0.25 (4 alleles).
a,b=np.array([2,2]),np.array([3,1])
raw=((.5-.25)**2-.5*.5/3-.25*.75/3)/(.5*.75+.25*.5)
assert raw < 0 and corrected_fst(a,b)==0
assert LEGACY_FST(np.array([6,0]),np.array([0,6]))==0 # retained edge case in the legacy FST option
assert corrected_fst(np.array([6,0]),np.array([0,6]))==.999999
assert corrected_fst(np.array([6,0]),np.array([6,0]))==0
assert np.isnan(corrected_fst(np.array([0,0]),np.array([6,0])))
a,b=np.array([5,1]),np.array([1,5])
expected=((1/6-5/6)**2-2*(1/6)*(5/6)/5)/((1/6)*(1/6)+(5/6)*(5/6))
assert np.isclose(corrected_fst(a,b),expected)
assert pbs_engine.pop_ac(np.array([0,1,2,-1]),[0,1,2,3]).tolist()==[3,3]
with tempfile.TemporaryDirectory() as d:
    d=Path(d)
    cols=['CHROM','POS','PBS_A','PBS_B','PBS_C']
    df=pd.DataFrame([['2',10,2,0,np.nan],['2',20,2,0,np.nan],['2',30,1,0,np.nan]],columns=cols)
    df.to_csv(d/'p.tsv',sep='\t',index=False)
    empirical(d/'p.tsv',d/'out.tsv')
    out=pd.read_csv(d/'out.tsv',sep='\t')
    assert out.Pval_PBS_A.tolist()==[1/3,1/3,1]
    assert out.Pval_PBS_B.tolist()==[1/3]*3
    assert out.Pval_PBS_C.isna().all()
    pd.DataFrame(columns=cols).to_csv(d/'10.tsv',sep='\t',index=False)
    merge_tables([d/'p.tsv',d/'10.tsv'],['2','10'],d/'merged.tsv')
    assert len(pd.read_csv(d/'merged.tsv',sep='\t'))==3
    fails(lambda:merge_tables([d/'missing.tsv'],['1'],d/'bad.tsv'))
    fails(lambda:merge_tables([d/'p.tsv'],['10'],d/'bad.tsv'))
    # Unsupported genotypes and ambiguous loci must fail before PBS allele counting.
    lines=Path('tests/data/full.vcf').read_text().splitlines()
    header=[x for x in lines if x.startswith('#')]
    records=[x for x in lines if not x.startswith('#')]
    (d/'duplicate.vcf').write_text('\n'.join(header+[records[0],records[0]])+'\n')
    fails(lambda:validate_records(d/'duplicate.vcf','1'))
    cells=records[0].split('\t')
    cells[9]='0'
    (d/'haploid.vcf').write_text('\n'.join(header+['\t'.join(cells)])+'\n')
    fails(lambda:validate_records(d/'haploid.vcf','1'))
    import subprocess
    cells[9]='0/0'
    cells[3]='C' # artificial reference consists entirely of A
    (d/'wrongref.vcf').write_text('\n'.join(header+['\t'.join(cells)])+'\n')
    result=subprocess.run(['bcftools','norm','-f','tests/data/reference.fa','-c','e',str(d/'wrongref.vcf')],capture_output=True,text=True)
    assert result.returncode!=0 and 'mismatch' in result.stderr.lower()
    cfg=yaml.safe_load(Path('tests/config.yaml').read_text())
    fails(lambda:metadata(cfg,'NE',['nonexistent']))
    fake=df.rename(columns={'PBS_A':'Pval_PBS_A_B_C'})
    fake['Pval_PBS_A_B_C']=.2
    fails(lambda:pval_report(fake,'A',[],d/'bad.tsv',.9,True))
print('PASS unit checks: FST algebra and clipping, allele counts, empirical ties/NaNs, merge order/empty inputs, missing files/samples/IDs')

# Command construction only: no VEP/cache/plugin execution is claimed by this test.
v={'conda_env':'VEP','cadd':False}
args=vep_command(v,{'cache_version':115,'cache_dir':'/existing/cache'},'/ref.fa','in.vcf.gz','out.vcf.gz',2)
assert '--offline' in args and '--cache' in args and '--check_ref' in args
assert args[args.index('--assembly')+1]=='GRCh37' and '--plugin' not in args
v.update(cadd=True,plugins_dir='/existing/Plugins',cadd_snv='/existing/SNVs.tsv.gz',cadd_indels='/existing/InDels.tsv.gz')
args=vep_command(v,{'cache_version':115,'cache_dir':'/existing/cache'},'/ref.fa','in','out',1)
assert args[args.index('--plugin')+1]=='CADD,snv=/existing/SNVs.tsv.gz,indels=/existing/InDels.tsv.gz'
print('PASS offline VEP/CADD argument construction (VEP execution is outside this unit test)')
