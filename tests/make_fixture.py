"""Generate a tiny artificial reference and VCF; no human sequence or WGS input."""
from pathlib import Path
import random
import subprocess
import yaml
p=Path('tests/data')
p.mkdir(parents=True,exist_ok=True)
fasta=''
fai=[]
for c in ['1','2']:
    fasta += f'>{c}\n'
    offset=len(fasta)
    fasta += ('A'*50+'\n')*20
    fai.append(f'{c}\t1000\t{offset}\t50\t51\n')
(p/'reference.fa').write_text(fasta)
(p/'reference.fa.fai').write_text(''.join(fai))
regions={'NE':'NE_Siberia','West':'West_Siberia','Central':'Central_Siberia'}
samples=[f'{r}{i}' for r in regions for i in range(1,4)]+['SW1','SW2','EA1','EA2']
rows=[(s,regions[s[:-1]]) for s in samples[:9]]+[(s,'SW_Europe' if s.startswith('SW') else 'Mainland_E_Asia') for s in samples[9:]]
for r,pop in regions.items():
    selected=[s for s,label in rows if label in (pop,'SW_Europe','Mainland_E_Asia')]
    (p/f'{r}.samples').write_text('\n'.join(selected)+'\n')
    (p/f'{r}.map.tsv').write_text(''.join(f'{s}\t{label}\n' for s,label in rows))
    (p/f'{r}.triplets.tsv').write_text(f'{pop}\tSW_Europe\tMainland_E_Asia\nSW_Europe\tMainland_E_Asia\t{pop}\nMainland_E_Asia\t{pop}\tSW_Europe\n')
(p/'siberians.samples').write_text('\n'.join(samples[:9])+'\n')
header='##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##contig=<ID=2,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##INFO=<ID=AC,Number=A,Type=Integer,Description="Stale count for test">\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Stale number for test">\n'
header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n'
rng=random.Random(12)
records=[]
for c in ['1','2']:
    for i in range(1,13):
        if i==1:
            gt=['0/0']*13 # MAC zero -> removed
        elif i==2:
            gt=['1/1']*3+['0/0']*3+['0/1']*3+['0/0']*2+['0/1']*2 # AC=9 -> retained
        elif i==3:
            gt=['./.']*3+['0/1']*6+['0/0']*2+['0/1']*2 # MAC 6 despite missing NE
        elif i==4:
            gt=['0/1']*5+['0/0']*8 # MAC exactly 5 -> retained
        else:
            gt=[rng.choice(['0/0','0/1','1/1']) for _ in samples]
        records.append('\t'.join([c,str(i*20),f'old_{c}_{i}','A','G','.','PASS','AC=0;AN=0','GT']+gt)+'\n')
(p/'full.vcf').write_text(header+''.join(records))
subprocess.run(['bcftools','view','-Oz','-o',str(p/'full.vcf.gz'),str(p/'full.vcf')],check=True)
subprocess.run(['bcftools','index','-t','-f',str(p/'full.vcf.gz')],check=True)
subprocess.run(['bcftools','view','-S',str(p/'siberians.samples'),'-Oz','-o',str(p/'siberians.vcf.gz'),str(p/'full.vcf.gz')],check=True)
subprocess.run(['bcftools','index','-t','-f',str(p/'siberians.vcf.gz')],check=True)
cfg=yaml.safe_load(Path('config.yaml').read_text())
cfg.update(input_vcf=str(p/'full.vcf.gz'),input_index=str(p/'full.vcf.gz.tbi'),mac_vcf=str(p/'siberians.vcf.gz'),mac_index=str(p/'siberians.vcf.gz.tbi'),reference_fasta=str(p/'reference.fa'),output_dir='results/test',chromosomes=['1','2'],test_mode=True)
for r in regions:
    cfg['regions'][r].update(samples=str(p/f'{r}.samples'),popmap=str(p/f'{r}.map.tsv'),triplets=str(p/f'{r}.triplets.tsv'))
cfg['vep']['enabled']=False
cfg['clumping']['candidate_p']=0.8
cfg['resources'].update(mem_mb=1024,runtime_minutes=5,tool_threads=1)
Path('tests/config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
empty=dict(cfg,mac_min=100,output_dir='results/test_empty')
Path('tests/config.empty.yaml').write_text(yaml.safe_dump(empty,sort_keys=False))
