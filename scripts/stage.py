"""Snakemake stage dispatcher; external commands use argument lists, never shell interpolation."""
import contextlib
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

sys.path.insert(0,str(Path(__file__).resolve().parent))
from core import *

sm = snakemake
cfg = sm.config
I,O = sm.input,sm.output
region = getattr(sm.wildcards,'region',None)
chrom = getattr(sm.wildcards,'chrom',None)


def execute(args, capture=False):
    args = list(map(str,args))
    print('+', ' '.join(args), flush=True)
    return subprocess.run(args,check=True,text=True,stdout=subprocess.PIPE if capture else LOG,stderr=LOG).stdout


def pipe(commands):
    procs = []
    previous = None
    try:
        for j,args in enumerate(commands):
            print('+', ' '.join(map(str,args)),flush=True)
            p = subprocess.Popen(list(map(str,args)),stdin=previous,stdout=subprocess.PIPE if j < len(commands)-1 else LOG,stderr=LOG)
            if previous is not None:
                previous.close()
            previous = p.stdout
            procs.append(p)
        codes = [p.wait() for p in procs]
        if any(codes):
            raise RuntimeError(f'External pipeline failed: {codes}')
    finally:
        for p in procs:
            if p.poll() is None:
                p.terminate()
                p.wait()


def index(path):
    execute(['bcftools','index','--tbi','--force',path])


def vcf_samples(path):
    v = VCF(str(path))
    try:
        return list(v.samples)
    finally:
        v.close()


def expanded(path):
    path = os.path.expanduser(os.path.expandvars(path))
    if '$' in path:
        raise ValueError(f'Unset environment variable in {path}')
    return Path(path)


def vep_preflight():
    v = cfg['vep']
    cache = expanded(v['cache_dir'])
    helptext = execute(['conda','run','--no-capture-output','-n',v['conda_env'],'vep','--help'],capture=True)
    match = re.search(r'ensembl-vep\s*[:=]\s*(\d+)',helptext,re.I)
    if not match:
        raise ValueError('Cannot determine installed VEP major version from vep --help')
    major = int(match.group(1))
    version = v.get('cache_version') or major
    if int(version) != major or not (cache/'homo_sapiens'/f'{version}_GRCh37').is_dir():
        raise ValueError('Require matching installed VEP version and homo_sapiens/<version>_GRCh37 cache')
    if v['cadd']:
        plugins = expanded(v['plugins_dir'])
        if not (plugins/'CADD.pm').is_file():
            raise ValueError('CADD.pm plugin missing')
        for key in ['cadd_snv','cadd_indels']:
            p = expanded(v[key])
            if not p.is_file() or not Path(str(p)+'.tbi').is_file():
                raise ValueError(f'CADD data or tabix index missing: {p}')
            with gzip.open(p,'rt') as f:
                lines=[]
                for line in f:
                    if not line.startswith('#'):
                        break
                    lines.append(line)
            if not any('GRCh37' in line for line in lines):
                raise ValueError(f'CADD header does not identify GRCh37: {p}')
        execute(['conda','run','--no-capture-output','-n',v['conda_env'],'tabix','--version'])
    return {'installed_major':major,'cache_version':int(version),'cache_dir':str(cache)}


def preflight():
    if cfg['reference_fasta'].endswith('.gz') and not Path(cfg['reference_fasta']+'.gzi').is_file():
        raise ValueError('Compressed reference requires existing .gzi index; inputs are never indexed in place')
    samples = vcf_samples(cfg['input_vcf'])
    mac_samples = vcf_samples(cfg['mac_vcf'])
    if not mac_samples or not set(mac_samples) <= set(samples):
        raise ValueError('Merged focal cohort is empty or contains samples absent from the full input VCF')
    focal_ids = set()
    for source in [cfg['input_vcf'],cfg['mac_vcf']]:
        v=VCF(str(source))
        try:
            if not set(cfg['chromosomes']) <= set(v.seqnames):
                raise ValueError('Source VCF lacks configured numeric chromosome contigs')
        finally:
            v.close()
    region_info = {}
    for r in cfg['regions']:
        mapping,trips,selected = metadata(cfg,r,samples)
        group=set(s for s in selected if mapping[s] == cfg['regions'][r]['focal'])
        if group & focal_ids:
            raise ValueError('A sample is assigned to multiple focal Siberian groups')
        focal_ids.update(group)
        region_info[r]={'triplets':trips,'n_selected':len(selected),'n_focal':sum(mapping[s]==cfg['regions'][r]['focal'] for s in selected)}
    if set(mac_samples) != focal_ids:
        raise ValueError('Merged Siberian VCF sample set differs from union of mapped focal Siberian samples; reconcile population maps before running')
    lengths = dict((r[0],int(r[1])) for r in table(cfg['reference_fasta']+'.fai',5))
    if cfg['reference_assembly'] != 'GRCh37':
        raise ValueError('Production pipeline requires GRCh37')
    for c in cfg['chromosomes']:
        if c not in lengths or (not cfg.get('test_mode',False) and lengths[c] != GRCH37[c]):
            raise ValueError(f'Reference contig {c} absent or inconsistent with GRCh37; require numeric contig names')
    for tool in ['bcftools','plink']:
        execute([tool,'--version'])
    info={'regions':region_info,'mac_samples':len(mac_samples),'config':cfg,'fst_policy':cfg['fst_policy'],'test_mode':cfg.get('test_mode',False)}
    if cfg['vep']['enabled']:
        if cfg.get('test_mode',False):
            raise ValueError('Synthetic test mode cannot claim real VEP annotation')
        info['vep']=vep_preflight()
    info['script_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.py')}
    Path(O[0]).write_text(json.dumps(info,indent=2))


def samples_stage():
    mapping,trips,selected = metadata(cfg,region,vcf_samples(cfg['input_vcf']))
    Path(O.samples).write_text('\n'.join(selected)+'\n')
    Path(O.focal).write_text('\n'.join(s for s in selected if mapping[s]==cfg['regions'][region]['focal'])+'\n')


def mac():
    threshold = cfg['mac_min']
    expr = f'INFO/AC>={threshold} && (INFO/AN-INFO/AC)>={threshold}'
    pipe([
        ['bcftools','view','-r',chrom,'-m2','-M2','-v','snps','-Ou',I.vcf],
        ['bcftools','+fill-tags','-Ou','--','-t','AC,AN'],
        ['bcftools','view','-i',expr,'-Ou'],
        ['bcftools','norm','-f',I.fasta,'-c','e','-Ou'],
        ['bcftools','annotate','--set-id','%CHROM:%POS','-Oz','-o',O.vcf]
    ])
    n = validate_records(O.vcf,chrom)
    index(O.vcf)
    positions=execute(['bcftools','query','-f','%CHROM\t%POS\n',O.vcf],capture=True)
    Path(O.positions).write_text(positions)
    Path(O.nrecords).write_text(str(n)+'\n')


def extract():
    expected=execute(['bcftools','query','-f','%CHROM\t%POS\t%REF\t%ALT\n',I.mac],capture=True)
    # Empty regions files are handled explicitly rather than sent as ambiguous empty -R input.
    commands=[['bcftools','view','-r',chrom,'-m2','-M2','-v','snps','-S',I.samples,'-Ou',I.vcf]]
    if Path(I.positions).stat().st_size:
        commands.append(['bcftools','view','-T',I.positions,'--targets-overlap','0','-Ou'])
    else:
        commands.append(['bcftools','view','-i','0','-Ou'])
    commands += [['bcftools','norm','-f',I.fasta,'-c','e','-Ou'],['bcftools','annotate','--set-id','%CHROM:%POS','-Oz','-o',O.vcf]]
    pipe(commands)
    validate_records(O.vcf,chrom)
    actual=execute(['bcftools','query','-f','%CHROM\t%POS\t%REF\t%ALT\n',O.vcf],capture=True)
    if actual != expected:
        raise ValueError('Triplet loci/alleles differ from MAC Siberian VCF; resolve duplicates or source inconsistencies')
    if set(vcf_samples(O.vcf)) != set(Path(I.samples).read_text().splitlines()):
        raise ValueError('Extracted sample set mismatch')
    index(O.vcf)


def pbs():
    n=run_pbs(I.vcf,I.popmap,I.triplets,O[0],cfg['fst_policy'],chrom)
    expected=int(Path(I.nrecords).read_text())
    if n != expected:
        raise ValueError('PBS row count differs from MAC reference')


def merge():
    merge_tables(I.parts,cfg['chromosomes'],O[0])


def concat():
    execute(['bcftools','concat','-Oz','-o',O.vcf,*I.vcfs])
    index(O.vcf)


def bfile():
    n=sum(1 for line in execute(['bcftools','query','-f','%ID\n',I.vcf],capture=True).splitlines())
    if n == 0:
        # A marked empty reference avoids PLINK's fatal "no variants" error.
        for f in [O.bed,O.bim,O.fam]:
            Path(f).touch()
        Path(O.status).write_text(json.dumps({'variants':0,'cohort':cfg['ld_samples'],'status':'empty_reference'}))
        return
    args=['plink','--vcf',I.vcf,'--double-id','--keep-allele-order','--make-bed','--threads',str(sm.threads),'--memory',str(sm.resources.mem_mb),'--out',str(O.bed)[:-4]]
    if region and cfg['ld_samples']=='focal':
        args += ['--keep',I.keep]
    execute(args)
    bim=pd.read_csv(O.bim,sep=r'\s+',header=None,dtype={0:str,1:str})
    ids=bim[0]+':'+bim[3].astype(str)
    if ids.duplicated().any() or not (ids==bim[1]).all():
        raise ValueError('PLINK BIM IDs do not match unique CHROM:POS keys')
    Path(O.status).write_text(json.dumps({'variants':n,'cohort':cfg['ld_samples']}))


def clump():
    df=pd.read_csv(I.pvals,sep='\t',dtype={'CHROM':str})
    bim=pd.read_csv(I.bim,sep=r'\s+',header=None) if Path(I.bim).stat().st_size else pd.DataFrame(columns=range(6))
    opts=cfg['clumping']
    chosen=pval_report(df,cfg['regions'][region]['focal'],bim[1].astype(str),O.pvals,opts['candidate_p'],opts['candidates_only'])
    focalcol=[c for c in df if c.startswith('Pval_PBS_'+cfg['regions'][region]['focal']+'_')][0]
    sig=df.loc[df[focalcol]<opts['candidate_p']]
    if sig.empty:
        Path(O.leads).write_text('')
        Path(O.clumped).write_text('CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2\n')
        Path(O.status).write_text(json.dumps({'status':'no_candidates','leads':0}))
        return
    prefix=str(O.clumped)[:-8]
    # Match the documented candidate-only reference. Optional full reference retains secondary neighbors.
    bprefix=str(I.bim)[:-4]
    if opts['candidates_only']:
        idsfile=Path(str(O.pvals)+'.ids')
        idsfile.write_text('\n'.join(chosen.SNP)+'\n')
        subset=prefix+'.candidates'
        execute(['plink','--bfile',bprefix,'--extract',idsfile,'--make-bed','--keep-allele-order','--out',subset,'--threads',sm.threads,'--memory',sm.resources.mem_mb])
        bprefix=subset
    execute(['plink','--bfile',bprefix,'--clump',O.pvals,'--clump-p1',opts['p1'],'--clump-p2',opts['p2'],'--clump-r2',opts['r2'],'--clump-kb',opts['kb'],'--clump-snp-field','SNP','--clump-field','P','--out',prefix,'--threads',sm.threads,'--memory',sm.resources.mem_mb])
    if not Path(O.clumped).is_file():
        raise ValueError('PLINK succeeded but did not produce clumps despite candidate variants; inspect log')
    lead_ids=[]
    # Whitespace parser handles PLINK report padding without fragile grep or line numbers.
    lines=[x.split() for x in Path(O.clumped).read_text().splitlines() if x.strip()]
    if not lines or 'SNP' not in lines[0]:
        raise ValueError('Unexpected PLINK clump header')
    ix=lines[0].index('SNP')
    lead_ids=[r[ix] for r in lines[1:]]
    if not set(lead_ids) <= set(chosen.SNP):
        raise ValueError('PLINK returned unmatched lead IDs')
    Path(O.leads).write_text(''.join(x+'\n' for x in lead_ids))
    Path(O.status).write_text(json.dumps({'status':'clumped','leads':len(lead_ids)}))


def leads():
    ids=Path(I.leads).read_text().splitlines()
    annotate_leads(I.pvals,ids,I.vcfs,O[0])


def vcf_stage():
    raw=Path(str(O.vcf)+'.raw.vcf')
    write_vcf(I.tsv,raw,I.fai,'synthetic_fixture' if cfg.get('test_mode',False) else 'GRCh37')
    try:
        # -c e rejects reference mismatch; never swaps or silently repairs REF/ALT.
        execute(['bcftools','norm','-f',I.fasta,'-c','e','-Oz','-o',O.vcf,raw])
        execute(['bcftools','view','-Ov','-o',os.devnull,O.vcf])
        index(O.vcf)
        n=len(pd.read_csv(I.tsv,sep='\t'))
        rows=execute(['bcftools','query','-f','%ID\n',O.vcf],capture=True).splitlines()
        if len(rows)!=n or len(set(rows))!=n:
            raise ValueError('Final VCF row count or ID mismatch')
        Path(O.validation).write_text(json.dumps({'assembly':'synthetic_fixture' if cfg.get('test_mode',False) else 'GRCh37','reference_checked':True,'records':n},indent=2))
    finally:
        raw.unlink(missing_ok=True)


def vep_stage():
    n=int(json.loads(Path(I.validation).read_text())['records'])
    if not n:
        shutil.copyfile(I.vcf,O.vcf)
        index(O.vcf)
        Path(O.status).write_text(json.dumps({'status':'skipped_no_leads'}))
        return
    report=json.loads(Path(I.preflight).read_text())['vep']
    v=cfg['vep']
    args=vep_command(v,report,I.fasta,I.vcf,O.vcf,sm.threads)
    execute(args)
    execute(['bcftools','view','-Ov','-o',os.devnull,O.vcf])
    before=execute(['bcftools','query','-f','%CHROM:%POS:%REF:%ALT\n',I.vcf],capture=True)
    after=execute(['bcftools','query','-f','%CHROM:%POS:%REF:%ALT\n',O.vcf],capture=True)
    if before != after:
        raise ValueError('VEP dropped or changed input variants; inspect VEP warnings')
    header=execute(['bcftools','view','-h',O.vcf],capture=True)
    if 'ID=CSQ,' not in header or (v['cadd'] and 'CADD_PHRED' not in header):
        raise ValueError('Expected VEP/CADD annotation fields missing')
    index(O.vcf)
    Path(O.status).write_text(json.dumps({'status':'annotated','records':n,'cadd':v['cadd'],'version':report['cache_version']},indent=2))


def focal_keep():
    ids=Path(I.samples).read_text().splitlines()
    Path(O[0]).write_text(''.join(f'{s}\t{s}\n' for s in ids))


actions={'preflight':preflight,'samples':samples_stage,'mac':mac,'extract':extract,'pbs':pbs,'merge':merge,'empirical':lambda:empirical(I[0],O[0]),'concat':concat,'bfile':bfile,'clump':clump,'leads':leads,'vcf':vcf_stage,'vep':vep_stage,'focal_keep':focal_keep}
for f in O:
    Path(f).parent.mkdir(parents=True,exist_ok=True)
Path(sm.log[0]).parent.mkdir(parents=True,exist_ok=True)
with open(sm.log[0],'w',buffering=1) as LOG:
    with contextlib.redirect_stdout(LOG),contextlib.redirect_stderr(LOG):
        try:
            actions[sm.params.action]()
        except Exception:
            traceback.print_exc()
            raise
