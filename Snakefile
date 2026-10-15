configfile: "config.yaml"
from pathlib import Path

OUT = config['output_dir'].rstrip('/')
if Path(OUT).is_absolute() or '..' in Path(OUT).parts or not OUT.startswith('results/'):
    raise ValueError('output_dir must be a new relative directory under results/; external input directories are read-only')
CHRS = [str(c) for c in config['chromosomes']]
if not CHRS or len(CHRS)!=len(set(CHRS)) or CHRS!=sorted(CHRS,key=int) or not set(CHRS)<=set(map(str,range(1,23))):
    raise ValueError('chromosomes must be unique numeric autosomes in ascending order')
REGIONS = list(config['regions'])
if set(REGIONS)!={'NE','West','Central'}:
    raise ValueError('Require NE, West and Central region definitions')
if config['ld_samples'] not in ('all_siberians','triplet','focal'):
    raise ValueError('Unknown LD sample policy')
if config['fst_policy'] not in ('corrected_clipped','archive_legacy'):
    raise ValueError('Unknown FST policy')
if type(config['mac_min']) is not int or config['mac_min']<1:
    raise ValueError('mac_min must be a positive integer')
CL=config['clumping']
if not (0<CL['candidate_p']<=1 and 0<CL['p1']<=1 and 0<CL['p2']<=1 and 0<=CL['r2']<=1 and CL['kb']>0):
    raise ValueError('Invalid clumping thresholds')
if CL['candidates_only'] and CL['p1']!=1:
    raise ValueError('Documented candidate-only mode prefilters P < candidate_p and uses p1=1')
if not CL['candidates_only'] and CL['p1']!=CL['candidate_p']:
    raise ValueError('Full-reference mode requires p1=candidate_p')
ENV='envs/analysis.yaml'
COMMON=['scripts/core.py','scripts/stage.py','scripts/pbs_calculations.py']
PREF=f'{OUT}/preflight.json'
VCF=config['input_vcf']
MAC=config['mac_vcf']
FA=config['reference_fasta']
MEM=config['resources']['mem_mb']
TIME=config['resources']['runtime_minutes']
T=config['resources']['tool_threads']
LD=lambda r: 'all_siberians' if config['ld_samples']=='all_siberians' else r

rule all:
    input:
        expand(f'{OUT}/final/{{region}}.validated.json',region=REGIONS),
        expand(f'{OUT}/vep/{{region}}.status.json',region=REGIONS) if config['vep']['enabled'] else []

rule preflight:
    input:
        vcf=VCF, tbi=config['input_index'], mac=MAC, mac_tbi=config['mac_index'],
        fasta=FA, fai=FA+'.fai',
        metadata=[config['regions'][r][k] for r in REGIONS for k in ('popmap','triplets')]+[config['regions'][r]['samples'] for r in REGIONS if config['regions'][r].get('samples')],
        helpers=COMMON
    output: PREF
    log: f'{OUT}/logs/preflight.log'
    params: action='preflight'
    threads: 1
    resources: mem_mb=4000, runtime=30
    conda: ENV
    script: 'scripts/stage.py'

rule sample_lists:
    input:
        preflight=PREF, vcf=VCF,
        popmap=lambda w:config['regions'][w.region]['popmap'],
        triplets=lambda w:config['regions'][w.region]['triplets'], helpers=COMMON
    output:
        samples=f'{OUT}/samples/{{region}}.txt', focal=f'{OUT}/samples/{{region}}.focal.txt'
    log: f'{OUT}/logs/samples_{{region}}.log'
    params: action='samples'
    threads: 1
    resources: mem_mb=1024, runtime=10
    conda: ENV
    script: 'scripts/stage.py'

rule mac_filter:
    input: preflight=PREF, vcf=MAC, tbi=config['mac_index'], fasta=FA, fai=FA+'.fai', helpers=COMMON
    output:
        vcf=f'{OUT}/mac/chr{{chrom}}.vcf.gz', tbi=f'{OUT}/mac/chr{{chrom}}.vcf.gz.tbi',
        positions=f'{OUT}/mac/chr{{chrom}}.positions.tsv', nrecords=f'{OUT}/mac/chr{{chrom}}.nrecords.txt'
    log: f'{OUT}/logs/mac_{{chrom}}.log'
    params: action='mac'
    threads: 1
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule extract_triplet:
    input:
        vcf=VCF, tbi=config['input_index'], samples=f'{OUT}/samples/{{region}}.txt',
        mac=f'{OUT}/mac/chr{{chrom}}.vcf.gz', positions=f'{OUT}/mac/chr{{chrom}}.positions.tsv',
        fasta=FA, fai=FA+'.fai', helpers=COMMON
    output: vcf=f'{OUT}/triplets/{{region}}/chr{{chrom}}.vcf.gz', tbi=f'{OUT}/triplets/{{region}}/chr{{chrom}}.vcf.gz.tbi'
    log: f'{OUT}/logs/extract_{{region}}_{{chrom}}.log'
    params: action='extract'
    threads: 1
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule pbs_chromosome:
    input:
        vcf=f'{OUT}/triplets/{{region}}/chr{{chrom}}.vcf.gz',
        nrecords=f'{OUT}/mac/chr{{chrom}}.nrecords.txt',
        popmap=lambda w:config['regions'][w.region]['popmap'],
        triplets=lambda w:config['regions'][w.region]['triplets'], helpers=COMMON
    output: f'{OUT}/pbs/{{region}}/chr{{chrom}}.tsv'
    log: f'{OUT}/logs/pbs_{{region}}_{{chrom}}.log'
    params: action='pbs'
    threads: 1
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule merge_chromosomes:
    input:
        parts=lambda w:expand(f'{OUT}/pbs/{{region}}/chr{{chrom}}.tsv',region=w.region,chrom=CHRS), helpers=COMMON
    output: f'{OUT}/pbs/{{region}}/merged.tsv'
    log: f'{OUT}/logs/merge_{{region}}.log'
    params: action='merge'
    threads: 1
    resources: mem_mb=2048, runtime=60
    conda: ENV
    script: 'scripts/stage.py'

rule empirical_pvalues:
    input: f'{OUT}/pbs/{{region}}/merged.tsv', helpers=COMMON
    output: f'{OUT}/pbs/{{region}}/merged_pval.tsv'
    log: f'{OUT}/logs/pvals_{{region}}.log'
    params: action='empirical'
    threads: 1
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule ld_vcf:
    input:
        vcfs=lambda w:expand(f'{OUT}/mac/chr{{chrom}}.vcf.gz',chrom=CHRS) if w.region=='all_siberians' else expand(f'{OUT}/triplets/{{region}}/chr{{chrom}}.vcf.gz',region=w.region,chrom=CHRS), helpers=COMMON
    output: vcf=f'{OUT}/ld/{{region}}.vcf.gz', tbi=f'{OUT}/ld/{{region}}.vcf.gz.tbi'
    log: f'{OUT}/logs/ld_vcf_{{region}}.log'
    params: action='concat'
    threads: 1
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule focal_keep:
    input: samples=f'{OUT}/samples/{{region}}.focal.txt', helpers=COMMON
    output: f'{OUT}/samples/{{region}}.keep'
    log: f'{OUT}/logs/keep_{{region}}.log'
    params: action='focal_keep'
    threads: 1
    resources: mem_mb=1024, runtime=5
    conda: ENV
    script: 'scripts/stage.py'

rule ld_bfile:
    input:
        vcf=f'{OUT}/ld/{{region}}.vcf.gz', tbi=f'{OUT}/ld/{{region}}.vcf.gz.tbi',
        keep=lambda w:f'{OUT}/samples/{w.region}.keep' if config['ld_samples']=='focal' else [], helpers=COMMON
    output:
        bed=f'{OUT}/ld/{{region}}.bed', bim=f'{OUT}/ld/{{region}}.bim', fam=f'{OUT}/ld/{{region}}.fam', status=f'{OUT}/ld/{{region}}.status.json'
    log: f'{OUT}/logs/ld_bfile_{{region}}.log'
    params: action='bfile'
    threads: T
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule clump:
    input:
        pvals=f'{OUT}/pbs/{{region}}/merged_pval.tsv',
        bed=lambda w:f'{OUT}/ld/{LD(w.region)}.bed',
        bim=lambda w:f'{OUT}/ld/{LD(w.region)}.bim',
        fam=lambda w:f'{OUT}/ld/{LD(w.region)}.fam', helpers=COMMON
    output:
        pvals=f'{OUT}/clumps/{{region}}.pvals.tsv', leads=f'{OUT}/clumps/{{region}}.leads.txt',
        clumped=f'{OUT}/clumps/{{region}}.clumped', status=f'{OUT}/clumps/{{region}}.status.json'
    log: f'{OUT}/logs/clump_{{region}}.log'
    params: action='clump'
    threads: T
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule annotate_alleles:
    input:
        pvals=f'{OUT}/pbs/{{region}}/merged_pval.tsv', leads=f'{OUT}/clumps/{{region}}.leads.txt',
        vcfs=lambda w:expand(f'{OUT}/mac/chr{{chrom}}.vcf.gz',chrom=CHRS), helpers=COMMON
    output: f'{OUT}/final/{{region}}.leads.tsv'
    log: f'{OUT}/logs/alleles_{{region}}.log'
    params: action='leads'
    threads: 1
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'

rule validated_vcf:
    input:
        tsv=f'{OUT}/final/{{region}}.leads.tsv', fasta=FA, fai=FA+'.fai', helpers=COMMON
    output:
        vcf=f'{OUT}/final/{{region}}.GRCh37.vcf.gz', tbi=f'{OUT}/final/{{region}}.GRCh37.vcf.gz.tbi', validation=f'{OUT}/final/{{region}}.validated.json'
    log: f'{OUT}/logs/vcf_{{region}}.log'
    params: action='vcf'
    threads: 1
    resources: mem_mb=4000, runtime=30
    conda: ENV
    script: 'scripts/stage.py'

rule offline_vep:
    input:
        vcf=f'{OUT}/final/{{region}}.GRCh37.vcf.gz', tbi=f'{OUT}/final/{{region}}.GRCh37.vcf.gz.tbi',
        validation=f'{OUT}/final/{{region}}.validated.json', preflight=PREF, fasta=FA, helpers=COMMON
    output:
        vcf=f'{OUT}/vep/{{region}}.vep.vcf.gz', tbi=f'{OUT}/vep/{{region}}.vep.vcf.gz.tbi', status=f'{OUT}/vep/{{region}}.status.json'
    log: f'{OUT}/logs/vep_{{region}}.log'
    params: action='vep'
    threads: config['resources']['vep_threads']
    resources: mem_mb=MEM, runtime=TIME
    conda: ENV
    script: 'scripts/stage.py'
