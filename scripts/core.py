"""Population-genomics calculations, input validation, and VCF utilities."""
import csv
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from cyvcf2 import VCF

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('pbs_engine', HERE/'pbs_calculations.py')
pbs_engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pbs_engine)
LEGACY_FST = pbs_engine.hudson_fst
GRCH37 = dict(zip(map(str, range(1, 23)), [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]))


def table(path, columns):
    with open(path, newline='') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    if not rows or any(len(r) != columns or any(not x.strip() for x in r) for r in rows):
        raise ValueError(f'{path}: require nonempty, headerless {columns}-column TSV')
    return rows


def metadata(cfg, region, vcf_samples):
    rc = cfg['regions'][region]
    rows = table(rc['popmap'], 2)
    ids = [x[0] for x in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f'{region}: duplicate sample IDs in population map')
    mapping = dict(rows)
    trips = table(rc['triplets'], 3)
    first = trips[0]
    expected = [first, [first[1],first[2],first[0]], [first[2],first[0],first[1]]]
    if first[0] != rc['focal'] or len(set(first)) != 3 or trips != expected:
        raise ValueError(f'{region}: triplets must contain exactly the three cyclic permutations, focal {rc["focal"]} first; found {trips}')
    if set(first[1:]) != {'SW_Europe','Mainland_E_Asia'}:
        raise ValueError(f'{region}: expected SW_Europe and Mainland_E_Asia comparison populations')
    selected = [s for s,p in rows if p in first]
    missing = set(selected) - set(vcf_samples)
    if missing:
        raise ValueError(f'{region}: mapped triplet samples absent from input VCF: {sorted(missing)}')
    if any(not any(mapping[s] == p for s in selected) for p in first):
        raise ValueError(f'{region}: population has no selected samples')
    # Optional sample lists must agree exactly with the population map.
    if rc.get('samples'):
        listed = [x[0] for x in table(rc['samples'], 1)]
        if len(listed) != len(set(listed)) or set(listed) != set(selected):
            raise ValueError(f'{region}: sample list differs from population map/triplets')
    return mapping, trips, selected


def schema(trips):
    return ['CHROM','POS'] + ['AF_'+p for p in trips[0]] + ['PBS_'+'_'.join(t) for t in trips]


def corrected_fst(ac1, ac2):
    """Hudson estimator with non-negative clipping and a finite upper bound."""
    n1,n2 = ac1.sum(), ac2.sum()
    if n1 < 2 or n2 < 2:
        return np.nan
    p1,p2 = ac1[1]/n1, ac2[1]/n2
    den = p1*(1-p2)+p2*(1-p1)
    if den == 0:
        return 0.0
    num = (p1-p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    return float(np.clip(num/den, 0.0, 0.999999))


def validate_records(path, chrom=None):
    """Reject ambiguous loci and unsupported genotypes before allele counting."""
    vcf = VCF(str(path))
    previous = 0
    count = 0
    try:
        for v in vcf:
            if chrom is not None and v.CHROM != str(chrom):
                raise ValueError('Unexpected chromosome')
            if v.POS <= previous:
                raise ValueError(f'{path}: unsorted or duplicate locus {v.CHROM}:{v.POS}; PBS output does not include allele keys')
            previous = v.POS
            if len(v.REF) != 1 or len(v.ALT) != 1 or len(v.ALT[0]) != 1 or v.REF not in 'ACGT' or v.ALT[0] not in 'ACGT':
                raise ValueError('PBS adaptation accepts biallelic A/C/G/T SNPs only')
            for g in v.genotypes:
                if len(g) != 3 or any(a not in (-1,0,1) for a in g[:2]):
                    raise ValueError(f'{path}: require diploid biallelic GT at {v.CHROM}:{v.POS}')
            count += 1
    finally:
        vcf.close()
    return count


def run_pbs(vcf, popmap, trips_path, output, policy, chrom):
    trips = table(trips_path, 3)
    n = validate_records(vcf, chrom)
    if policy not in ('corrected_clipped','archive_legacy'):
        raise ValueError('Unknown FST policy')
    pbs_engine.hudson_fst = corrected_fst if policy == 'corrected_clipped' else LEGACY_FST
    if n:
        pbs_engine.main(str(vcf), str(popmap), str(trips_path), str(output), 0, n)
        df = pd.read_csv(output, sep='\t')
        if list(df.columns) != schema(trips) or len(df) != n:
            raise ValueError('Unexpected PBS output schema or row count')
    else:
        pd.DataFrame(columns=schema(trips)).to_csv(output, sep='\t', index=False)
    return n


def merge_tables(paths, chromosomes, output):
    """Stream numeric chromosome order, including real headers for empty inputs."""
    if len(paths) != len(chromosomes):
        raise ValueError('Missing chromosome input')
    header = None
    total = 0
    with open(output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        for path,chrom in zip(paths, chromosomes):
            with open(path, newline='') as f:
                reader = csv.reader(f, delimiter='\t')
                current = next(reader, None)
                if not current or current[:2] != ['CHROM','POS']:
                    raise ValueError(f'{path}: invalid header')
                if header is None:
                    header = current
                    writer.writerow(header)
                elif current != header:
                    raise ValueError(f'{path}: inconsistent header')
                prev = 0
                for row in reader:
                    if len(row) != len(header) or row[0] != str(chrom) or int(row[1]) <= prev:
                        raise ValueError(f'{path}: wrong chromosome, duplicate, or unsorted row')
                    prev = int(row[1])
                    writer.writerow(row)
                    total += 1
    return total


def empirical(input_path, output):
    df = pd.read_csv(input_path, sep='\t', dtype={'CHROM':str})
    cols = [c for c in df if c.startswith('PBS_')]
    if len(cols) != 3:
        raise ValueError('Expected three PBS columns')
    for col in cols:
        df[col] = pd.to_numeric(df[col], errors='raise').astype(float)
        valid = df[col].notna()
        if not np.isfinite(df.loc[valid,col]).all():
            raise ValueError('Infinite PBS')
        n = int(valid.sum())
        df['Pval_'+col] = np.nan
        if n:
            df.loc[valid,'Pval_'+col] = df.loc[valid,col].rank(method='min',ascending=False)/n
    df.to_csv(output,sep='\t',index=False)


def pval_report(df, focal, all_ids, output, threshold, candidates_only):
    col = f'Pval_PBS_{focal}_SW_Europe_Mainland_E_Asia'
    # Accept either cyclic first-triplet comparison order by selecting the exact focal PBS column.
    candidates = [c for c in df if c.startswith('Pval_PBS_'+focal+'_')]
    if len(candidates) != 1:
        raise ValueError('Ambiguous focal PBS column')
    col = candidates[0]
    df = df.copy()
    df['SNP'] = df.CHROM.astype(str) + ':' + df.POS.astype(str)
    valid = df[col].notna()
    if candidates_only:
        valid &= df[col] < threshold
    chosen = df.loc[valid,['SNP',col]].rename(columns={col:'P'})
    if not set(chosen.SNP) <= set(all_ids):
        raise ValueError('P-value IDs absent from LD BIM reference')
    chosen.to_csv(output,sep='\t',index=False)
    return chosen


def annotate_leads(pval_path, ids, vcf_paths, output):
    df = pd.read_csv(pval_path,sep='\t',dtype={'CHROM':str})
    df['ID'] = df.CHROM + ':' + df.POS.astype(str)
    if df.ID.duplicated().any() or len(ids) != len(set(ids)) or not set(ids) <= set(df.ID):
        raise ValueError('Duplicate or unmatched lead ID')
    df = df[df.ID.isin(ids)].copy()
    needed = set(ids)
    alleles = {}
    for path in vcf_paths:
        vcf = VCF(str(path))
        try:
            for v in vcf:
                key = f'{v.CHROM}:{v.POS}'
                if key in needed:
                    if key in alleles or len(v.ALT) != 1:
                        raise ValueError('Ambiguous allele join')
                    if v.ID != key:
                        raise ValueError('VCF variant ID mismatch')
                    alleles[key] = (v.REF,v.ALT[0])
        finally:
            vcf.close()
    if set(alleles) != needed:
        raise ValueError('Missing lead alleles')
    df['REF'] = [alleles[k][0] for k in df.ID]
    df['ALT'] = [alleles[k][1] for k in df.ID]
    columns=['CHROM','POS','ID','REF','ALT']+[c for c in df if c not in ['CHROM','POS','ID','REF','ALT']]
    df[columns].to_csv(output,sep='\t',index=False)


def write_vcf(input_path, output, fai, assembly):
    df = pd.read_csv(input_path,sep='\t',dtype={'CHROM':str})
    contigs = dict((r[0],int(r[1])) for r in table(fai,5))
    df = df.sort_values(['CHROM','POS'],key=lambda x: x.astype(int))
    with open(output,'w') as f:
        f.write('##fileformat=VCFv4.2\n##source=pbs_selection_snakemake\n')
        f.write(f'##reference={assembly}\n')
        for c,n in contigs.items():
            f.write(f'##contig=<ID={c},length={n}>\n')
        f.write('##FILTER=<ID=PASS,Description="All workflow checks passed">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for r in df.itertuples():
            if r.CHROM not in contigs or not 0 < r.POS <= contigs[r.CHROM] or r.REF not in 'ACGT' or r.ALT not in 'ACGT':
                raise ValueError('Invalid SNP VCF record')
            f.write(f'{r.CHROM}\t{r.POS}\t{r.ID}\t{r.REF}\t{r.ALT}\t.\tPASS\t.\n')


def vep_command(vep_cfg, report, fasta, input_vcf, output_vcf, threads):
    """Build an offline-only command for the existing conda VEP installation."""
    import os
    def expand(value):
        value=os.path.expandvars(os.path.expanduser(str(value)))
        if '$' in value:
            raise ValueError('Unset environment variable in VEP path')
        return value
    args=['conda','run','--no-capture-output','-n',vep_cfg['conda_env'],'vep',
          '--offline','--cache','--species','homo_sapiens','--assembly','GRCh37',
          '--cache_version',str(report['cache_version']),'--dir_cache',str(report['cache_dir']),
          '--fasta',str(fasta),'--check_ref','--format','vcf','--vcf','--symbol',
          '--fork',str(threads),'--compress_output','bgzip','--input_file',str(input_vcf),
          '--output_file',str(output_vcf),'--force_overwrite','--no_stats']
    if vep_cfg['cadd']:
        args += ['--dir_plugins',expand(vep_cfg['plugins_dir']),'--plugin',
                 f"CADD,snv={expand(vep_cfg['cadd_snv'])},indels={expand(vep_cfg['cadd_indels'])}"]
    return args
