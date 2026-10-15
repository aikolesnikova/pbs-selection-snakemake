#!/usr/bin/env python

import cyvcf2
import pandas as pd
import numpy as np
import argparse

def load_sample_pop_map(pop_file):
    df = pd.read_csv(pop_file, sep='\t', header=None, names=['ID', 'POP'])
    return df.set_index('ID')['POP'].to_dict()

def load_population_indices(samples, sample_to_pop, target_pop):
    return [i for i, s in enumerate(samples) if sample_to_pop.get(s) == target_pop]

def compute_pbs_from_fst(fst12, fst13, fst23):
    try:
        T12 = -np.log(1 - fst12) if fst12 < 1 else np.nan
        T13 = -np.log(1 - fst13) if fst13 < 1 else np.nan
        T23 = -np.log(1 - fst23) if fst23 < 1 else np.nan
        return (T12 + T13 - T23) / 2
    except:
        return np.nan

def hudson_fst(ac1, ac2):
    n1 = ac1.sum()
    n2 = ac2.sum()
    if n1 == 0 or n2 == 0:
        return np.nan
    p1 = ac1[1] / n1 if n1 > 0 else 0
    p2 = ac2[1] / n2 if n2 > 0 else 0
    # If both populations are monomorphic (all 0 or all 1), FST is defined as 0
    if (ac1[0] == n1 or ac1[1] == n1) and (ac2[0] == n2 or ac2[1] == n2):
        return 0.0
    num = (p1 - p2) ** 2 - (p1 * (1 - p1) / (n1 - 1) + p2 * (1 - p2) / (n2 - 1) if n1 > 1 and n2 > 1 else 0)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num / den if den != 0 else 0.0

def pop_ac(gtc, idx):
    # Get allele counts for 0, 1, 2 (number of alternate alleles)
    vals = [gtc[i] for i in idx if gtc[i] >= 0]
    if len(vals) == 0:
        return None  # all missing
    # Each value is 0, 1, or 2 (number of alt alleles for diploid)
    # Let's count ref and alt alleles
    n_ref = 2 * vals.count(0) + vals.count(1)
    n_alt = 2 * vals.count(2) + vals.count(1)
    total = n_ref + n_alt
    # Return array: [ref count, alt count]
    return np.array([n_ref, n_alt])

def pop_af(ac):
    n = ac.sum()
    if n == 0:
        return np.nan
    return ac[1] / n

def main(vcf_file, pop_file, triplet_file, output_file, start_idx, stop_idx):
    print(f"Processing SNPs from {start_idx} to {stop_idx}...")

    sample_to_pop = load_sample_pop_map(pop_file)
    triplets = pd.read_csv(triplet_file, sep='\t', header=None, names=['POP1', 'POP2', 'POP3'])
    vcf = cyvcf2.VCF(vcf_file)
    samples = list(vcf.samples)
    print(f"Loaded {len(samples)} samples from VCF.")

    # Choose the main triplet for frequency reporting (first triplet in file)
    freq_pops = triplets.iloc[0][['POP1','POP2','POP3']].values.tolist()

    out_rows = []
    snp_counter = 0
    for i, variant in enumerate(vcf):
        if i < start_idx:
            continue
        if i >= stop_idx:
            break
        gt = variant.genotypes
        # sum alleles per sample, -1 for missing
        gtc = np.array([
            x[0] + x[1] if x[0] >= 0 and x[1] >= 0 else -1
            for x in gt
        ])  # 0, 1, 2 or -1 for missing

        chrom = variant.CHROM
        pos = variant.POS

        row_data = {
            'CHROM': chrom,
            'POS': pos
        }

        # Calculate and add AF columns for the "main" triplet
        afs = []
        for pidx, pop in enumerate(freq_pops):
            idx = load_population_indices(samples, sample_to_pop, pop)
            ac = pop_ac(gtc, idx)
            if ac is not None:
                af = pop_af(ac)
            else:
                af = np.nan
            row_data[f'AF_{pop}'] = af
            afs.append(af)

        # PBS for all triplets as before
        for _, triplet in triplets.iterrows():
            pop1, pop2, pop3 = triplet['POP1'], triplet['POP2'], triplet['POP3']
            idx1 = load_population_indices(samples, sample_to_pop, pop1)
            idx2 = load_population_indices(samples, sample_to_pop, pop2)
            idx3 = load_population_indices(samples, sample_to_pop, pop3)
            if not idx1 or not idx2 or not idx3:
                row_data[f'PBS_{pop1}_{pop2}_{pop3}'] = np.nan
                continue
            ac1 = pop_ac(gtc, idx1)
            ac2 = pop_ac(gtc, idx2)
            ac3 = pop_ac(gtc, idx3)
            if ac1 is None or ac2 is None or ac3 is None:
                row_data[f'PBS_{pop1}_{pop2}_{pop3}'] = np.nan
                continue
            fst12 = hudson_fst(ac1, ac2)
            fst13 = hudson_fst(ac1, ac3)
            fst23 = hudson_fst(ac2, ac3)
            pbs_val = compute_pbs_from_fst(fst12, fst13, fst23)
            row_data[f'PBS_{pop1}_{pop2}_{pop3}'] = pbs_val

        out_rows.append(row_data)
        snp_counter += 1

        if snp_counter % 1000 == 0:
            print(f"  Processed {snp_counter} SNPs")

    df_out = pd.DataFrame(out_rows)
    df_out.to_csv(output_file, sep='\t', index=False)
    print(f"Finished. Wrote {len(df_out)} rows to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute per-site population branch statistics from VCF genotypes.")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--popfile", required=True, help="Sample to population file")
    parser.add_argument("--triplets", required=True, help="Triplet populations file")
    parser.add_argument("--out", required=True, help="Output file")
    parser.add_argument("--start", type=int, required=True, help="Start SNP index (inclusive)")
    parser.add_argument("--stop", type=int, required=True, help="Stop SNP index (exclusive)")
    args = parser.parse_args()
    main(args.vcf, args.popfile, args.triplets, args.out, args.start, args.stop)
