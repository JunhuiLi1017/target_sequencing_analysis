import argparse
import pandas as pd
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description='Extract ref/alt read counts and nucleotide counts from mpileup for variants.')
    parser.add_argument('--input', required=True, help='Input mpileup file from samtools mpileup.')
    parser.add_argument('--bed', required=True, help='BED file with variants (chrom, start, end, ref, alt).')
    parser.add_argument('--output', required=True, help='Output file for read counts (tab-separated).')
    return parser.parse_args()

def load_bed(bed_file):
    variants = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) >= 5:
                chrom = cols[0]
                start = int(cols[1])  # 0-based
                ref = cols[3].upper()
                alt = cols[4].upper()
                variants[(chrom, start + 1)] = {'ref': ref, 'alt': alt, 'end': int(cols[2])}
    return variants

def parse_pileup(pileup_file, variants):
    results = []
    pileup_data = {}
    with open(pileup_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue
            chrom = cols[0]
            pos = int(cols[1])
            pileup_data[(chrom, pos)] = {'ref_base': cols[2].upper(), 'bases': cols[4]}

    for (chrom, pos), variant in variants.items():
        ref_allele = variant['ref']
        alt_allele = variant['alt']
        end_pos = variant['end']  # 1-based

        ref_count = 0
        alt_count = 0

        a_count = 0
        c_count = 0
        g_count = 0
        t_count = 0

        is_indel = len(ref_allele) != len(alt_allele)

        if (chrom, pos) not in pileup_data:
            continue
        bases = pileup_data[(chrom, pos)]['bases']
        ref_base = pileup_data[(chrom, pos)]['ref_base']

        tokens = []
        i = 0
        while i < len(bases):
            if bases[i] in ['.', ',']:
                tokens.append('REF')
                i += 1
            elif bases[i] in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']:
                nt = bases[i].upper()
                tokens.append(nt)
                if nt == 'A':
                    a_count += 1
                elif nt == 'C':
                    c_count += 1
                elif nt == 'G':
                    g_count += 1
                elif nt == 'T':
                    t_count += 1
                i += 1
            elif bases[i] in ['+', '-']:
                sign = bases[i]
                i += 1
                length_str = ''
                while i < len(bases) and bases[i].isdigit():
                    length_str += bases[i]
                    i += 1
                length = int(length_str) if length_str else 0
                seq = bases[i:i+length].upper()
                i += length
                tokens.append(f'{sign}{length}{seq}')
            elif bases[i] in ['*', '^', '$']:
                if bases[i] == '^':
                    i += 2
                else:
                    i += 1
                tokens.append('SKIP')
            else:
                i += 1

        if not is_indel:
            for token in tokens:
                if token == 'REF' and ref_base == ref_allele:
                    ref_count += 1
                elif token == alt_allele:
                    alt_count += 1
        else:
            if len(ref_allele) > len(alt_allele):
                del_len = len(ref_allele) - len(alt_allele)
                del_seq = ref_allele[len(alt_allele):]
                del_pattern = f'-{del_len}{del_seq}'
                for token in tokens:
                    if token == 'REF':
                        ref_count += 1
                    elif token == del_pattern:
                        alt_count += 1
                # For multi-base deletions, check subsequent positions for *
                for check_pos in range(pos + 1, pos + del_len):
                    if (chrom, check_pos) in pileup_data:
                        bases2 = pileup_data[(chrom, check_pos)]['bases']
                        alt_count += sum(1 for b in bases2 if b == '*')
                        ref_count -= sum(1 for b in bases2 if b == '*')
                        if ref_count < 0:
                            ref_count = 0
            elif len(alt_allele) > len(ref_allele):
                ins_len = len(alt_allele) - len(ref_allele)
                ins_seq = alt_allele[len(ref_allele):]
                ins_pattern = f'+{ins_len}{ins_seq}'
                for token in tokens:
                    if token == 'REF':
                        ref_count += 1
                    elif token == ins_pattern:
                        alt_count += 1

        results.append({
            'chrom': chrom,
            'pos': pos,
            'ref': ref_allele,
            'alt': alt_allele,
            'ref_count': ref_count,
            'alt_count': alt_count,
            'A': a_count,
            'C': c_count,
            'G': g_count,
            'T': t_count
        })

    return results

def main():
    args = parse_arguments()
    variants = load_bed(args.bed)
    results = parse_pileup(args.input, variants)
    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False, sep='\t')

if __name__ == '__main__':
    main()