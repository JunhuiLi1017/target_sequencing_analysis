import argparse
import pandas as pd
import re

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Extract ref/alt read counts from mpileup for SNVs and indels.')
    parser.add_argument('--input', required=True, help='Input mpileup file from samtools mpileup.')
    parser.add_argument('--bed', required=True, help='BED file with variants (chrom, start, end, ref, alt).')
    parser.add_argument('--output', required=True, help='Output file for read counts (tab-separated).')
    return parser.parse_args()

def load_bed(bed_file):
    """Load BED file into a dictionary of variants."""
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
                # Store as 1-based position for mpileup, include end for multi-base deletions
                variants[(chrom, start + 1)] = {'ref': ref, 'alt': alt, 'end': int(cols[2])}
    return variants

def parse_pileup(pileup_file, variants):
    """Parse mpileup output and count ref/alt reads for SNVs and indels."""
    results = []
    
    # Load mpileup into memory for multi-position access
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
        end_pos = variant['end']  # 1-based end position
        
        ref_count = 0
        alt_count = 0
        
        # Determine variant type
        is_indel = len(ref_allele) != len(alt_allele)
        
        if not is_indel:
            # SNV
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
                    tokens.append(bases[i].upper())
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
            
            for token in tokens:
                if token == 'REF' and ref_base == ref_allele:
                    ref_count += 1
                elif token == alt_allele:
                    alt_count += 1
        else:
            # Indel
            if len(ref_allele) > len(alt_allele):
                # Deletion
                del_len = len(ref_allele) - len(alt_allele)
                del_seq = ref_allele[len(alt_allele):]
                
                # Check first position for -N[sequence]
                if (chrom, pos) not in pileup_data:
                    continue
                bases = pileup_data[(chrom, pos)]['bases']
                
                tokens = []
                i = 0
                while i < len(bases):
                    if bases[i] in ['.', ',']:
                        tokens.append('REF')
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
                
                # Count deletion at start position
                del_pattern = f'-{del_len}{del_seq}'
                for token in tokens:
                    if token == 'REF':
                        ref_count += 1
                    elif token == del_pattern:
                        alt_count += 1
                
                # For multi-base deletions, check subsequent positions for *
                for check_pos in range(pos + 1, pos + del_len):
                    if (chrom, check_pos) in pileup_data:
                        bases = pileup_data[(chrom, check_pos)]['bases']
                        # Count * as evidence of deletion continuation
                        alt_count += sum(1 for b in bases if b == '*')
                        # Adjust ref_count to exclude * reads
                        ref_count -= sum(1 for b in bases if b == '*')
                        if ref_count < 0:
                            ref_count = 0
            elif len(alt_allele) > len(ref_allele):
                # Insertion
                if (chrom, pos) not in pileup_data:
                    continue
                bases = pileup_data[(chrom, pos)]['bases']
                ins_len = len(alt_allele) - len(ref_allele)
                ins_seq = alt_allele[len(ref_allele):]
                
                tokens = []
                i = 0
                while i < len(bases):
                    if bases[i] in ['.', ',']:
                        tokens.append('REF')
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
            'alt_count': alt_count
        })
    
    return results

def main():
    """Main function to process variants and write output."""
    args = parse_arguments()
    variants = load_bed(args.bed)
    results = parse_pileup(args.input, variants)
    
    # Save results to a tab-separated file
    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False, sep='\t')

if __name__ == '__main__':
    main()