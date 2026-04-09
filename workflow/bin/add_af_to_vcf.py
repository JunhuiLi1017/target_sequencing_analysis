#!/usr/bin/env python3
"""
add_af_to_vcf.py
Add AF field between DP and GQ in FORMAT and sample columns of a VCF.
"""

import sys
import gzip
import argparse
from collections import OrderedDict

def open_vcf(file_path):
    """Open .vcf or .vcf.gz"""
    return gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')

def write_vcf_line(out_file, fields):
    """Write a tab-separated line"""
    out_file.write('\t'.join(fields) + '\n')

def compute_af(ad_str):
    """Compute AF from AD (e.g., '0,5' → 5/(0+5) = 1.0)"""
    if not ad_str or ad_str in {'.', ''}:
        return '.'
    try:
        ref, alt = map(int, ad_str.split(','))
        total = ref + alt
        return '.' if total == 0 else f"{alt / total:.4f}"
    except:
        return '.'

def main():
    parser = argparse.ArgumentParser(description="Add AF field between DP and GQ in VCF")
    parser.add_argument("-i", "--input", required=True, help="Input VCF or VCF.GZ")
    parser.add_argument("-o", "--output", required=True, help="Output VCF.GZ")
    args = parser.parse_args()

    in_vcf = args.input
    out_vcf = args.output

    # Open input and output
    with open_vcf(in_vcf) as infile, \
         gzip.open(out_vcf, 'wt') as outfile:

        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('#'):
                # Write header unchanged (except last #CHROM line will be modified)
                if line.startswith('#CHROM'):
                    header_fields = line.split('\t')
                    format_idx = header_fields.index('FORMAT')
                    sample_start = format_idx + 1
                    # Modify FORMAT later per line
                    write_vcf_line(outfile, header_fields)
                else:
                    outfile.write(line + '\n')
                continue

            # --- Variant line ---
            fields = line.split('\t')
            format_str = fields[8]
            format_fields = format_str.split(':')

            # --- Insert AF after DP in FORMAT ---
            try:
                dp_idx = format_fields.index('DP')
            except ValueError:
                # If DP not present, append AF at end
                format_fields.append('AF')
                new_format = ':'.join(format_fields)
            else:
                # Insert after DP
                format_fields.insert(dp_idx + 1, 'AF')
                new_format = ':'.join(format_fields)

            fields[8] = new_format

            # --- Process each sample ---
            for i in range(9, len(fields)):
                geno_str = fields[i]
                if geno_str in {'.', ''}:
                    fields[i] = '.' * len(format_fields)
                    continue

                geno_fields = geno_str.split(':')
                # Map old indices to values
                old_dict = dict(zip(format_str.split(':'), geno_fields))

                # Compute AF from AD
                ad = old_dict.get('AD', '.')
                af = compute_af(ad)

                # Build new genotype string in new order
                new_geno = []
                for f in format_fields:
                    if f == 'AF':
                        new_geno.append(af)
                    else:
                        # Copy from old, preserve order
                        new_geno.append(old_dict.get(f, '.'))
                fields[i] = ':'.join(new_geno)

            # Write modified line
            write_vcf_line(outfile, fields)

    print(f"Done! Output written to: {out_vcf}")

if __name__ == "__main__":
    main()
