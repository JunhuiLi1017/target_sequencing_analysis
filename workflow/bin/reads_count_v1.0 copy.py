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
                # Store as 1-based position for mpileup compatibility
                variants[(chrom, start + 1)] = {'ref': ref, 'alt': alt}
    return variants

def parse_pileup(pileup_file, variants):
    """Parse mpileup output and count ref/alt reads for SNVs and indels."""
    results = []
    
    with open(pileup_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue
            chrom = cols[0]
            pos = int(cols[1])  # 1-based
            ref_base = cols[2].upper()
            bases = cols[4]
            
            # Check if position is in the BED file
            if (chrom, pos) not in variants:
                continue
            
            ref_allele = variants[(chrom, pos)]['ref']
            alt_allele = variants[(chrom, pos)]['alt']
            
            # Initialize counts
            ref_count = 0
            alt_count = 0
            
            # Determine variant type (SNV or indel)
            is_indel = len(ref_allele) != len(alt_allele)
            
            if not is_indel:
                # SNV: Count single bases
                # Remove indels and special characters
                cleaned_bases = re.sub(r'[\+\-][0-9]+[ACGTNacgtn]+', '', bases)
                cleaned_bases = re.sub(r'[\*\^\$]', '', cleaned_bases)
                
                for base in cleaned_bases:
                    if base in ['.', ',']:  # Reference match
                        if ref_base == ref_allele:  # Ensure pileup ref matches BED ref
                            ref_count += 1
                    elif base.upper() == alt_allele:
                        alt_count += 1
            else:
                # Indel: Parse +N[sequence] and -N[sequence]
                # Reference count: Reads with no indel (.,,)
                cleaned_bases = bases
                indel_pattern = r'([\+\-])[0-9]+[ACGTNacgtn]+'
                indels = re.findall(indel_pattern, bases)
                
                # Count reference reads (no indel)
                ref_bases = re.sub(indel_pattern, '', cleaned_bases)
                ref_bases = re.sub(r'[\*\^\$]', '', ref_bases)
                ref_count = sum(1 for base in ref_bases if base in ['.', ','])
                
                # Count alternate reads (indels)
                if len(ref_allele) > len(alt_allele):
                    # Deletion: Look for -N[sequence]
                    del_len = len(ref_allele) - len(alt_allele)
                    del_seq = ref_allele[len(alt_allele):]
                    pattern = rf'-[{del_len}]{del_seq.upper()}'
                    alt_count = len(re.findall(pattern, bases, re.IGNORECASE))
                elif len(alt_allele) > len(ref_allele):
                    # Insertion: Look for +N[sequence]
                    ins_len = len(alt_allele) - len(ref_allele)
                    ins_seq = alt_allele[len(ref_allele):]
                    pattern = rf'\+[{ins_len}]{ins_seq.upper()}'
                    alt_count = len(re.findall(pattern, bases, re.IGNORECASE))
            
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