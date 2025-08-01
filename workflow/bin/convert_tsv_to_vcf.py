import sys

def convert_tsv_to_vcf(input_file_path, output_file_path):
    """
    Converts a custom tab-delimited variant file to VCF format.

    Args:
        input_file_path (str): Path to the input tab-delimited file.
        output_file_path (str): Path to the output VCF file.
    """

    # Define the VCF Header lines, including INFO and FILTER field definitions.
    # These describe the custom fields from your input that will go into the VCF INFO column.
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
        "##INFO=<ID=DP_MQ,Number=1,Type=Integer,Description=\"Depth with mapping quality filter\">",
        "##INFO=<ID=BCNT,Number=1,Type=Integer,Description=\"Base Count\">",
        "##INFO=<ID=ENS,Number=1,Type=Float,Description=\"Estimated Non-Silent Event Count\">",
        "##INFO=<ID=EBAF,Number=1,Type=Float,Description=\"Estimated B-Allele Frequency\">",
        "##INFO=<ID=LMP,Number=1,Type=Float,Description=\"Log Mutation Probability\">",
        "##INFO=<ID=LEP,Number=1,Type=Float,Description=\"Log Error Probability\">",
        "##INFO=<ID=LMQF,Number=1,Type=Float,Description=\"Low Mapping Quality Fraction\">",
        "##INFO=<ID=PGF,Number=1,Type=Float,Description=\"Proximity Gap Fraction\">",
        "##INFO=<ID=PGF_MUT,Number=1,Type=Float,Description=\"Proximity Gap Fraction (Mut)\">",
        "##INFO=<ID=MNMF,Number=1,Type=Float,Description=\"Mismatch Fraction (Non-Mappable)\">",
        "##INFO=<ID=MMMQ,Number=1,Type=Integer,Description=\"Mean Mismatch Mapping Quality\">",
        "##INFO=<ID=FOXOG,Number=1,Type=Float,Description=\"FoxoG score\">",
        "##INFO=<ID=LR,Number=1,Type=Float,Description=\"Likelihood Ratio\">",
        # Define FILTER fields based on your input's 'Filter' column values
        "##FILTER=<ID=normalFilter,Description=\"Passed normal filtering criteria (original: normalFilter;)\">",
        "##FILTER=<ID=background_error,Description=\"Filtered due to background error\">",
        "##FILTER=<ID=proximityGap_rep0,Description=\"Filtered due to proximity gap (replicate 0)\">",
        "##FILTER=<ID=nmFilter,Description=\"Filtered due to NM (Normalized Mismatch) criteria\">",
        "##FILTER=<ID=low_nonSECnt,Description=\"Filtered due to low non-silent event count\">",
        # The main VCF column header line
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ]

    # Map your input column names (keys) to standard VCF INFO field IDs (values)
    # This dictionary helps in constructing the INFO field for each variant.
    info_field_map = {
        "depth": "DP",
        "depth_mQ": "DP_MQ",
        "bCnt": "BCNT",
        "estimatedNonSECnt": "ENS",
        "estimatedBAF": "EBAF",
        "logMutProb": "LMP",
        "logErrProb": "LEP",
        "lowMQfrac": "LMQF",
        "proximityGapFrac": "PGF",
        "proximityGapFrac_mut": "PGF_MUT",
        "misNMfrac": "MNMF",
        "meanMisMQ": "MMMQ",
        "FoxoG": "FOXOG",
        "likelihoodRatio": "LR"
    }

    try:
        with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
            # Print the VCF header to the output file
            for line in vcf_header:
                outfile.write(line + '\n')

            # Read the header line from the input file
            header_line = infile.readline().strip()
            # Clean up the header line: remove '#' and any leading/trailing spaces, then split by tab
            header_cols = [col.strip() for col in header_line.lstrip('#').split('\t')]
            
            # Create a dictionary to easily get column index by name
            col_indices = {col: i for i, col in enumerate(header_cols)}

            # Validate that all essential VCF columns are present in the input
            essential_cols = ["#chr", "pos", "ID", "ref", "alt", "Filter"]
            for col in essential_cols:
                # Check for both '#chr' and 'chr' as column names
                if col not in col_indices and col.lstrip('#') not in col_indices:
                    raise ValueError(f"Missing essential column in input file: {col}")
            
            # Adjust the column name for 'chr' if it was '#chr' in the input header
            if '#chr' in col_indices:
                col_indices['chr'] = col_indices.pop('#chr')
            elif 'chr' in col_indices: # Ensure 'chr' is the key if it was already 'chr'
                pass
            else:
                raise ValueError("Could not find 'chr' or '#chr' column in input file header.")

            # Process each data line in the input file
            for line in infile:
                line = line.strip()
                # Skip empty lines or lines that start with '#' (comments or repeated header)
                if not line or line.startswith('#'):
                    continue

                # Split the line into parts based on tab delimiter and strip whitespace
                parts = [p.strip() for p in line.split('\t')]

                # Extract the mandatory VCF fields using their identified column indices
                chrom = parts[col_indices['chr']]
                pos = parts[col_indices['pos']]
                variant_id = parts[col_indices['ID']]
                ref = parts[col_indices['ref']]
                alt = parts[col_indices['alt']]
                
                # The QUAL (Quality) field is not explicitly in your input, so set it to '.' (missing)
                qual = "." 

                # Process the FILTER field
                raw_filter = parts[col_indices['Filter']].strip()
                # If the raw filter is "PASS" or "normalFilter;", map it to "PASS" in VCF.
                # Otherwise, use the raw filter string. Multiple filters should be separated by ';'
                if raw_filter == "PASS" or raw_filter == "normalFilter;":
                    vcf_filter = "PASS"
                else:
                    # Remove any trailing semicolon and ensure multiple filters are ';' delimited
                    vcf_filter = raw_filter.rstrip(';')

                # Construct the INFO field by iterating through the info_field_map
                info_parts = []
                for col_name, info_id in info_field_map.items():
                    # Check if the column exists in the current line's parts
                    if col_name in col_indices:
                        value = parts[col_indices[col_name]]
                        # Replace '-' (used for missing data in your input) with '.' (VCF standard for missing)
                        if value == "-":
                            value = "."
                        info_parts.append(f"{info_id}={value}")
                
                # Join all INFO key=value pairs with a semicolon
                vcf_info = ";".join(info_parts)

                # Write the complete VCF record to the output file
                outfile.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{vcf_filter}\t{vcf_info}\n")

    except FileNotFoundError:
        print(f"Error: One of the files not found. Input: '{input_file_path}', Output: '{output_file_path}'", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)
    except IndexError:
        print(f"Error: Malformed line or incorrect column count in file. Please check line: {line}", file=sys.stderr)
        sys.exit(1)

# This block allows the script to be run from the command line.
# Example usage:
# Save the code above as a Python file (e.g., `tsv_to_vcf.py`).
# Then, run it from your terminal like this:
# python tsv_to_vcf.py <input_tsv_file> <output_vcf_file>
# Replace `<input_tsv_file>` with the actual path to your tab-delimited file.
# Replace `<output_vcf_file>` with the desired name for your VCF output file.
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python tsv_to_vcf.py <input_tsv_file> <output_vcf_file>", file=sys.stderr)
        sys.exit(1)
    
    input_tsv_file = sys.argv[1]
    output_vcf_file = sys.argv[2]
    convert_tsv_to_vcf(input_tsv_file, output_vcf_file)