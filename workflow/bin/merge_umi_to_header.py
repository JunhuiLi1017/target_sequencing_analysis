script_name = "merge_umi_to_header"
Version = "1.0"

import gzip
import argparse
from Bio import SeqIO

def merge_umi_to_header(umi_fastq, in_fastq, out_fastq):
    # extract umi seq and its header into a dict as key-value
    umi_dict = {}
    with gzip.open(umi_fastq, "rt") as umi_handle:
        for record in SeqIO.parse(umi_handle, "fastq"):
            header_without_umi = record.id
            umi_dict[header_without_umi] = str(record.seq)

    with gzip.open(in_fastq, "rt") as in_handle, gzip.open(out_fastq, "wt") as out_handle:
        for record1 in SeqIO.parse(in_handle, "fastq"):
            header_without_umi1 = record1.id
            header_without_umi2 = record1.description.split(" ")[1]
            umi = umi_dict.get(header_without_umi1, "")
            record1.id = f"{header_without_umi1}:{umi}"
            record1.description = f"{header_without_umi2}"
            SeqIO.write(record1, out_handle, "fastq")

if __name__ == "__main__":
    """ merge umi from a fastq file to the header of R1 and/or R2 """
    user_input = argparse.ArgumentParser(description = "merge umi from a fastq file to the header of R1 and/or R2")
    user_input.add_argument("-u", "--umi_fastq", required = True, help="path to the UMI fastq file(gzipped)")
    user_input.add_argument("-i", "--in_fastq", required = True, help = "path to the input fastq file(gzipped)")
    user_input.add_argument("-o", "--out_fastq", required = True, help = "path to the output fastq file(gzipped)")
    args=user_input.parse_args()
    merge_umi_to_header(args.umi_fastq, args.in_fastq, args.out_fastq)
