#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Extract ref/alt read counts and nucleotide counts from mpileup for variants. "
            "Also outputs strand-specific ref/alt counts and Phred qualities for ALT-supporting reads."
        )
    )
    parser.add_argument("--input", required=True, help="Input mpileup file from samtools mpileup (plain text).")
    parser.add_argument("--bed", required=True, help="BED file with variants (chrom, start, end, ref, alt).")
    parser.add_argument("--output", required=True, help="Output file (tab-separated).")
    return parser.parse_args()


def load_bed(bed_file):
    """
    Expect at least 5 columns:
      chrom  start(0-based)  end  ref  alt
    """
    variants = {}
    with open(bed_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom = cols[0]
            start = int(cols[1])  # 0-based
            end = int(cols[2])
            ref = cols[3].upper()
            alt = cols[4].upper()
            variants[(chrom, start + 1)] = {"ref": ref, "alt": alt, "end": end}
    return variants


def phred(ch: str) -> int:
    return ord(ch) - 33


def parse_bases_and_quals_with_strand(bases: str, quals: str):
    """
    Parse mpileup bases string aligned to qualities, tracking strand.

    Returns list of (call, qscore, strand) where:
      call in {'REF','A','C','G','T','N','DEL','REFSKIP'}
      strand in {'+','-'} when applicable; otherwise None

    mpileup rules:
      - '^' + next char: start of read + mapping quality (no base, no qual consumed)
      - '$': end of read (no base, no qual consumed)
      - '+<len><seq>' / '-<len><seq>': indel annotation (no base, no qual consumed)
      - '.', ',': ref match ('.' forward, ',' reverse) consumes 1 qual
      - A/C/G/T: forward; a/c/g/t: reverse consumes 1 qual
      - '*': deletion placeholder consumes 1 qual (strand not represented -> None)
      - '>'/'<': ref skip consumes 1 qual (strand not represented -> None)
      - N/n: ambiguous base (N forward, n reverse) consumes 1 qual
    """
    obs = []
    i = 0
    q = 0

    while i < len(bases):
        b = bases[i]

        # read segment start marker
        if b == "^":
            i += 2
            continue

        # read segment end marker
        if b == "$":
            i += 1
            continue

        # indel annotation; does not consume quality
        if b in ["+", "-"]:
            i += 1
            nstr = ""
            while i < len(bases) and bases[i].isdigit():
                nstr += bases[i]
                i += 1
            n = int(nstr) if nstr else 0
            i += n
            continue

        # remaining symbols consume one base-quality
        if q >= len(quals):
            break
        qscore = phred(quals[q])
        q += 1

        # ref match
        if b == ".":
            obs.append(("REF", qscore, "+"))
            i += 1
            continue
        if b == ",":
            obs.append(("REF", qscore, "-"))
            i += 1
            continue

        # base calls (strand from case)
        if b in ["A", "C", "G", "T"]:
            obs.append((b, qscore, "+"))
            i += 1
            continue
        if b in ["a", "c", "g", "t"]:
            obs.append((b.upper(), qscore, "-"))
            i += 1
            continue

        # ambiguous base
        if b == "N":
            obs.append(("N", qscore, "+"))
            i += 1
            continue
        if b == "n":
            obs.append(("N", qscore, "-"))
            i += 1
            continue

        # deletion placeholder / refskip (strand not encoded)
        if b == "*":
            obs.append(("DEL", qscore, None))
            i += 1
            continue
        if b in [">", "<"]:
            obs.append(("REFSKIP", qscore, None))
            i += 1
            continue

        # unknown
        i += 1

    return obs


def count_indel_support(bases: str, ref_allele: str, alt_allele: str, depth_mpileup: int):
    """
    Count reads supporting an indel (insertion/deletion) from the mpileup bases string.

    We assume a standard VCF-style representation:
      - Insertion: len(alt) > len(ref), alt = ref + inserted_seq
      - Deletion:  len(ref) > len(alt), ref = alt + deleted_seq

    mpileup encodes indels at a position as:
      +<len><seq>  insertion of <seq>
      -<len><seq>  deletion of <seq>

    We:
      - count occurrences of the exact +/− pattern matching this variant as alt-support
      - approximate ref_count as depth_mpileup - alt_count
    """
    alt_count = 0

    # Determine type and the sequence that should appear in mpileup
    if len(alt_allele) > len(ref_allele):
        # insertion
        ins_seq = alt_allele[1:]  # everything after the first ref base
        ins_len = len(ins_seq)
        target_sign = "+"
        target_len = ins_len
        target_seq = ins_seq.upper()
    elif len(ref_allele) > len(alt_allele):
        # deletion
        del_seq = ref_allele[1:]  # deleted sequence after leftmost base
        del_len = len(del_seq)
        target_sign = "-"
        target_len = del_len
        target_seq = del_seq.upper()
    else:
        # Not an indel (lengths equal) -> nothing to do here
        return 0, 0

    i = 0
    n_bases = len(bases)
    while i < n_bases:
        b = bases[i]

        # Skip read start/end markers
        if b == "^":
            i += 2
            continue
        if b == "$":
            i += 1
            continue

        # Look for the specific +lenSEQ / -lenSEQ pattern
        if b in ["+", "-"]:
            sign = b
            i += 1
            # parse length
            len_str = ""
            while i < n_bases and bases[i].isdigit():
                len_str += bases[i]
                i += 1
            if not len_str:
                continue
            try:
                indel_len = int(len_str)
            except ValueError:
                continue

            # extract sequence of that length
            seq = bases[i : i + indel_len]
            i += indel_len

            if (
                sign == target_sign
                and indel_len == target_len
                and seq.upper() == target_seq
            ):
                alt_count += 1
            continue

        # All other symbols: advance by one
        i += 1

    # Approximate ref-supporting reads as reads without this exact indel
    ref_count = max(depth_mpileup - alt_count, 0)
    return ref_count, alt_count


def parse_pileup(pileup_file, variants):
    # Load mpileup lines into dict for fast lookup
    pileup_data = {}
    with open(pileup_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            # expected: chrom pos ref depth bases quals
            if len(cols) < 6:
                continue
            chrom = cols[0]
            pos = int(cols[1])
            pileup_data[(chrom, pos)] = {
                "ref_base": cols[2].upper(),
                "depth": int(cols[3]),
                "bases": cols[4],
                "quals": cols[5],
            }

    results = []

    for (chrom, pos), variant in variants.items():
        if (chrom, pos) not in pileup_data:
            continue

        ref_allele = variant["ref"]
        alt_allele = variant["alt"]

        ref_base = pileup_data[(chrom, pos)]["ref_base"]
        depth_mpileup = pileup_data[(chrom, pos)]["depth"]
        bases = pileup_data[(chrom, pos)]["bases"]
        quals = pileup_data[(chrom, pos)]["quals"]

        is_indel = len(ref_allele) != len(alt_allele)

        obs = parse_bases_and_quals_with_strand(bases, quals)

        # nucleotide totals INCLUDING ref matches '.'/','
        A = C = G = T = 0
        Nn = DEL = REFSKIP = 0

        # ref/alt counts + strand breakdown
        ref_count = alt_count = 0
        ref_fwd = ref_rev = 0
        alt_fwd = alt_rev = 0

        alt_qscores = []
        alt_qscores_fwd = []
        alt_qscores_rev = []

        for call, qscore, strand in obs:
            # ---- nucleotide totals ----
            if call == "REF":
                if ref_base == "A":
                    A += 1
                elif ref_base == "C":
                    C += 1
                elif ref_base == "G":
                    G += 1
                elif ref_base == "T":
                    T += 1
            elif call == "A":
                A += 1
            elif call == "C":
                C += 1
            elif call == "G":
                G += 1
            elif call == "T":
                T += 1
            elif call == "N":
                Nn += 1
            elif call == "DEL":
                DEL += 1
            elif call == "REFSKIP":
                REFSKIP += 1

            # ---- ref/alt counts + strand-specific counts (SNV only) ----
            if not is_indel:
                if call == "REF" and ref_base == ref_allele:
                    ref_count += 1
                    if strand == "+":
                        ref_fwd += 1
                    elif strand == "-":
                        ref_rev += 1
                elif call == alt_allele:
                    alt_count += 1
                    alt_qscores.append(qscore)
                    if strand == "+":
                        alt_fwd += 1
                        alt_qscores_fwd.append(qscore)
                    elif strand == "-":
                        alt_rev += 1
                        alt_qscores_rev.append(qscore)

        # For indels, override ref/alt counts using the explicit indel parser
        if is_indel:
            ref_count, alt_count = count_indel_support(
                bases, ref_allele, alt_allele, depth_mpileup
            )
            # We do not currently derive strand/quality metrics for indels,
            # so ref_fwd/ref_rev/alt_fwd/alt_rev and alt_qscores* remain 0/empty.

        effective_depth = A + C + G + T
        alt_qual_mean = (sum(alt_qscores) / len(alt_qscores)) if alt_qscores else ""
        alt_qual_mean_fwd = (sum(alt_qscores_fwd) / len(alt_qscores_fwd)) if alt_qscores_fwd else ""
        alt_qual_mean_rev = (sum(alt_qscores_rev) / len(alt_qscores_rev)) if alt_qscores_rev else ""

        results.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref_allele,
                "alt": alt_allele,
                "depth_mpileup": depth_mpileup,
                "effective_depth": effective_depth,
                "ref_count": ref_count,
                "alt_count": alt_count,
                "ref_fwd": ref_fwd,
                "ref_rev": ref_rev,
                "alt_fwd": alt_fwd,
                "alt_rev": alt_rev,
                "A": A,
                "C": C,
                "G": G,
                "T": T,
                "N": Nn,
                "DEL": DEL,
                "REFSKIP": REFSKIP,
                "alt_quals": ",".join(map(str, alt_qscores)),
                "alt_qual_mean": alt_qual_mean,
                "alt_qual_mean_fwd": alt_qual_mean_fwd,
                "alt_qual_mean_rev": alt_qual_mean_rev,
            }
        )

    return results


def main():
    args = parse_arguments()
    variants = load_bed(args.bed)
    results = parse_pileup(args.input, variants)
    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()

