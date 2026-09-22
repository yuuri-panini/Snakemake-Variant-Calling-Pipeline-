"""
Extract the variants that fall inside a region of interest (e.g. one gene)
from a VCF file.

The original VCF header is copied across, so the output is still a valid
VCF with the real sample name, contig lines and INFO definitions, and can
be read by bedtools, bcftools, IGV, etc.

Usage (as in the course, plus the optional --chrom and -o):
    python parseVCF.py <input.vcf> <start> <end> [--chrom 17] [-o parsed_vcf_file.vcf]

start and end are 1-based and inclusive, like VCF positions. Leave out
--chrom to match the positions on every chromosome.
"""
import argparse
import gzip
import sys


def open_vcf(path):
    """Open a plain or gzipped VCF for reading as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def in_region(chrom, pos, args):
    """True if a variant at chrom:pos is inside the requested region."""
    if args.chrom is not None and chrom != args.chrom:
        return False
    if args.start is not None and pos < args.start:
        return False
    if args.end is not None and pos > args.end:
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="Keep the VCF variants inside a region.")
    parser.add_argument("vcf_file", help="input VCF (.vcf or .vcf.gz)")
    parser.add_argument("start", type=int, help="first position of the region")
    parser.add_argument("end", type=int, help="last position of the region")
    parser.add_argument("--chrom", help="chromosome/contig name exactly as in the VCF, e.g. 17 or chr17")
    parser.add_argument("-o", "--output", dest="output_file", default="parsed_vcf_file.vcf",
                        help="output VCF (default: parsed_vcf_file.vcf)")
    args = parser.parse_args()

    if args.start is not None and args.end is not None and args.start > args.end:
        sys.exit(f"--start ({args.start}) is after --end ({args.end})")

    kept = 0
    chroms_seen = set()

    with open_vcf(args.vcf_file) as infile, open(args.output_file, "w") as outfile:
        for line in infile:
            # Copy the header lines (## meta lines and the #CHROM line) unchanged
            if line.startswith("#"):
                outfile.write(line)
                continue

            # Variant lines: CHROM is column 0 and POS is column 1
            split_line = line.split("\t")
            chrom = split_line[0]
            pos = int(split_line[1])
            chroms_seen.add(chrom)

            if in_region(chrom, pos, args):
                outfile.write(line)
                kept += 1

    print(f"Kept {kept} variants from {args.vcf_file}", file=sys.stderr)

    # An empty result is often a chromosome naming mismatch (17 vs chr17)
    if args.chrom is not None and chroms_seen and args.chrom not in chroms_seen:
        print(
            f"Warning: no variants on chromosome '{args.chrom}'. "
            f"Chromosomes in this VCF: {', '.join(sorted(chroms_seen))}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
