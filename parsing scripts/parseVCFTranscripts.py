"""
Report the missense variants inside a region of interest, one row per
affected transcript, from a snpEff-annotated VCF.

snpEff puts its annotations in the ANN field of INFO. Each annotation
(one per allele x transcript) is separated by a comma, and its fields are
separated by a pipe '|':

    Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID |
    Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c |
    HGVS.p | cDNA.pos/cDNA.length | CDS.pos/CDS.length | AA.pos/AA.length |
    Distance | ERRORS/WARNINGS/INFO

Usage (as in the course, plus the optional --chrom):
    python parseVCFTranscripts.py <annotated.vcf> <start> <end> [--chrom 17] > report.tsv

The report is tab-separated and printed to stdout.
"""
import argparse
import gzip
import sys

# Positions of the fields we need inside one ANN annotation
ANNOTATION = 1
GENE_NAME = 3
FEATURE_ID = 6  # the transcript ID for transcript annotations
HGVS_C = 9      # nucleotide change
HGVS_P = 10     # amino acid change

HEADER = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER",
          "ALLELE", "GENE", "TRANSCRIPT", "NUC_CHANGE", "AA_CHANGE"]


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


def ann_annotations(info):
    """Yield each snpEff annotation in an INFO column as a list of fields.

    The INFO column is split into its ';' separated keys first, so other
    INFO fields (DP, AF, LOF, NMD, ...) can't be mistaken for ANN fields.
    """
    for field in info.split(";"):
        if field.startswith("ANN="):
            for annotation in field[len("ANN="):].split(","):
                yield annotation.split("|")


def main():
    parser = argparse.ArgumentParser(description="List missense variants per transcript.")
    parser.add_argument("vcf_file", help="snpEff-annotated VCF (.vcf or .vcf.gz)")
    parser.add_argument("start", type=int, help="first position of the region")
    parser.add_argument("end", type=int, help="last position of the region")
    parser.add_argument("--chrom", help="chromosome/contig name exactly as in the VCF, e.g. 17 or chr17")
    args = parser.parse_args()

    if args.start is not None and args.end is not None and args.start > args.end:
        sys.exit(f"--start ({args.start}) is after --end ({args.end})")

    print("\t".join(HEADER))

    rows = 0
    variants_in_region = 0
    variants_without_ann = 0
    chroms_seen = set()

    with open_vcf(args.vcf_file) as infile:
        for line in infile:
            # Only interested in the variant lines, not the header lines
            if line.startswith("#"):
                continue

            split_line = line.rstrip("\n").split("\t")
            chrom, pos, var_id, ref, alt, _qual, filt, info = split_line[:8]
            chroms_seen.add(chrom)

            if not in_region(chrom, int(pos), args):
                continue
            variants_in_region += 1

            annotations = list(ann_annotations(info))
            if not annotations:
                variants_without_ann += 1

            for ann in annotations:
                if len(ann) <= HGVS_P:
                    continue  # malformed annotation
                # A variant can have several effects on one transcript joined
                # by '&', e.g. missense_variant&splice_region_variant
                effects = ann[ANNOTATION].split("&")
                if "missense_variant" in effects:
                    print("\t".join([
                        chrom, pos, var_id, ref, alt, filt,
                        ann[0], ann[GENE_NAME], ann[FEATURE_ID], ann[HGVS_C], ann[HGVS_P],
                    ]))
                    rows += 1

    print(f"{variants_in_region} variants in region, {rows} missense transcript rows",
          file=sys.stderr)

    if variants_in_region and variants_without_ann == variants_in_region:
        print("Warning: no ANN field found. Has this VCF been annotated with snpEff?",
              file=sys.stderr)
    if args.chrom is not None and chroms_seen and args.chrom not in chroms_seen:
        print(
            f"Warning: no variants on chromosome '{args.chrom}'. "
            f"Chromosomes in this VCF: {', '.join(sorted(chroms_seen))}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
