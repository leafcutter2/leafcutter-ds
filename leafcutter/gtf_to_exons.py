import argparse
import gzip
import re
import sys

import pandas as pd

parser = argparse.ArgumentParser(
    description="Extract exon annotations from a GTF file for use with leafcutter-ds --exon_file."
)
parser.add_argument("input_gtf", help="Input GTF file (plain or .gz).")
parser.add_argument("output_file", help="Output exon table (plain or .gz).")

args = parser.parse_args()

def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

sys.stderr.write(f"Reading {args.input_gtf}\n")

GTF_COLS = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

rows = []
with _open(args.input_gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        if parts[2] != "exon":
            continue
        rows.append(parts)

gtf = pd.DataFrame(rows, columns=GTF_COLS)

sys.stderr.write("Processing...\n")

def _extract_attr(attr_str, key):
    m = re.search(rf'{key} "([^"]+)"', attr_str)
    return m.group(1) if m else ""

gtf["gene_name"] = gtf["attributes"].apply(lambda x: _extract_attr(x, "gene_name"))

empty = gtf["gene_name"] == ""
if empty.any():
    sys.stderr.write("Warning: empty 'gene_name' attributes found; using 'gene_id' for those rows\n")
    gtf.loc[empty, "gene_name"] = gtf.loc[empty, "attributes"].apply(lambda x: _extract_attr(x, "gene_id"))

exons = gtf[["chr", "start", "end", "strand", "gene_name"]].drop_duplicates()

sys.stderr.write(f"Saving exons to {args.output_file}\n")
exons.to_csv(args.output_file, sep="\t", index=False)
