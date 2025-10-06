import argparse
import pandas as pd
import gzip
from Bio import SeqIO
import urllib.request
import os
def main():
    parser = argparse.ArgumentParser(description="Calculate short sequences from PDB fasta")
    parser.add_argument("-i", "--input", required=True, help="Input PDB fasta file (.gz)")
    parser.add_argument("-l", "--length", type=int, default=100, help="Maximum length to define short sequences")
    args = parser.parse_args()
    fastafile = args.input
    if not os.path.exists(fastafile):
        print(f"{fastafile} not found, downloading from PDB...")
        urllib.request.urlretrieve(
            "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz",
            fastafile
        )
    allrecords = 0
    count = 0
    with gzip.open(fastafile, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            allrecords += 1
            if len(record.seq) <= args.length:
                count += 1
    percentage = (count / allrecords) * 100 if allrecords > 0 else 0
    print(f"The rate of seq.length <= {args.length} : {percentage:.2f}%")            