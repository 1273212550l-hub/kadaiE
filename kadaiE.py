import gzip
import urllib.request
from mimetypes import guess_type
import Bio.PDB
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
import os
from Bio import SeqIO
if not os.path.exists("pdb_seqres.txt.gz"):
    urllib.request.urlretrieve(
        "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz",
        "pdb_seqres.txt.gz",
    )
fastafile = "pdb_seqres.txt.gz"
allrecords = 0
count = 0
with gzip.open(fastafile, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        allrecords += 1
        if len(record.seq) <= 100:
            count += 1
percentage = (count / allrecords) * 100
print(f"The rate of seq.length <= 100 : {percentage:.2f}%")

def _open(file):
    encoding = guess_type(file)[1]
    if encoding == "gzip":
        return gzip.open(file, mode="rt")
    else:
        return open(file)
urllib.request.urlretrieve("https://files.rcsb.org/download/1alk.cif.gz", "1alk.cif.gz")
ciffile = "1alk.cif.gz"
pdb_parser = FastMMCIFParser(QUIET=True)
with _open(ciffile) as handle:
    struc = pdb_parser.get_structure("1alk", handle)

atom_list = Bio.PDB.Selection.unfold_entities(struc, "A")
target_residue = struc[0]["A"][("H_PO4", 453, " ")]
ns = NeighborSearch(atom_list)
nearby_residues = set()
for atom in target_residue.get_atoms():
    nearby_residues = nearby_residues.union(ns.search(atom.get_coord(), 5, level="R"))
for res in nearby_residues:
    if res.get_resname() != "HOH":
        print(f"{res.get_resname()}-{res.id[1]}")