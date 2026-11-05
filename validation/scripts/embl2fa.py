#!/usr/bin/env python3
"""Convert a FamDB EMBL export to FASTA keyed by family name (NM field).

Both createFakeSequence.pl (--tecons_file) and rusted-garlic emit (--te_cons) then use
the same TE consensus set. Needed because GARLIC 1.4's EMBL parser keys records by the
ID field, which in Dfam 4.0 is the accession (DF...), not the family name. When a name
appears twice (RepBase-derived stub + Dfam record), the Dfam (DF) record wins.

Usage:
  famdb.py -i <Dfam 4.0 FamDB> families -ad --curated --format embl 'homo sapiens' > lib.embl
  python3 embl2fa.py lib.embl > dfam-4.0-homo-sapiens.fa
"""
import re
import sys

recs = {}  # name -> (is_dfam, seq); insertion order preserved
acc = name = seq = None


def flush():
    if name and seq:
        is_dfam = acc.startswith("DF")
        if name not in recs or (is_dfam and not recs[name][0]):
            recs[name] = (is_dfam, seq)


for line in open(sys.argv[1]):
    if line.startswith("ID "):
        acc = line.split()[1].rstrip(";")
        name = acc.rstrip("_")
        seq = ""
    elif line.startswith("NM "):
        name = line.split()[1]
    elif line.startswith("     ") and seq is not None:
        seq += re.sub(r"[\s\d]", "", line)
    elif line.startswith("//"):
        flush()
        acc = name = seq = None
flush()

for n, (_, s) in recs.items():
    sys.stdout.write(f">{n}\n{s}\n")
