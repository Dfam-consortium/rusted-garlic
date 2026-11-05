"""Build a FASTA of test windows for comparing complexity implementations.

Usage: python3 make_windows.py chr22.fa simulated.fa > windows.fa
Each record is one window: 150 real chr22 windows (some with N), 40 simulated
windows (upper-cased), and hand-built extremes / short records.
"""
import random, sys
def fasta(path):
    out=[]; buf=[]
    for l in open(path):
        if l.startswith('>'):
            if buf: out.append(''.join(buf)); buf=[]
        else: buf.append(l.strip())
    if buf: out.append(''.join(buf))
    return out
rng=random.Random(1)
recs=[]
ref=''.join(fasta(sys.argv[1])).upper()
clean=withn=0
while clean<140 or withn<10:
    p=rng.randrange(0,len(ref)-1000); w=ref[p:p+1000]; n=w.count('N')
    if n==0 and clean<140: recs.append((f"chr22_{p}",w)); clean+=1
    elif 0<n<=300 and withn<10: recs.append((f"chr22N_{p}",w)); withn+=1
sim=''.join(fasta(sys.argv[2])).upper()
for i in range(40):
    p=rng.randrange(0,len(sim)-1000); recs.append((f"sim_{p}",sim[p:p+1000]))
recs+= [("polyA","A"*1000),("CA_repeat","CA"*500),("ACGT_repeat","ACGT"*250),
        ("CAG_repeat",("CAG"*334)[:1000]),("A999C","A"*999+"C"),
        ("quarterN",''.join(rng.choice("ACGT") for _ in range(750))+"N"*250),
        ("short50",''.join(rng.choice("ACGT") for _ in range(50))),
        ("short7","ACGTACG"),("one","A")]
for name,s in recs:
    print(">"+name)
    for i in range(0,len(s),60): print(s[i:i+60])
