#!/usr/bin/env python3
"""Reconstruct final alignment from archived cleaned gene alignments.
This does not reconstruct the earlier MACSE alignment/cleaning operation.
"""
from pathlib import Path
import json
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
GENES = 'atp6 atp8 cox1 cox2 cox3 cytb nad1 nad2 nad3 nad4 nad4l nad5 nad6'.split()
def canon(acc):
    return acc.rsplit('.',1)[0] if acc.rsplit('.',1)[-1].isdigit() else acc
root=Path(__file__).resolve().parents[1]
out=root/'results/phylogeny';out.mkdir(parents=True,exist_ok=True)
manifest=pd.read_csv(root/'data/accessions_30.tsv',sep='\t')
selected={canon(a) for a in manifest.accession}|{'MN793322'}
parts={};ids=None;start=1;partition=['#nexus','begin sets;']
for gene in GENES:
    records={r.id:str(r.seq).upper() for r in SeqIO.parse(root/f'data/clean_gene_alignments/{gene}_clean_NT.fasta','fasta') if canon(r.id) in selected}
    assert len(records)==31 and len({len(s) for s in records.values()})==1
    length=len(next(iter(records.values())))
    assert length%3==0 and all('!' not in s for s in records.values())
    if ids is None:ids=set(records)
    assert set(records)==ids
    parts[gene]=records
    partition.append(f'charset {gene} = {start}-{start+length-1};');start+=length
alignment=[SeqRecord(Seq(''.join(parts[g][a] for g in GENES)),id=a,description='') for a in sorted(ids)]
assert start-1==11229
SeqIO.write(alignment,out/'PCG13_31_aligned.fasta','fasta')
partition.append('end;');(out/'PCG13_31_partitions.nex').write_text('\n'.join(partition)+'\n')
old={r.id:str(r.seq) for r in SeqIO.parse(root/'data/reference/phylogeny_31/PCG13_31_aligned.fasta','fasta')}
assert {r.id:str(r.seq) for r in alignment}==old
# Derive a fresh 2-million-generation file from the archived initial MrBayes file.
# Only ngen, append flag and output filename are changed; models/partitions remain.
import re
text=(root/'data/reference/phylogeny_31/PCG13_31_BI.nex').read_text()
text=re.sub(r'ngen=\d+','ngen=2000000',text)
text=text.replace('append=yes','append=no').replace('filename=PCG13_31_BI','filename=BI_new')
(out/'BI_fresh.nex').write_text(text)
(out/'alignment_validation.json').write_text(json.dumps(dict(sequences=31,positions=11229,archived_alignment_identical=True),indent=2))
print('Verified 31 sequences x 11229 positions; fresh BI configuration prepared.')
