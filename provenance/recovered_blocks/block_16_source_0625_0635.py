from pathlib import Path
from Bio import SeqIO
for p in Path(".").rglob("*"):
    if p.suffix.lower() not in {".fasta", ".fa", ".fas", ".fna"}:
        continue
    try:
        registros = list(SeqIO.parse(p, "fasta"))
    except Exception:
        continue
    if len(registros) == 31 and len({len(r.seq) for r in registros}) == 1:
        print(p.resolve(), "longitud:", len(registros[0].seq))
