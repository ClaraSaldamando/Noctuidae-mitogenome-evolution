from Bio import SeqIO
from pathlib import Path
carpeta = Path("04_alignments/gene_alignments_final40_clean")
genes = "atp6 atp8 cox1 cox2 cox3 cytb nad1 nad2 nad3 nad4 nad4l nad5 nad6".split()
total = 0
for gen in genes:
    r = list(SeqIO.parse(carpeta / f"{gen}_clean_NT.fasta", "fasta"))
    largo = len(r[0].seq)
    total += largo
    externo = any(x.id == "MN793322.1" for x in r)
    print(gen, largo, "bp; Ostrinia:", externo)
print("Longitud concatenada:", total, "bp")
