from pathlib import Path
from Bio import SeqIO
origen = Path("~/noctuidae_mitogenomes_v2/04_alignments/gene_alignments_final40_clean").expanduser()
genes = "atp6 atp8 cox1 cox2 cox3 cytb nad1 nad2 nad3 nad4 nad4l nad5 nad6".split()
inicio = 1
lineas = ["#nexus", "begin sets;"]
for gen in genes:
    largo = len(next(SeqIO.parse(origen / f"{gen}_clean_NT.fasta", "fasta")).seq)
    fin = inicio + largo - 1
    lineas.append(f"    charset {gen} = {inicio}-{fin};")
    inicio = fin + 1
assert inicio - 1 == 11229
lineas.append("end;")
Path("PCG13_31_partitions.nex").write_text("\n".join(lineas) + "\n")
print("Particiones:", len(genes), "| Longitud total:", inicio - 1)
