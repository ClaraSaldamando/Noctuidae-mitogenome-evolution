from Bio import SeqIO
import pandas as pd
ruta = "04_alignments/concatenated_final40/PCG13_concatenated_40mitogenomes.fasta"
r = list(SeqIO.parse(ruta, "fasta"))
print("Concatenado original:", len(r), "secuencias")
print("Longitudes:", sorted({len(x.seq) for x in r}))
print("Mapa taxonómico:")
m = pd.read_csv("06_codon_usage_final40/accession_species_final40.csv")
print(m[m.astype(str).apply(
    lambda c: c.str.contains("Ostrinia", case=False, na=False)
).any(axis=1)].to_string(index=False))
