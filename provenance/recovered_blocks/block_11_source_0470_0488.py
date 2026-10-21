from pathlib import Path
import pandas as pd
from Bio import SeqIO, Phylo
base = Path.cwd()
tabla = pd.read_csv(base / "revision_30_noctuidos/Supplementary_ENC_GC3s_30.tsv", sep="\t")
elegidas = set(tabla["accession"].astype(str))
arbol = base / "04_alignments/concatenated_final40/PCG13_final40_by_gene.treefile"
puntas = {x.name for x in Phylo.read(arbol, "newick").get_terminals()}
print("Accesiones elegidas:", len(elegidas))
print("Presentes en árbol:", len(elegidas & puntas))
print("Faltan en árbol:", sorted(elegidas - puntas))
genes = "atp6 atp8 cox1 cox2 cox3 cytb nad1 nad2 nad3 nad4 nad4l nad5 nad6".split()
for gen in genes:
    ruta = base / f"04_alignments/gene_alignments_final40_clean/{gen}_clean_NT.fasta"
    secuencias = list(SeqIO.parse(ruta, "fasta"))
    ids = {s.id for s in secuencias}
    print(gen, "total:", len(secuencias),
          "elegidas:", len(ids & elegidas),
          "faltan:", sorted(elegidas - ids))
