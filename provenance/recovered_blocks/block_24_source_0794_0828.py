from pathlib import Path
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
d = pd.read_csv("../Supplementary_ENC_GC3s_30.tsv", sep="\t")
nombres = dict(zip(d.accession.astype(str).str.split(".").str[0], d.species))
nombres["MN793322"] = "Ostrinia nubilalis"
t = Phylo.read("PCG13_31_ML.treefile", "newick")
t.root_with_outgroup("MN793322.1")
print("Descendientes de la raíz:", len(t.root.clades))
print("Politomías internas:", sum(len(c.clades) > 2 for c in t.get_nonterminals()))
assert all(len(c.clades) <= 2 for c in t.get_nonterminals())
for c in t.find_clades():
    if c is not t.root:
        c.branch_length = 1
def etiqueta(c):
    if not c.is_terminal():
        return None
    especie = nombres[c.name.split(".")[0]].replace(" ", r"\ ")
    return rf"$\it{{{especie}}}$  {c.name}"
fig, ax = plt.subplots(figsize=(15, 13))
Phylo.draw(t, axes=ax, do_show=False, label_func=etiqueta,
           branch_labels=lambda c: c.name if not c.is_terminal() else None)
ax.set_xlabel("Topology shown as a cladogram; branch lengths are not to scale")
ax.set_ylabel("")
ax.set_yticks([])
fig.tight_layout()
destino = Path("/mnt/e/Mitogenoma/nuevas figuras")
for ext in ("png", "pdf"):
    archivo = Path(f"Figure1_phylogeny_ML_31_cladogram.{ext}")
    fig.savefig(archivo, dpi=300, bbox_inches="tight")
    shutil.copy2(archivo, destino / archivo.name)
    print("Copiado:", destino / archivo.name)
plt.close(fig)
