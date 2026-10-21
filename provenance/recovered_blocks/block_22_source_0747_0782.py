from pathlib import Path
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
datos = pd.read_csv("../Supplementary_ENC_GC3s_30.tsv", sep="\t")
def canon(s):
    return str(s).split(".")[0]
nombres = dict(zip(datos.accession.map(canon), datos.species))
nombres["MN793322"] = "Ostrinia nubilalis"
arbol = Phylo.read("PCG13_31_ML.treefile", "newick")
assert len(arbol.get_terminals()) == 31
def etiqueta(c):
    if not c.is_terminal():
        return None
    especie = nombres[canon(c.name)].replace(" ", r"\ ")
    return rf"$\it{{{especie}}}$  {c.name}"
def apoyo(c):
    if c.is_terminal():
        return None
    return c.name or (f"{c.confidence:g}" if c.confidence is not None else None)
fig, ax = plt.subplots(figsize=(15, 12))
Phylo.draw(arbol, axes=ax, do_show=False,
           label_func=etiqueta, branch_labels=apoyo)
ax.set_xlabel("Substitutions per site")
ax.set_ylabel("")
ax.set_yticks([])
fig.tight_layout()
destino = Path("/mnt/e/Mitogenoma/nuevas figuras")
destino.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    archivo = Path(f"Figure1_phylogeny_ML_31.{ext}")
    fig.savefig(archivo, dpi=300, bbox_inches="tight")
    shutil.copy2(archivo, destino / archivo.name)
    print("Copiado:", destino / archivo.name)
plt.close(fig)
