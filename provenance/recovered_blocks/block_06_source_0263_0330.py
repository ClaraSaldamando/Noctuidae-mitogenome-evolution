from pathlib import Path
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
carpeta = Path.cwd()
tax = pd.read_csv("Supplementary_ENC_GC3s_30.tsv", sep="\t")
especie = dict(zip(tax["accession"].astype(str), tax["species"]))
filas = []
for rec in SeqIO.parse("codonw_input_30.fasta", "fasta"):
    acc = rec.id
    if acc not in especie:
        raise ValueError(f"Accesión no encontrada en la tabla: {acc}")
    seq = str(rec.seq).upper()
    if len(seq) % 3:
        raise ValueError(f"Secuencia fuera de marco: {acc}")
    p1, p2, p3 = seq[0::3], seq[1::3], seq[2::3]
    if set(seq) - set("ATGC"):
        raise ValueError(f"Bases ambiguas en: {acc}")
    a, t, g, c = (p3.count(base) for base in "ATGC")
    filas.append({
        "accession": acc,
        "species": especie[acc],
        "codons": len(p3),
        "GC1": 100*(p1.count("G")+p1.count("C"))/len(p1),
        "GC2": 100*(p2.count("G")+p2.count("C"))/len(p2),
        "GC3": 100*(g+c)/len(p3),
        "GC12": 50*((p1.count("G")+p1.count("C"))/len(p1)
                    +(p2.count("G")+p2.count("C"))/len(p2)),
        "PR2_A3_A3T3": a/(a+t),
        "PR2_G3_G3C3": g/(g+c)
    })
df = pd.DataFrame(filas)
assert len(df) == 30, f"Se esperaban 30 mitogenomas; se encontraron {len(df)}"
assert df["species"].nunique() == 14, "El número de especies no es 14"
df.to_csv("Codon_PR2_Neutrality_30mitogenomes.tsv", sep="\t", index=False)
medias = df.groupby("species", as_index=False).agg(
    n=("accession", "size"),
    PR2_AT=("PR2_A3_A3T3", "mean"),
    PR2_GC=("PR2_G3_G3C3", "mean")
)
medias.to_csv("Figure4_PR2_14_species.tsv", sep="\t", index=False)
fig, ax = plt.subplots(figsize=(8, 6))
for i, row in medias.iterrows():
    ax.scatter(row["PR2_GC"], row["PR2_AT"], s=70,
               color=plt.cm.tab20(i), edgecolor="black", linewidth=0.4,
               label=row["species"])
ax.axhline(0.5, color="gray", linestyle="--", linewidth=1)
ax.axvline(0.5, color="gray", linestyle="--", linewidth=1)
ax.set(xlim=(0, 1), ylim=(0, 1),
       xlabel="G3 / (G3 + C3)", ylabel="A3 / (A3 + T3)")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left",
          fontsize=8, frameon=False)
fig.tight_layout()
destino = Path("/mnt/e/Mitogenoma/nuevas figuras")
destino.mkdir(parents=True, exist_ok=True)
for extension in ("png", "pdf"):
    archivo = carpeta / f"Figure4_PR2_14_species.{extension}"
    fig.savefig(archivo, dpi=300, bbox_inches="tight")
    shutil.copy2(archivo, destino / archivo.name)
    print("Copiado:", destino / archivo.name)
plt.close(fig)
print("Mitogenomas:", len(df), "| Especies:", len(medias),
      "| Codones:", df["codons"].sum())
print("Rango PR2 A/T:", round(df["PR2_A3_A3T3"].min(), 3),
      "a", round(df["PR2_A3_A3T3"].max(), 3))
print("Rango PR2 G/C:", round(df["PR2_G3_G3C3"].min(), 3),
      "a", round(df["PR2_G3_G3C3"].max(), 3))
