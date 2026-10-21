from pathlib import Path
import shutil
import pandas as pd
import matplotlib.pyplot as plt
base = Path.cwd()
d = pd.read_csv(base / "dnds_M0_30/dNdS_M0_13genes_final30.csv")
assert len(d) == 13 and d["status"].eq("OK").all()
grupos = {
    "atp6": "ATP synthase", "atp8": "ATP synthase",
    "cox1": "Cytochrome c oxidase", "cox2": "Cytochrome c oxidase",
    "cox3": "Cytochrome c oxidase", "cytb": "Cytochrome b",
}
d["group"] = d["Gene"].map(lambda g: grupos.get(g, "NADH dehydrogenase"))
colores = {
    "ATP synthase": "#D08B3E",
    "Cytochrome c oxidase": "#247BA0",
    "Cytochrome b": "#54A24B",
    "NADH dehydrogenase": "#B65F73",
}
fig, ax = plt.subplots(figsize=(9, 5))
barras = ax.bar(d["Gene"], d["omega"],
                color=d["group"].map(colores), edgecolor="black",
                linewidth=0.4)
ax.bar_label(barras, fmt="%.4f", padding=3, fontsize=8, rotation=90)
ax.set_ylabel("dN/dS (ω), CODEML M0")
ax.set_xlabel("Mitochondrial protein-coding gene")
ax.set_ylim(0, d["omega"].max()*1.35)
ax.spines[["top", "right"]].set_visible(False)
for nombre, color in colores.items():
    ax.bar([], [], color=color, label=nombre)
ax.legend(frameon=False, ncol=2, fontsize=8, loc="upper left")
fig.tight_layout()
destino = Path("/mnt/e/Mitogenoma/nuevas figuras")
destino.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    archivo = base / f"Figure6_dNdS_M0_30.{ext}"
    fig.savefig(archivo, dpi=300, bbox_inches="tight")
    shutil.copy2(archivo, destino / archivo.name)
    print("Copiado:", destino / archivo.name)
tabla = destino / "dNdS_M0_13genes_30.csv"
shutil.copy2(base / "dnds_M0_30/dNdS_M0_13genes_final30.csv", tabla)
print("Copiado:", tabla)
plt.close(fig)
