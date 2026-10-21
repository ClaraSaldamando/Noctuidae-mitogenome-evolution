from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
carpeta = Path.cwd()
d = pd.read_csv("Figure5_neutrality_14_species.tsv", sep="\t")
assert len(d) == 14
r = linregress(d["GC3"], d["GC12"])
fig, ax = plt.subplots(figsize=(8, 6))
for i, row in d.iterrows():
    ax.scatter(row["GC3"], row["GC12"], s=75,
               color=plt.cm.tab20(i), edgecolor="black",
               linewidth=0.4, label=row["species"], zorder=3)
x = np.linspace(d["GC3"].min(), d["GC3"].max(), 100)
ax.plot(x, r.intercept + r.slope*x, color="#15616D", linewidth=2)
ax.set_xlabel("GC3 (%)")
ax.set_ylabel("GC12 (%)")
ax.text(0.03, 0.97,
        f"GC12 = {r.intercept:.3f} + {r.slope:.3f} × GC3\n"
        f"R² = {r.rvalue**2:.3f}; p = {r.pvalue:.4f}; n = {len(d)}",
        transform=ax.transAxes, va="top", fontsize=10)
ax.grid(alpha=0.2)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left",
          fontsize=8, frameon=False)
fig.tight_layout()
destino = Path("/mnt/e/Mitogenoma/nuevas figuras")
destino.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    archivo = carpeta / f"Figure5_neutrality_14_species.{ext}"
    fig.savefig(archivo, dpi=300, bbox_inches="tight")
    shutil.copy2(archivo, destino / archivo.name)
    print("Copiado:", destino / archivo.name)
plt.close(fig)
