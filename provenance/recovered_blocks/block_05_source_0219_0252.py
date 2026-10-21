from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
folder = Path.cwd()
df = pd.read_csv(folder / "Figure3_species_means_14.tsv", sep="\t")
reg = linregress(df["GC3s"], df["ENC"])
x = np.linspace(0.035, 0.090, 300)
enc_esperado = 2 + x + 29 / (x**2 + (1-x)**2)
fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(x, enc_esperado, "--", color="#777777",
        linewidth=2, label="Expected ENC (GC3s only)")
ax.plot(x, reg.intercept + reg.slope*x, color="#15616D",
        linewidth=2, label=f"Observed regression (R² = {reg.rvalue**2:.3f})")
colores = plt.cm.tab20(np.linspace(0, 1, len(df)))
for (_, row), color in zip(df.iterrows(), colores):
    ax.scatter(row["GC3s"], row["ENC"], s=65, color=color,
               edgecolor="black", linewidth=0.4, label=row["species"], zorder=3)
ax.set(xlim=(0.035, 0.090), ylim=(28, 38),
       xlabel="GC3s", ylabel="Effective number of codons (ENC)")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8,
          frameon=False)
ax.grid(alpha=0.2)
fig.tight_layout()
destino = Path("/mnt/e/Mitogenoma/nuevas figuras")
destino.mkdir(parents=True, exist_ok=True)
for extension in ("png", "pdf"):
    archivo = folder / f"Figure3_ENC_GC3s_regression_14_species.{extension}"
    fig.savefig(archivo, dpi=300, bbox_inches="tight")
    shutil.copy2(archivo, destino / archivo.name)
    print(destino / archivo.name)
plt.close(fig)
