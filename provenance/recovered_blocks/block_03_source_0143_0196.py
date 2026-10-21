from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
folder = Path.cwd()
data = pd.read_csv(folder / "Supplementary_ENC_GC3s_30.tsv", sep="\t")
data["species"] = data["species"].where(
    ~data["species"].str.startswith("Helicoverpa armigera"),
    "Helicoverpa armigera"
)
means = data.groupby("species", as_index=False)[["ENC", "GC3s"]].mean()
if len(means) != 14:
    raise ValueError(f"Se esperaban 14 especies; hay {len(means)}")
means.to_csv(folder / "Figure3_species_means_14.tsv",
             sep="\t", index=False, float_format="%.6f")
colors = {
    "Agrotis": "#1f77b4",
    "Athetis": "#e69f00",
    "Helicoverpa": "#009e73",
    "Mamestra": "#d55e00",
    "Mythimna": "#56b4e9",
    "Sesamia": "#795548",
    "Spodoptera": "#c62828",
    "Trichoplusia": "#424242",
}
markers = ["o", "s", "^", "D", "P", "X", "v", "<",
           ">", "h", "8", "p", "*", "H"]
fig, ax = plt.subplots(figsize=(9, 6))
x = np.linspace(0.035, 0.09, 400)
expected = 2 + x + 29 / (x*x + (1-x)*(1-x))
ax.plot(x, expected, color="black", linestyle="--", linewidth=1.8,
        label="Expected ENC curve")
for i, row in means.iterrows():
    genus = row["species"].split()[0]
    ax.scatter(row["GC3s"], row["ENC"], s=85,
               color=colors[genus], marker=markers[i],
               edgecolor="black", linewidth=0.5,
               label=row["species"], zorder=3)
ax.set_xlim(0.035, 0.09)
ax.set_ylim(28, 38)
ax.set_xlabel("GC3s")
ax.set_ylabel("Effective number of codons (ENC)")
ax.grid(alpha=0.2)
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
          fontsize=8, frameon=False)
fig.tight_layout()
for extension in ("png", "pdf"):
    fig.savefig(folder / f"Figure3_ENC_GC3s_14_species_v2.{extension}",
                dpi=300, bbox_inches="tight")
plt.close(fig)
print("Especies:", len(means))
print("Figura:", folder / "Figure3_ENC_GC3s_14_species_v2.png")
