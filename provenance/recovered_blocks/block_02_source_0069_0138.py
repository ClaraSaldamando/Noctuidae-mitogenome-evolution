from pathlib import Path
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
folder = Path.cwd()
enc = pd.read_csv(folder / "enc_gc3s_30.tsv", sep=r"\s+")
qc = pd.read_csv(folder / "codon_qc_30.tsv", sep="\t")
enc = enc.rename(columns={"title": "accession", "Nc": "ENC"})
data = qc[["accession", "species"]].merge(
    enc[["accession", "ENC", "GC3s"]],
    on="accession", validate="one_to_one"
)
if len(data) != 30 or data[["ENC", "GC3s"]].isna().any().any():
    raise ValueError("La tabla no contiene 30 valores completos")
if not data["GC3s"].between(0, 1).all():
    raise ValueError("GC3s está fuera del intervalo 0–1")
data["genus"] = data["species"].str.split().str[0]
data = data.sort_values(["species", "accession"])
data[["accession", "species", "ENC", "GC3s"]].to_csv(
    folder / "Supplementary_ENC_GC3s_30.tsv",
    sep="\t", index=False, float_format="%.6f"
)
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
x = np.linspace(0, 1, 500)
expected = 2 + x + 29 / (x*x + (1-x)*(1-x))
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
for ax in axes:
    ax.plot(x, expected, color="black", linewidth=1.5,
            label="Expected ENC")
    for genus, group in data.groupby("genus"):
        ax.scatter(group["GC3s"], group["ENC"], s=55,
                   color=colors[genus], edgecolor="white",
                   linewidth=0.5, label=genus, zorder=3)
    ax.set_xlabel("GC3s")
    ax.set_ylabel("Effective number of codons (ENC)")
    ax.grid(alpha=0.2)
axes[0].set_title("Full expected curve")
axes[0].set_xlim(0, 1)
axes[0].set_ylim(20, 62)
pad_x = max(0.015, (data.GC3s.max() - data.GC3s.min()) * 0.2)
axes[1].set_title("Observed mitogenomes")
axes[1].set_xlim(max(0, data.GC3s.min()-pad_x),
                 data.GC3s.max()+pad_x)
axes[1].set_ylim(data.ENC.min()-2, data.ENC.max()+2)
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=5,
           frameon=False, bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(rect=(0, 0.09, 1, 1))
for extension in ("png", "pdf"):
    fig.savefig(folder / f"Figure3_ENC_GC3s_30.{extension}",
                dpi=300, bbox_inches="tight")
plt.close(fig)
print("Mitogenomas:", len(data))
print("ENC:", round(data.ENC.min(), 2), "a",
      round(data.ENC.max(), 2))
print("GC3s:", round(data.GC3s.min(), 3), "a",
      round(data.GC3s.max(), 3))
print("Figura:", folder / "Figure3_ENC_GC3s_30.png")
