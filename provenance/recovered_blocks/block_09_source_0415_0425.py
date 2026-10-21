import pandas as pd
from scipy.stats import linregress
d = pd.read_csv("Codon_PR2_Neutrality_30mitogenomes.tsv", sep="\t")
m = d.groupby("species", as_index=False)[["GC3", "GC12"]].mean()
r = linregress(m["GC3"], m["GC12"])
m.to_csv("Figure5_neutrality_14_species.tsv", sep="\t", index=False)
print(f"Especies: {len(m)}")
print(f"Pendiente: {r.slope:.4f}")
print(f"Intercepto: {r.intercept:.4f}")
print(f"R²: {r.rvalue**2:.4f}")
print(f"p: {r.pvalue:.4g}")
