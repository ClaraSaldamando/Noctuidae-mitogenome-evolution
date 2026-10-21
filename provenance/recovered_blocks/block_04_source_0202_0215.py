import pandas as pd
from scipy.stats import linregress, pearsonr, spearmanr
for archivo, etiqueta in [
    ("Supplementary_ENC_GC3s_30.tsv", "30 mitogenomas"),
    ("Figure3_species_means_14.tsv", "14 medias por especie"),
]:
    df = pd.read_csv(archivo, sep="\t")
    print(f"\n{etiqueta}: {archivo}")
    print("Columnas:", ", ".join(df.columns))
    x, y = df["GC3s"], df["ENC"]
    reg = linregress(x, y)
    print(f"n={len(df)}; pendiente={reg.slope:.4f}; intercepto={reg.intercept:.4f}; R²={reg.rvalue**2:.4f}; p_regresión={reg.pvalue:.4g}")
    print(f"Pearson: r={pearsonr(x, y).statistic:.4f}; p={pearsonr(x, y).pvalue:.4g}")
    print(f"Spearman: rho={spearmanr(x, y).statistic:.4f}; p={spearmanr(x, y).pvalue:.4g}")
