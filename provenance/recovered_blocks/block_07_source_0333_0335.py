import pandas as pd
d = pd.read_csv("Supplementary_ENC_GC3s_30.tsv", sep="\t")
print(d.groupby("species").size().to_string())
