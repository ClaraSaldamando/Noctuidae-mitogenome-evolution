from Bio import Phylo
t = Phylo.read("PCG13_31_ML.treefile", "newick")
internos = t.get_nonterminals()
print("Nodos con más de dos descendientes:",
      [(len(c.clades), c.name) for c in internos if len(c.clades) > 2])
print("Ramas internas de longitud cero o casi cero:",
      sum((c.branch_length or 0) < 1e-7 for c in internos))
