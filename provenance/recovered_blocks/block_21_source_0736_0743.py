from Bio import SeqIO, Phylo
aln = list(SeqIO.parse("PCG13_31_aligned.fasta", "fasta"))
tips = [x.name for x in Phylo.read("PCG13_31_ML.treefile", "newick").get_terminals()]
print("Alineamiento:", len(aln), "secuencias;", len(aln[0].seq), "bp")
print("Árbol ML:", len(tips), "puntas")
print("Identificadores coinciden:", {x.id for x in aln} == set(tips))
print("\nEsquema seleccionado:")
print(open("PCG13_31_ML.best_scheme.nex").read()[:2500])
