from Bio import SeqIO, Phylo
a = "04_alignments/gene_alignments_final40_clean/cox1_clean_NT.fasta"
t = "04_alignments/concatenated_final40/PCG13_final40_by_gene.treefile"
print("Ejemplos del alineamiento:")
print(*[r.id for r in list(SeqIO.parse(a, "fasta"))[:20]], sep="\n")
print("\nEjemplos del árbol:")
print(*[x.name for x in Phylo.read(t, "newick").get_terminals()[:20]], sep="\n")
