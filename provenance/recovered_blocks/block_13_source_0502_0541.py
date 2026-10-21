from pathlib import Path
import pandas as pd
from Bio import SeqIO, Phylo
base = Path.cwd()
salida = base / "revision_30_noctuidos/dnds_input_30"
origen = base / "04_alignments/gene_alignments_final40_clean"
arbol_origen = base / "04_alignments/concatenated_final40/PCG13_final40_by_gene.treefile"
def sin_version(nombre):
    return nombre.rsplit(".", 1)[0] if nombre.rsplit(".", 1)[-1].isdigit() else nombre
tabla = pd.read_csv(base / "revision_30_noctuidos/Supplementary_ENC_GC3s_30.tsv",
                    sep="\t")
elegidas = {sin_version(str(x)) for x in tabla["accession"]}
assert len(elegidas) == 30
arbol = Phylo.read(arbol_origen, "newick")
puntas = list(arbol.get_terminals())
assert len({sin_version(x.name) for x in puntas}) == len(puntas)
assert elegidas <= {sin_version(x.name) for x in puntas}
genes = "atp6 atp8 cox1 cox2 cox3 cytb nad1 nad2 nad3 nad4 nad4l nad5 nad6".split()
seleccionados = {}
for gen in genes:
    registros = list(SeqIO.parse(origen / f"{gen}_clean_NT.fasta", "fasta"))
    ids = {sin_version(r.id) for r in registros}
    assert len(ids) == len(registros), f"Identificadores repetidos: {gen}"
    assert elegidas <= ids, f"Faltan accesiones en {gen}: {sorted(elegidas-ids)}"
    elegidos = [r for r in registros if sin_version(r.id) in elegidas]
    assert len(elegidos) == 30
    assert len({len(r.seq) for r in elegidos}) == 1
    assert len(elegidos[0].seq) % 3 == 0
    seleccionados[gen] = elegidos
salida.mkdir(parents=True, exist_ok=True)
for gen, registros in seleccionados.items():
    SeqIO.write(registros, salida / f"{gen}_30_NT.fasta", "fasta")
for punta in list(arbol.get_terminals()):
    if sin_version(punta.name) not in elegidas:
        arbol.prune(punta)
assert len(arbol.get_terminals()) == 30
Phylo.write(arbol, salida / "PCG13_pruned_30.nwk", "newick")
print("Preparados:", len(genes), "alineamientos de 30 secuencias")
print("Árbol:", len(arbol.get_terminals()), "puntas")
print("Carpeta:", salida)
