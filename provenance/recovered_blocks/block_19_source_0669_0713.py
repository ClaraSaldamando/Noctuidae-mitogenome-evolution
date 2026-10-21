from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
base = Path.cwd()
origen = base / "04_alignments/gene_alignments_final40_clean"
salida = base / "revision_30_noctuidos/phylogeny_31"
salida.mkdir(parents=True, exist_ok=True)
tabla = pd.read_csv(base / "revision_30_noctuidos/Supplementary_ENC_GC3s_30.tsv",
                    sep="\t")
def canon(x):
    return x.rsplit(".", 1)[0] if x.rsplit(".", 1)[-1].isdigit() else x
elegidas = {canon(str(x)) for x in tabla["accession"]}
elegidas.add(canon("MN793322.1"))
assert len(elegidas) == 31
genes = "atp6 atp8 cox1 cox2 cox3 cytb nad1 nad2 nad3 nad4 nad4l nad5 nad6".split()
partes = {}
nombres = None
for gen in genes:
    registros = list(SeqIO.parse(origen / f"{gen}_clean_NT.fasta", "fasta"))
    seleccion = {r.id: str(r.seq).upper() for r in registros
                 if canon(r.id) in elegidas}
    assert len(seleccion) == 31, f"{gen}: {len(seleccion)} secuencias"
    assert len({len(s) for s in seleccion.values()}) == 1
    assert len(next(iter(seleccion.values()))) % 3 == 0
    assert all("!" not in s for s in seleccion.values()), f"Frameshift: {gen}"
    if nombres is None:
        nombres = set(seleccion)
    assert set(seleccion) == nombres, f"Nombres distintos en {gen}"
    partes[gen] = seleccion
orden = sorted(nombres)
secuencias = [
    SeqRecord(Seq("".join(partes[gen][acc] for gen in genes)),
              id=acc, description="")
    for acc in orden
]
archivo = salida / "PCG13_31_aligned.fasta"
SeqIO.write(secuencias, archivo, "fasta")
assert len(secuencias) == 31
assert {len(r.seq) for r in secuencias} == {11229}
print("Alineamiento:", archivo)
print("Secuencias:", len(secuencias))
print("Longitud:", len(secuencias[0].seq), "bp")
print("Grupo externo:", "MN793322.1" in orden)
