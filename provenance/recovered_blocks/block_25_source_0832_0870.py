from Bio import SeqIO
from pathlib import Path
r = list(SeqIO.parse("PCG13_31_aligned.fasta", "fasta"))
assert len(r) == 31 and {len(x.seq) for x in r} == {11229}
assert len({x.id for x in r}) == 31
assert any(x.id == "MN793322.1" for x in r)
lineas = [
    "#NEXUS",
    "begin data;",
    "  dimensions ntax=31 nchar=11229;",
    "  format datatype=dna missing=? gap=-;",
    "  matrix",
]
lineas += [f"  {x.id}  {str(x.seq).upper()}" for x in r]
lineas += [
    "  ;",
    "end;",
    "begin mrbayes;",
    "  set autoclose=yes nowarn=yes;",
    "  outgroup MN793322.1;",
    "  charset P1 = 1-678 679-843 5941-6951 6952-7305 10693-11229;",
    "  charset P2 = 844-2373 2374-3054 3055-3843 3844-4995;",
    "  charset P3 = 4996-5940 7306-8643 8644-8940 8941-10692;",
    "  partition genes = 3: P1, P2, P3;",
    "  set partition=genes;",
    "  lset applyto=(all) nst=6 rates=invgamma;",
    "  unlink revmat=(all) statefreq=(all) shape=(all) pinvar=(all);",
    "  prset applyto=(all) ratepr=variable;",
    "  mcmcp ngen=1000000 samplefreq=1000 printfreq=10000",
    "        diagnfreq=10000 nruns=2 nchains=4",
    "        filename=PCG13_31_BI;",
    "  mcmc;",
    "  sump;",
    "  sumt;",
    "end;",
]
archivo = Path("PCG13_31_BI.nex")
archivo.write_text("\n".join(lineas) + "\n")
print("Archivo MrBayes:", archivo, "| 31 secuencias | 11229 bp | 3 particiones")
