from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
base = Path.cwd()
source = base / "01_raw_data/final_genbank"
out = base / "revision_30_noctuidos"
duplicates = {
    "NC_022689.1", "NC_022185.1", "NC_023118.1",
    "NC_053623", "NC_061705.1", "NC_054179",
    "NC_061562", "NC_015835.1", "NC_019622"
}
sense = set(CodonTable.unambiguous_dna_by_id[5].forward_table)
records = []
for path in sorted(source.glob("*.gb")):
    accession = path.stem
    if accession in duplicates or accession == "MN793322.1":
        continue
    record = SeqIO.read(path, "genbank")
    codons = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        sequence = str(feature.extract(record.seq)).upper()
        codons.extend(sequence[i:i+3]
                      for i in range(0, len(sequence) - 2, 3)
                      if sequence[i:i+3] in sense)
    records.append((accession, "".join(codons)))
assert len(records) == 30
assert sum(len(sequence) // 3 for _, sequence in records) == 111739
with (out / "codonw_input_30.fasta").open("w") as file:
    for accession, sequence in records:
        file.write(f">{accession}\n")
        for i in range(0, len(sequence), 80):
            file.write(sequence[i:i+80] + "\n")
print("Entrada CodonW:", len(records), "mitogenomas")
