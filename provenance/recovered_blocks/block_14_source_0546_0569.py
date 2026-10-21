from pathlib import Path
origen = Path("run_codeml_M0_final40.py")
destino = Path("revision_30_noctuidos/run_codeml_M0_30.py")
s = origen.read_text()
cambios = [
    ('project\n    / "04_alignments"\n    / "gene_alignments_final40_clean"',
     'project\n    / "revision_30_noctuidos"\n    / "dnds_input_30"'),
    ('project\n    / "04_alignments"\n    / "concatenated_final40"\n    / "PCG13_final40_by_gene.treefile"',
     'project\n    / "revision_30_noctuidos"\n    / "dnds_input_30"\n    / "PCG13_pruned_30.nwk"'),
    ('project / "07_dnds_M0_final40"',
     'project / "revision_30_noctuidos" / "dnds_M0_30"'),
    ('f"{gene}_clean_NT.fasta"', 'f"{gene}_30_NT.fasta"'),
]
for viejo, nuevo in cambios:
    if s.count(viejo) != 1:
        raise SystemExit(f"No se encontró exactamente una vez: {viejo}")
    s = s.replace(viejo, nuevo)
s = s.replace("final40", "final30")
s = s.replace("40 mitogenomas", "30 mitogenomas")
s = s.replace("!= 40", "!= 30")
s = s.replace("no 40", "no 30")
s = s.replace('f"40 1\\n', 'f"30 1\\n')
destino.write_text(s)
print("Programa creado:", destino)
