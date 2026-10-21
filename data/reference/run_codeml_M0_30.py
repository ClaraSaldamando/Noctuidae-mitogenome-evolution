#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO, Phylo
import re
import shutil
import subprocess
import sys
import csv

# ============================================================
# CODEML M0 – 13 genes mitocondriales, 30 mitogenomas
# Código genético mitocondrial de invertebrados: icode = 4
# ============================================================

project = Path.cwd()

alignment_dir = (
    project
    / "revision_30_noctuidos"
    / "dnds_input_30"
)

tree_source = (
    project
    / "revision_30_noctuidos"
    / "dnds_input_30"
    / "PCG13_pruned_30.nwk"
)

output_dir = project / "revision_30_noctuidos" / "dnds_M0_30"

genes = [
    "atp6",
    "atp8",
    "cox1",
    "cox2",
    "cox3",
    "cytb",
    "nad1",
    "nad2",
    "nad3",
    "nad4",
    "nad4l",
    "nad5",
    "nad6",
]

# ------------------------------------------------------------
# Comprobar codeml
# ------------------------------------------------------------

codeml_path = shutil.which("codeml")

if codeml_path is None:
    sys.exit("ERROR: codeml no está disponible en PATH.")

print("codeml:", codeml_path)

# ------------------------------------------------------------
# Leer y limpiar árbol
# ------------------------------------------------------------

if not tree_source.exists():
    sys.exit(f"ERROR: no se encontró el árbol: {tree_source}")

tree_text = tree_source.read_text().strip()

# Eliminar etiquetas internas de soporte:
# Ejemplo: )99.5/100:0.01  ->  ):0.01
tree_clean = re.sub(
    r"\)([^():,;]+):",
    "):",
    tree_text
)

# Asegurar punto y coma final
if not tree_clean.endswith(";"):
    tree_clean += ";"

# ------------------------------------------------------------
# Obtener identificadores del árbol
# ------------------------------------------------------------

temporary_tree = project / "tree_temp_clean.nwk"
temporary_tree.write_text(tree_clean + "\n")

tree = Phylo.read(temporary_tree, "newick")
tree_ids = {tip.name for tip in tree.get_terminals()}

temporary_tree.unlink()

print("Número de puntas del árbol:", len(tree_ids))

if len(tree_ids) != 30:
    sys.exit(
        f"ERROR: el árbol tiene {len(tree_ids)} puntas, no 30."
    )

# ------------------------------------------------------------
# Crear directorio principal
# ------------------------------------------------------------

output_dir.mkdir(parents=True, exist_ok=True)

results = []

# ------------------------------------------------------------
# Procesar cada gen
# ------------------------------------------------------------

for gene in genes:

    fasta_file = alignment_dir / f"{gene}_30_NT.fasta"

    if not fasta_file.exists():
        sys.exit(f"ERROR: no se encontró {fasta_file}")

    records = list(SeqIO.parse(fasta_file, "fasta"))

    if len(records) != 30:
        sys.exit(
            f"ERROR: {gene} contiene {len(records)} secuencias."
        )

    alignment_ids = {record.id for record in records}

    missing_in_tree = sorted(alignment_ids - tree_ids)
    missing_in_alignment = sorted(tree_ids - alignment_ids)

    if missing_in_tree or missing_in_alignment:
        print(f"\nERROR de identificadores en {gene}")

        if missing_in_tree:
            print(
                "Presentes en alineamiento pero no en árbol:",
                ", ".join(missing_in_tree)
            )

        if missing_in_alignment:
            print(
                "Presentes en árbol pero no en alineamiento:",
                ", ".join(missing_in_alignment)
            )

        sys.exit(1)

    lengths = {len(record.seq) for record in records}

    if len(lengths) != 1:
        sys.exit(
            f"ERROR: {gene} tiene longitudes diferentes: {lengths}"
        )

    seq_length = next(iter(lengths))

    if seq_length % 3 != 0:
        sys.exit(
            f"ERROR: la longitud de {gene} no es múltiplo de 3."
        )

    gene_dir = output_dir / gene
    gene_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Crear archivo de secuencias en formato PAML
    # --------------------------------------------------------

    paml_seq = gene_dir / f"{gene}.phy"

    with paml_seq.open("w") as handle:
        handle.write(f"{len(records)} {seq_length}\n")

        for record in records:
            sequence = str(record.seq).upper()

            # MACSE puede usar ! para frameshifts.
            # En alineamientos limpios no debería haberlos.
            if "!" in sequence:
                sys.exit(
                    f"ERROR: se encontró ! en {gene}, {record.id}"
                )

            handle.write(f"{record.id}\n")
            handle.write(f"{sequence}\n")

    # --------------------------------------------------------
    # Crear árbol PAML
    # --------------------------------------------------------

    tree_file = gene_dir / "tree_final30.nwk"

    tree_file.write_text(
        f"30 1\n{tree_clean}\n"
    )

    # --------------------------------------------------------
    # Crear archivo codeml.ctl
    # --------------------------------------------------------

    ctl_file = gene_dir / "codeml.ctl"

    ctl_text = f"""
      seqfile = {gene}.phy
     treefile = tree_final30.nwk
      outfile = {gene}_M0.mlc

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0

        model = 0
      NSsites = 0

        icode = 4
    fix_kappa = 0
        kappa = 2

    fix_omega = 0
        omega = 0.1

    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 8

        getSE = 0
 RateAncestor = 0

   Small_Diff = 0.5e-6

    cleandata = 1
       method = 0
"""

    ctl_file.write_text(ctl_text.strip() + "\n")

    # --------------------------------------------------------
    # Ejecutar codeml
    # --------------------------------------------------------

    print(
        f"\nEjecutando {gene}: "
        f"n={len(records)}, longitud={seq_length}"
    )

    log_file = gene_dir / f"{gene}_codeml.log"

    with log_file.open("w") as log:
        process = subprocess.run(
            [codeml_path, "codeml.ctl"],
            cwd=gene_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )

    if process.returncode != 0:
        print(
            f"ERROR: codeml falló para {gene}. "
            f"Revisa {log_file}"
        )

        results.append({
            "Gene": gene,
            "n_sequences": len(records),
            "alignment_length_nt": seq_length,
            "omega": None,
            "lnL": None,
            "status": "ERROR"
        })

        continue

    mlc_file = gene_dir / f"{gene}_M0.mlc"

    if not mlc_file.exists():
        print(f"ERROR: no se creó {mlc_file}")

        results.append({
            "Gene": gene,
            "n_sequences": len(records),
            "alignment_length_nt": seq_length,
            "omega": None,
            "lnL": None,
            "status": "NO_MLC"
        })

        continue

    mlc_text = mlc_file.read_text(
        errors="replace"
    )

    # --------------------------------------------------------
    # Extraer omega
    # --------------------------------------------------------

    omega_matches = re.findall(
        r"omega\s+\(dN/dS\)\s*=\s*([0-9.eE+-]+)",
        mlc_text
    )

    if not omega_matches:
        omega_matches = re.findall(
            r"\bomega\s*=\s*([0-9.eE+-]+)",
            mlc_text
        )

    omega = (
        float(omega_matches[-1])
        if omega_matches
        else None
    )

    # --------------------------------------------------------
    # Extraer log-likelihood
    # --------------------------------------------------------

    lnl_matches = re.findall(
        r"lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*"
        r"([\-0-9.]+)",
        mlc_text
    )

    lnL = (
        float(lnl_matches[-1])
        if lnl_matches
        else None
    )

    status = "OK" if omega is not None else "OMEGA_NOT_FOUND"

    results.append({
        "Gene": gene,
        "n_sequences": len(records),
        "alignment_length_nt": seq_length,
        "omega": omega,
        "lnL": lnL,
        "status": status
    })

    print(
        f"Resultado {gene}: "
        f"omega={omega}, lnL={lnL}, estado={status}"
    )

# ------------------------------------------------------------
# Guardar tabla final
# ------------------------------------------------------------

results_file = (
    output_dir
    / "dNdS_M0_13genes_final30.csv"
)

with results_file.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "Gene",
            "n_sequences",
            "alignment_length_nt",
            "omega",
            "lnL",
            "status"
        ]
    )

    writer.writeheader()
    writer.writerows(results)

print("\n========================================")
print("Análisis terminado")
print("========================================")
print("Resultados:", results_file)

print("\nResumen:")

for row in results:
    print(
        f"{row['Gene']:6s} "
        f"omega={str(row['omega']):>12s} "
        f"estado={row['status']}"
    )
