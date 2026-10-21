# Noctuidae mitochondrial genome evolution

Reproducibility materials for **Mitochondrial Genome Evolution in Noctuidae: Codon Usage Bias, Nucleotide Composition, and Phylogenetic Relationships Among Major Agricultural Pests**, Clara I. Saldamando-Benjumea, submitted to Mitochondrial DNA Part A, GDNA-FLRP-2026-0036.R2.

## Scope and provenance
This package combines archived inputs/outputs and original terminal code with **newly reconstructed, tested Python scripts**. Reconstructed scripts are explicitly identified in their docstrings; they must not be described as the exact original scripts. All 30 noctuid accessions are specified in data/accessions_30.tsv. The phylogenetic outgroup is Ostrinia nubilalis MN793322.1. The original 40 GenBank records are preserved so that the nine duplicate exclusions can be checked.

## Quick start
Use Python 3.12 in a new environment, then run from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/01_reproduce_tables.py
python scripts/02_prepare_phylogeny.py
python scripts/03_plot_results.py
```

Outputs are written under results/. Original inputs and reference outputs under data/ are not edited. The plotter makes new renderings of the analyses, not byte-identical copies of historical figures. Archived original figures are retained in data/reference/.

## What was verified
- 390 annotated CDS from 30 mitogenomes and 111,739 valid sense codons, using NCBI table 5.
- All retained CDS annotations have codon_start=1; applying codon_start explicitly preserves the original sequences.
- Reconstructed RSCU agrees with all 1,860 archived accession/codon values within 0.000001; all 868 species/codon means agree within the same tolerance.
- Reconstructed CodonW input sequences match the archived input exactly after normalizing accession version suffixes.
- Nine excluded records have full nucleotide sequences identical to retained records.
- Neutrality regression from reconstructed sequences: GC12 = 27.21927049 + 0.11379699*GC3; R squared 0.32014902; p 0.03494254.
- Gene category summaries and correlations reproduce the rounded supplementary values. RSCU is undefined for absent amino acid families; these undefined values are omitted from the gene-level mean of A/T-ending codons, rather than assigned zero. This matters for atp8.
- Final concatenation reconstructed from the archived cleaned gene alignments matches the archived 31-sequence, 11,229-position alignment exactly.

Detailed results are in results/validation.json, neutrality_diagnostics.json, neutrality_leave_one_out.tsv, Supplementary_S1.tsv to Supplementary_S4.tsv and phylogeny/alignment_validation.json.

## Interpretation of metrics
Species RSCU is the unweighted mean of accession-specific RSCU, not RSCU calculated from pooled species codons. Gene-level RSCU pools counts over 30 accessions. Gene-level GC3s is the unweighted mean over accessions at positions permitting a synonymous third-base change under table 5. These gene-level values are distinct from the imported CodonW species-level ENC and GC3s outputs. PR2 uses all retained third positions. Ambiguous and stop triplets and incomplete terminal triplets are excluded before pooled codon counting. Codon_start is applied before grouping triplets.

## External programs
The original log records IQ-TREE 2.0.7, seed 20260927, ModelFinder MFP+MERGE, 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates. Saved MrBayes configurations/logs are provided. CODEML controls use model=0, NSsites=0, CodonFreq=2 and icode=4. Optional reruns require programs already installed on PATH:

```bash
python scripts/04_run_external.py codonw
python scripts/04_run_external.py ml
python scripts/04_run_external.py bi
python scripts/04_run_external.py m0
```

Run the quick-start preparation first. New external outputs use separate names. MrBayes receives a fresh two-million-generation configuration derived from the archived initial file; the original continuation files are preserved separately. These external programs were **not rerun** during package preparation. ENC values and omega used in the Python pipeline are imported archived estimates. Running the external programs does not automatically replace them; compare new estimates with archived outputs explicitly.

## Remaining upstream documentation
The pipeline starts phylogenetic reconstruction from archived, cleaned MACSE gene alignments. The complete original MACSE command and cleaning decisions were not recovered. Do not claim an end-to-end reproduction from unaligned GenBank CDS to final phylogeny until this upstream step is documented and validated. Initial download command provenance and full original software-version evidence also remain incomplete. The included GenBank snapshots make the downstream analyses independent of online changes.

The recovered terminal plotting code uses Python/matplotlib. The manuscript's descriptions of R plotting should be reconciled with the actual production files before editor notification. Genetic code numbering is program-specific: NCBI table 5 and PAML icode=4 are not interchangeable labels; verify CodonW's own -code 4 setting in its original program version.

## Layout
- scripts/: portable, reconstructed and tested analysis/preparation/plotting code, plus optional external runners.
- data/genbank/: original annotated sequence snapshots.
- data/accessions_30.tsv: final noctuid accession manifest.
- data/clean_gene_alignments/: archived cleaned gene alignments.
- data/reference/: archived final tables, configurations, alignments, trees, logs and figures.
- results/: newly generated tables, validation records and figures.
- provenance/recovered_blocks/: complete original Python blocks recovered from the supplied terminal export; source line ranges and working directories are indexed. These include diagnostics and intermediate attempts and are not a batch workflow.
- Supplementary_analytical_scripts.pdf and .md: documentation and the same portable source code.

## Reuse and citation
Author: Clara I. Saldamando-Benjumea. Cite the manuscript title and its accession list when reusing study materials. Public sequence records retain their original provenance. A software license and permanent release identifier should be selected by the author before a FAIR release; this review package does not invent either. No DOI or repository URL has yet been assigned.
