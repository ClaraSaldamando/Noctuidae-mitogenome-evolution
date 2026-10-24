#!/usr/bin/env python3
"""Reconstructed, auditable CDS extraction and sense-codon filtering.

Input: 30 accession IDs (one per line) and corresponding GenBank files.
Output: CDS-level QC and filtered sequences. Requires Biopython.
This script is a documented reconstruction, not claimed as historical code.
"""
import argparse
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable

GENES = {'atp6','atp8','cox1','cox2','cox3','cytb','nad1','nad2',
         'nad3','nad4','nad4l','nad5','nad6'}
SENSE = set(CodonTable.unambiguous_dna_by_id[5].forward_table)

def gene_name(feature):
    name = feature.qualifiers.get('gene', [''])[0].lower().replace('-', '')
    if not name:
        product=feature.qualifiers.get('product',[''])[0].lower()
        if 'nadh' in product: name='nad'+product.split()[-1]
        elif 'oxidase' in product: name={'i':'cox1','ii':'cox2','iii':'cox3','1':'cox1','2':'cox2','3':'cox3'}.get(product.split()[-1],'')
        elif 'atp' in product: name='atp'+product.split()[-1]
        elif 'cytochrome b' in product: name='cytb'
    aliases = {'cob': 'cytb', 'cytb': 'cytb', 'coi': 'cox1', 'coii': 'cox2', 'coiii': 'cox3'}
    name = aliases.get(name, name)
    if name in ['co1','co2','co3']: name='cox'+name[-1]
    return 'nad'+name[2:] if name.startswith('nd') else name

def process_feature(feature, genome):
    # Biopython respects complement(), join(), and annotated strand.
    seq = str(feature.extract(genome.seq)).upper().replace('U', 'T')
    start = int(feature.qualifiers.get('codon_start', ['1'])[0])
    if start not in (1, 2, 3):
        raise ValueError(f'Invalid codon_start: {start}')
    seq = seq[start-1:]
    triplets = [seq[i:i+3] for i in range(0, len(seq)//3*3, 3)]
    valid = [c for c in triplets if c in SENSE]
    return start, valid, sum(c in CodonTable.unambiguous_dna_by_id[5].stop_codons for c in triplets), sum(
        len(c) == 3 and c not in SENSE and c not in CodonTable.unambiguous_dna_by_id[5].stop_codons for c in triplets), len(seq)%3

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genbank-dir', type=Path, required=True)
    p.add_argument('--accessions', type=Path, required=True,
                   help='TSV with accession column, or plain text one accession per line')
    p.add_argument('--out-dir', type=Path, default=Path('results/cds_filtering'))
    args = p.parse_args()
    lines = [s.strip() for s in args.accessions.read_text().splitlines() if s.strip()]
    accessions = [s.split('\t')[0] for s in lines]
    if accessions[0].lower() == 'accession':
        accessions.pop(0)
    if len(accessions) != 30 or len(set(accessions)) != 30:
        raise ValueError('Expected 30 unique accessions')
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows, sequences, n_codons = [], [], 0
    for accession in accessions:
        candidates = list(args.genbank_dir.glob(accession + '.*'))
        if len(candidates) != 1:
            raise FileNotFoundError(f'{accession}: expected one GenBank file; found {len(candidates)}')
        record = SeqIO.read(candidates[0], 'genbank')
        found = []
        for f in record.features:
            if f.type != 'CDS':
                continue
            gene = gene_name(f)
            if gene not in GENES:
                raise ValueError(f'Unexpected CDS gene {accession}: {gene}')
            start, valid, stops, ambiguous, remainder = process_feature(f, record)
            # Stop detection uses table 5, not the standard nuclear code.
            rows.append(dict(accession=accession, gene=gene, strand=f.location.strand,
                             codon_start=start, sense_codons=len(valid),
                             excluded_stop=stops, excluded_ambiguous=ambiguous,
                             excluded_terminal_bases=remainder))
            sequences.append((f'{accession}|{gene}', ''.join(valid)))
            n_codons += len(valid)
            found.append(gene)
        if len(found) != 13 or set(found) != GENES or len(set(found)) != 13:
            raise ValueError(f'{accession}: missing/duplicate genes: {found}')
    if len(rows) != 390 or n_codons != 111739:
        raise ValueError(f'QC mismatch: {len(rows)} CDS, {n_codons} sense codons')
    if any(r['codon_start'] != 1 for r in rows):
        raise ValueError('Expected codon_start=1 for all 390 CDS')
    with (args.out_dir/'cds_filtering_qc.tsv').open('w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)
    with (args.out_dir/'filtered_CDS_390.fasta').open('w') as h:
        for name, seq in sequences:
            h.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                h.write(seq[i:i+80]+'\n')
    print(f'Validated {len(rows)} CDS; {n_codons} sense codons; table 5')

if __name__ == '__main__':
    main()
