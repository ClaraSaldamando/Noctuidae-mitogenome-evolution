#!/usr/bin/env python3
"""Reconstructed reproducibility implementation, not recovered original code.

Run from any directory. Read preserved GenBank records and accession manifest;
write new results without modifying historical reference files. NCBI table 5.
RSCU species means average accession-specific RSCU; gene summaries pool codons.
CODEML omega values are imported from the archived M0 summary, not newly estimated.
"""
from pathlib import Path
from collections import Counter, defaultdict
import argparse, json, math, platform
import numpy as np
import pandas as pd
from scipy import stats
from Bio import SeqIO, __version__ as biopython_version
from Bio.Data import CodonTable

GENES = 'atp6 atp8 cox1 cox2 cox3 cytb nad1 nad2 nad3 nad4 nad4l nad5 nad6'.split()
CODE = CodonTable.unambiguous_dna_by_id[5].forward_table
FAMILIES = defaultdict(list)
for codon, aa in CODE.items(): FAMILIES[aa].append(codon)

def canon(acc):
    return acc.rsplit('.', 1)[0] if acc.rsplit('.', 1)[-1].isdigit() else acc

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

def rscu(counts):
    values = {}
    for aa, codons in FAMILIES.items():
        n = sum(counts[c] for c in codons)
        values.update({c: counts[c] * len(codons) / n if n else float('nan') for c in codons})
    return values

def composition(sequence):
    p = [sequence[i::3] for i in range(3)]
    gc = [100 * (v.count('G') + v.count('C')) / len(v) for v in p]
    a,t,g,c = (p[2].count(b) for b in 'ATGC')
    return dict(GC1=gc[0], GC2=gc[1], GC3=gc[2], GC12=(gc[0]+gc[1])/2,
                PR2_A3_A3T3=a/(a+t), PR2_G3_G3C3=g/(g+c))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=Path, default=None)
    args=parser.parse_args(); root=args.root.resolve(); out=args.out or root/'results'
    out.mkdir(parents=True, exist_ok=True)
    ref=root/'data/reference'; raw=root/'data/genbank'
    manifest=pd.read_csv(root/'data/accessions_30.tsv', sep='\t')
    assert len(manifest)==30 and manifest.accession.nunique()==30
    source={canon(f.stem):f for f in raw.glob('*.gb')}
    records=[]; long=[]; qc=[]; comps=[]; genome=[]; gene_rows=[]; gene_counts=defaultdict(Counter)
    total=0
    for row in manifest.itertuples():
        acc=str(row.accession); species=str(row.species)
        species='Helicoverpa armigera' if species.startswith('Helicoverpa armigera') else species
        rec=SeqIO.read(source[canon(acc)], 'genbank'); pooled=Counter(); seqparts=[]; genes_seen=[]
        for feature in rec.features:
            if feature.type!='CDS': continue
            gene=gene_name(feature)
            if gene not in GENES: raise ValueError(f'Unrecognized CDS: {acc} {gene}')
            offset=int(feature.qualifiers.get('codon_start',['1'])[0])-1
            seq=str(feature.extract(rec.seq)).upper()[offset:]
            valid=[seq[i:i+3] for i in range(0,len(seq)-2,3) if seq[i:i+3] in CODE]
            counts=Counter(valid); pooled.update(counts); seqparts.extend(valid); genes_seen.append(gene)
            # Third-position synonyms are evaluated within the first two bases.
            syn=[c for c in valid if any(CODE.get(c[:2]+b)==CODE[c] for b in 'ACGT' if b!=c[2])]
            gc3s=100*sum(c[2] in 'GC' for c in syn)/len(syn)
            gene_rows.append(dict(accession=acc,gene=gene,codons=len(valid),codon_start=offset+1,GC3s=gc3s))
            gene_counts[gene].update(counts)
        assert sorted(genes_seen)==sorted(GENES), (acc,genes_seen)
        sequence=''.join(seqparts); records.append((acc,sequence)); total+=len(seqparts)
        qc.append(dict(accession=acc,species=species,CDS=len(genes_seen),sense_codons=len(seqparts)))
        comps.append(dict(accession=acc,species=species,**composition(sequence)))
        for c,value in sorted(rscu(pooled).items()):
            long.append(dict(accession=acc,species=species,amino_acid=CODE[c],codon=c,count=pooled[c],RSCU=value))
        count=Counter(str(rec.seq).upper()); n=sum(count[b] for b in 'ACGT')
        genome.append(dict(accession=acc,species=species,length=len(rec),unambiguous_bases=n,
                           **{b:100*count[b]/n for b in 'ACGT'},
                           AT_skew=(count['A']-count['T'])/(count['A']+count['T']),
                           GC_skew=(count['G']-count['C'])/(count['G']+count['C'])))
    assert total==111739 and len(gene_rows)==390
    df=pd.DataFrame(long); df.to_csv(out/'rscu_30.tsv',sep='\t',index=False)
    means=df.groupby(['species','amino_acid','codon'],as_index=False).agg(genomes=('accession','size'),mean_RSCU=('RSCU','mean'))
    means.to_csv(out/'rscu_medias_14_especies.tsv',sep='\t',index=False)
    assert means.species.nunique()==14 and len(means)==14*62
    pd.DataFrame(genome).to_csv(out/'genome_composition.tsv',sep='\t',index=False)
    pd.DataFrame(qc).to_csv(out/'codon_qc_30.tsv',sep='\t',index=False)
    comp=pd.DataFrame(comps); comp.to_csv(out/'Codon_PR2_Neutrality_30mitogenomes.tsv',sep='\t',index=False)
    species=comp.groupby('species',as_index=False)[['GC3','GC12','PR2_A3_A3T3','PR2_G3_G3C3']].mean()
    species.to_csv(out/'species_composition.tsv',sep='\t',index=False)
    x=species.GC3.to_numpy();y=species.GC12.to_numpy();fit=stats.linregress(x,y)
    residuals=y-(fit.intercept+fit.slope*x); crit=stats.t.ppf(.975,len(x)-2)
    diagnostics=dict(n=len(x),slope=fit.slope,intercept=fit.intercept,r_squared=fit.rvalue**2,
                     p_value=fit.pvalue,slope_ci95=[fit.slope-crit*fit.stderr,fit.slope+crit*fit.stderr],
                     shapiro_p=stats.shapiro(residuals).pvalue,
                     abs_residual_GC3_pearson_r=stats.pearsonr(abs(residuals),x).statistic,
                     abs_residual_GC3_pearson_p=stats.pearsonr(abs(residuals),x).pvalue)
    (out/'neutrality_diagnostics.json').write_text(json.dumps(diagnostics,indent=2))
    loo=[]
    for i,name in enumerate(species.species):
        f=stats.linregress(np.delete(x,i),np.delete(y,i))
        loo.append(dict(omitted_species=name,slope=f.slope,intercept=f.intercept,r_squared=f.rvalue**2,p_value=f.pvalue))
    pd.DataFrame(loo).to_csv(out/'neutrality_leave_one_out.tsv',sep='\t',index=False)
    # Read archived estimates; do not treat this as a CODEML rerun.
    omega=pd.read_csv(ref/'dnds_M0_30/dNdS_M0_13genes_final30.csv').set_index('Gene').omega
    gdf=pd.DataFrame(gene_rows);gdf.to_csv(out/'gene_accession_codon_metrics.tsv',sep='\t',index=False)
    summary=[]
    for gene in GENES:
        sub=gdf[gdf.gene==gene];rv=rscu(gene_counts[gene]);group=next((g for g in ['cox','nad','atp'] if gene.startswith(g)),'cytb')
        summary.append(dict(gene=gene,group=group,CDS=len(sub),sense_codons=int(sub.codons.sum()),
                            mean_GC3s=sub.GC3s.mean(),mean_RSCU_AT_end=float(np.nanmean([v for c,v in rv.items() if c[-1] in 'AT'])),omega=float(omega[gene])))
    gs=pd.DataFrame(summary);gs.to_csv(out/'Supplementary_S1.tsv',sep='\t',index=False)
    groups=gs.groupby('group',as_index=False).agg(genes=('gene','size'),CDS=('CDS','sum'),sense_codons=('sense_codons','sum'),
          mean_GC3s=('mean_GC3s','mean'),mean_RSCU_AT_end=('mean_RSCU_AT_end','mean'),mean_omega=('omega','mean'))
    groups.to_csv(out/'Supplementary_S2.tsv',sep='\t',index=False)
    corr=[]
    for variable in ['mean_GC3s','mean_RSCU_AT_end']:
        pear=stats.pearsonr(gs[variable],gs.omega);spear=stats.spearmanr(gs[variable],gs.omega)
        corr.append(dict(variable=variable,n=13,Pearson_r=pear.statistic,Pearson_p=pear.pvalue,Spearman_rho=spear.statistic,Spearman_p=spear.pvalue))
    pd.DataFrame(corr).to_csv(out/'Supplementary_S3.tsv',sep='\t',index=False)
    # Identity-based audit uses full nucleotide sequences, no species-name guessing.
    allrec=[(f,SeqIO.read(f,'genbank')) for f in sorted(raw.glob('*.gb'))]
    chosen={canon(a) for a in manifest.accession}|{'MN793322'};duplicates=[]
    for f,rec in allrec:
        if canon(f.stem) in chosen:continue
        matches=[(ff,rr) for ff,rr in allrec if canon(ff.stem) in chosen and str(rr.seq).upper()==str(rec.seq).upper()]
        for ff,rr in matches:duplicates.append(dict(excluded=f.stem,retained=ff.stem,length=len(rec),species=rec.annotations.get('organism','')))
    pd.DataFrame(duplicates).to_csv(out/'Supplementary_S4.tsv',sep='\t',index=False)
    # Compare reconstructed data with archived outputs using normalized accession versions.
    old=pd.read_csv(ref/'rscu_30.tsv',sep='\t');old['key']=old.accession.map(canon);df['key']=df.accession.map(canon)
    match=df.merge(old,on=['key','codon'],suffixes=('_new','_old'),validate='one_to_one')
    assert len(match)==30*62
    maxdiff=float(abs(match.RSCU_new-match.RSCU_old).max());assert maxdiff<1e-6
    oldmeans=pd.read_csv(ref/'rscu_medias_14_especies.tsv',sep='\t')
    m=means.merge(oldmeans,on=['species','codon'],suffixes=('_new','_old'),validate='one_to_one')
    assert len(m)==14*62 and abs(m.mean_RSCU_new-m.mean_RSCU_old).max()<1e-6
    with (out/'codonw_input_30.fasta').open('w') as f:
        for acc,seq in records:f.write(f'>{acc}\n{seq}\n')
    historical={canon(r.id):str(r.seq) for r in SeqIO.parse(ref/'codonw_input_30.fasta','fasta')}
    assert all(historical[canon(a)]==s for a,s in records)
    checks=dict(CDS_records=len(gene_rows),sense_codons=total,RSCU_records=len(match),RSCU_max_abs_difference=maxdiff,
                species_RSCU_records=len(m),CodonW_input_matches=True,duplicate_pairs=len(duplicates),
                codon_start_values=sorted(set(gdf.codon_start)),Python=platform.python_version(),Biopython=biopython_version,
                external_phylogeny_and_CODEML_rerun=False)
    (out/'validation.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2));print(gs.to_string(index=False));print(diagnostics)
if __name__=='__main__': main()
