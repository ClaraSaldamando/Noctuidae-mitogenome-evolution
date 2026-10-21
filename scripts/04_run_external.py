#!/usr/bin/env python3
"""Run optional external analyses in new output directories.
Executables must already be installed on PATH. These expensive analyses were
not rerun during package preparation. Run 02_prepare_phylogeny.py first for ML/BI.
"""
import argparse, shutil, subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('analysis',choices=['ml','bi','m0','codonw']);args=p.parse_args()
if args.analysis in ['ml','bi']:
    wd=root/'results/phylogeny'
    if not (wd/'PCG13_31_aligned.fasta').exists():p.error('Run 02_prepare_phylogeny.py first.')
    if args.analysis=='ml':
        if (wd/'ML_new.log').exists():p.error('ML_new already exists; use a fresh working copy.')
        cmd=['iqtree2','-s','PCG13_31_aligned.fasta','-p','PCG13_31_partitions.nex','-m','MFP+MERGE','-B','1000','--alrt','1000','-T','AUTO','-o','MN793322.1','--seed','20260927','--prefix','ML_new']
    else:
        if list(wd.glob('BI_new*')):p.error('BI_new already exists; use a fresh working copy.')
        cmd=['mb','BI_fresh.nex']
    subprocess.run(cmd,cwd=wd,check=True)
elif args.analysis=='m0':
    dst=root/'results/M0_new'
    if dst.exists():p.error('M0_new already exists; use a fresh working copy.')
    shutil.copytree(root/'data/reference/dnds_M0_30',dst)
    for ctl in sorted(dst.glob('*/codeml.ctl')):
        subprocess.run(['codeml','codeml.ctl'],cwd=ctl.parent,check=True)
else:
    wd=root/'results'
    if not (wd/'codonw_input_30.fasta').exists():p.error('Run 01_reproduce_tables.py first.')
    if (wd/'codonw_new.tsv').exists():p.error('codonw_new.tsv already exists.')
    subprocess.run(['codonw','codonw_input_30.fasta','codonw_new.tsv','codonw_new.blk','-code','4','-enc','-gc3s','-machine','-nomenu','-silent','-noblk'],cwd=wd,check=True)
