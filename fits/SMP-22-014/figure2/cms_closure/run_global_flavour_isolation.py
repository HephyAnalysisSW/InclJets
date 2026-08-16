#!/usr/bin/env python3
"""Attribute the remaining scale-aligned global POD residual flavour by flavour.

At the ADbar reference point, every hybrid uses the same direct LHAPDF input
table at Q0=1.65 GeV and the original 11-dataset, 198-fixed-nuisance global
likelihood.  A hybrid differs only by replacing the selected input flavour(s)
with their full-POD projection.  No fit, profiling, or nuisance adjustment is
performed.
"""
from __future__ import annotations
import argparse
import os
import shutil
import sys
from pathlib import Path
import numpy as np
import yaml
from run_gluon_isolation import parse_rows, write_lhapdf_set
import run_scale_aligned_global_decomposition as global_decomposition

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--source',required=True,type=Path); p.add_argument('--scale-scan',required=True,type=Path); p.add_argument('--global-scan',required=True,type=Path); p.add_argument('--output',required=True,type=Path); a=p.parse_args()
    here=Path(__file__).resolve().parent; figure_dir,project_root=here.parent,here.parents[3]
    sys.path.insert(0,str(figure_dir/'likelihood_scans'))
    from scan_tools import parse_likelihood, likelihood_groups, run_xfitter
    global_decomposition.parse_likelihood = parse_likelihood
    global_decomposition.likelihood_groups = likelihood_groups
    source,scale,global_scan,out=a.source.resolve(),a.scale_scan.resolve(),a.global_scan.resolve(),a.output.resolve()
    if out.exists(): raise SystemExit(f'Refusing to overwrite {out}')
    out.mkdir(parents=True)
    direct_source=source/'runs/_reference/direct'; projection=source/'runs/_reference/full_pod_projection.npz'
    with np.load(projection,allow_pickle=False) as d: flavours,x,projected=d['flavors'].astype(int),d['x_grid'],d['projected_grid']
    replacement={int(pid):projected[i] for i,pid in enumerate(flavours)}
    direct_grid=scale/'reference/direct_export_direct_reference/output/direct_reference'
    if not direct_grid.is_dir(): raise RuntimeError(f'Missing direct grid {direct_grid}')
    lhapdf=out/'lhapdf'; cases=[('direct_input',{})]+[(f'flavour_{pid:+d}',{int(pid):replacement[int(pid)]}) for pid in flavours]+[('all_flavours',replacement)]
    for label,replaced in cases: write_lhapdf_set(direct_grid,lhapdf/f'global_{label}',f'global_{label}',replaced or None,x)
    nuisances=yaml.safe_load((figure_dir/'reference_fit/nuisances.yaml').read_text())['nuisances']
    direct_like=parse_likelihood(global_scan/'reference/matched_direct/output/likelihood.txt'); direct_rows=parse_rows(global_scan/'reference/matched_direct/output/fittedresults.txt')
    likelihoods={'direct_input':direct_like}; rows={'direct_input':direct_rows}
    old=os.environ.get('LHAPDF_DATA_PATH','')
    try:
      os.environ['LHAPDF_DATA_PATH']=str(lhapdf)+(':'+old if old else '')
      for number,(label,_) in enumerate(cases[1:],1):
        target=out/label; like=global_decomposition.evaluate(target,direct_source,project_root,figure_dir,run_xfitter,global_decomposition.table_parameters(direct_source,f'global_{label}'))
        if like['free_parameter_count']!=0 or like['nuisance_treatment']!='fixed' or like['nuisance_count']!=198: raise RuntimeError(f'Fixed-global contract failed for {label}: {like}')
        likelihoods[label]=like; rows[label]=parse_rows(target/'output/fittedresults.txt')
        delta=global_decomposition.totals(like,nuisances); base=global_decomposition.totals(direct_like,nuisances)
        print(f'[{number}/{len(cases)-1}] {label}: HERA {delta["HERA"]-base["HERA"]:+.6f}, CMS {delta["CMS"]-base["CMS"]:+.6f}')
    finally:
      if old: os.environ['LHAPDF_DATA_PATH']=old
      else: os.environ.pop('LHAPDF_DATA_PATH',None)
    if any(not np.array_equal(direct_rows['pt'],record['pt']) for record in rows.values()): raise RuntimeError('CMS rows do not align')
    np.savez_compressed(out/'flavour_isolation.npz',**{f'{label}_{key}':value for label,record in rows.items() for key,value in record.items()})
    base=global_decomposition.totals(direct_like,nuisances); shifts={label:{key:global_decomposition.totals(like,nuisances)[key]-base[key] for key in ('HERA','CMS','sum')} for label,like in likelihoods.items() if label!='direct_input'}
    (out/'summary.yaml').write_text(yaml.safe_dump({'description':__doc__,'direct_input':base,'shifts':shifts,'flavours':[int(x) for x in flavours]},sort_keys=False))

if __name__=='__main__': main()
