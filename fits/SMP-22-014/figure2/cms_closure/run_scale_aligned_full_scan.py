#!/usr/bin/env python3
"""Resumable scale-aligned HERA+CMS scan for every stored PDF parameter point.

Each point exports its native direct PDF to a full LHAPDF grid at Q0=1.65 GeV,
then evaluates matched-direct and stored full-POD PDFs through the same QCDNUM
route, 11 datasets, and frozen 198 nuisance shifts. No minimization/profile.
"""
from __future__ import annotations
import argparse, os, sys, shutil
from pathlib import Path
import numpy as np, yaml
from run_gluon_isolation import export_full_direct_grid, write_lhapdf_set
import run_scale_aligned_global_decomposition as gd

def points(source):
 out=[]
 for p in source.glob('runs/*/sigma_*/direct/evaluation.yaml'):
  root=p.parent.parent; ev=yaml.safe_load(p.read_text()); name=next(iter(ev['parameter_overrides']));
  tag=root.name; coord=float(tag.replace('sigma_p','').replace('sigma_m','-').replace('d','.'))
  out.append((f'{name}_{tag}',root,name,coord))
 return out
def main():
 a=argparse.ArgumentParser(description=__doc__);a.add_argument('--source',required=True,type=Path);a.add_argument('--output',required=True,type=Path);args=a.parse_args()
 here=Path(__file__).resolve().parent;fig,project=here.parent,here.parents[3];sys.path.insert(0,str(fig/'likelihood_scans'))
 from scan_tools import parse_likelihood,likelihood_groups,run_xfitter
 gd.parse_likelihood=parse_likelihood;gd.likelihood_groups=likelihood_groups
 source,args.output=args.source.resolve(),args.output.resolve();args.output.mkdir(parents=True,exist_ok=True)
 nuis=yaml.safe_load((fig/'reference_fit/nuisances.yaml').read_text())['nuisances']; old=os.environ.get('LHAPDF_DATA_PATH',''); rows=[]; todo=points(source)
 try:
  for n,(label,root,name,coord) in enumerate(todo,1):
   point=args.output/label; saved=point/'result.yaml'
   if saved.is_file(): rows.append(yaml.safe_load(saved.read_text()));print(f'[{n}/{len(todo)}] cached {label}');continue
   point.mkdir(parents=True,exist_ok=True); direct,pod=root/'direct',root/'full_pod'; setname=f'aligned_{label}'; grid=point/f'direct_export_{setname}'/'output'/setname
   if not grid.is_dir(): grid=export_full_direct_grid(direct,point,project,fig,run_xfitter,1.65,setname)
   pdfroot=point/'lhapdf';write_lhapdf_set(grid,pdfroot/setname,setname);os.environ['LHAPDF_DATA_PATH']=str(pdfroot)+(':'+old if old else '')
   dl=gd.evaluate(point/'matched_direct',direct,project,fig,run_xfitter,gd.table_parameters(direct,setname));pl=gd.evaluate(point/'full_pod',pod,project,fig,run_xfitter)
   d,p=gd.totals(dl,nuis),gd.totals(pl,nuis); ev=yaml.safe_load((direct/'evaluation.yaml').read_text());r={'parameter_name':name,'scan_coordinate_sigma':coord,'parameter_value':float(ev['parameter_values'][name]),'direct':d,'pod':p};saved.write_text(yaml.safe_dump(r,sort_keys=False));rows.append(r);print(f'[{n}/{len(todo)}] {label}: delta={p["sum"]-d["sum"]:+.4f}')
 finally:
  if old:os.environ['LHAPDF_DATA_PATH']=old
  else:os.environ.pop('LHAPDF_DATA_PATH',None)
 rows.sort(key=lambda r:(r['parameter_name'],r['scan_coordinate_sigma'])); arr={'parameter_name':np.asarray([r['parameter_name'] for r in rows]),'scan_coordinate_sigma':np.asarray([r['scan_coordinate_sigma'] for r in rows]),'parameter_value':np.asarray([r['parameter_value'] for r in rows])}
 for route,key in (('direct','direct'),('full_pod','pod')):
  for g in ('HERA','CMS','sum'):arr[f'{route}_chi2_{g}']=np.asarray([r[key][g] for r in rows])
 np.savez_compressed(args.output/'scan_results.npz',**arr)
if __name__=='__main__':main()
