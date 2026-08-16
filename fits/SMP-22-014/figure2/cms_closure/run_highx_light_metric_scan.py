#!/usr/bin/env python3
"""Test protected high-x light-combination projection metrics at the reference point.

The existing relative-gluon, valence, and F2 terms are retained.  This scan
adds stable relative valence and charge-weighted-F2 terms at 0.05<=x<=0.7 and evaluates the
same fixed-global HERA+CMS likelihood.  It is a projection-metric scan, not a
fit to experimental data.
"""
from __future__ import annotations
import argparse, os, sys
from pathlib import Path
import numpy as np
import yaml
from run_gluon_isolation import write_lhapdf_set
import run_scale_aligned_global_decomposition as global_decomposition

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--source',required=True,type=Path);p.add_argument('--scale-scan',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 here=Path(__file__).resolve().parent; fig,root=here.parent,here.parents[3];sys.path.insert(0,str(fig/'likelihood_scans'));sys.path.insert(0,str(root))
 from scan_tools import parse_likelihood, likelihood_groups, run_xfitter
 from projection_metrics import Figure2ProjectionOperator
 global_decomposition.parse_likelihood=parse_likelihood;global_decomposition.likelihood_groups=likelihood_groups
 from pod_projection.pod_projection import LHAPDF_XGRID
 source,scale,out=a.source.resolve(),a.scale_scan.resolve(),a.output.resolve()
 if out.exists():raise SystemExit(f'Refusing to overwrite {out}')
 out.mkdir(parents=True)
 config=yaml.safe_load((fig/'likelihood_scans/scan_config.yaml').read_text())['projection']; start,stop=config['x_slice'];x=np.asarray(LHAPDF_XGRID[start:stop]);flavours=tuple(int(v) for v in config['flavors'])
 with np.load(source/'runs/_reference/full_pod_projection.npz',allow_pickle=False) as d:target=d['target_grid']
 if target.shape!=(len(flavours),len(x)):raise RuntimeError(f'Unexpected projection grid {target.shape}')
 common=dict(relative_weight=config['relative_weight'],relative_x_range=tuple(config['relative_x_range']),relative_floor=config['relative_floor'],relative_valence_weight=config['relative_valence_weight'],relative_valence_x_range=tuple(config['relative_valence_x_range']),relative_valence_floor=config['relative_valence_floor'],relative_f2_weight=config['relative_f2_weight'],relative_f2_x_range=tuple(config['relative_f2_x_range']),relative_f2_floor=config['relative_f2_floor'],relative_light_x_range=(.05,.99),relative_light_floor=1e-12)
 direct=source/'runs/_reference/direct';grid=scale/'reference/direct_export_direct_reference/output/direct_reference';nuisances=yaml.safe_load((fig/'reference_fit/nuisances.yaml').read_text())['nuisances'];base=parse_likelihood(scale/'global_fixed/reference/matched_direct/output/likelihood.txt');base_totals=global_decomposition.totals(base,nuisances)
 old=os.environ.get('LHAPDF_DATA_PATH','');records=[]
 try:
  for valence_weight, f2_weight in ((.01,0.),(.03,0.),(.1,0.),(.01,.01),(.03,.03)):
   label=f'val_{valence_weight:g}_f2_{f2_weight:g}'.replace('.','d');op=Figure2ProjectionOperator.build(config['basis_set'],100,flavours,x,1.65,'relative_gluon_light',relative_highx_valence_weight=valence_weight,relative_highx_f2_weight=f2_weight,**common);projected,coeff,residual=op.project_grid(target)
   name=f'highx_{label}';write_lhapdf_set(grid,out/'lhapdf'/name,name,{pid:projected[i] for i,pid in enumerate(flavours)},x);os.environ['LHAPDF_DATA_PATH']=str(out/'lhapdf')+(':'+old if old else '')
   like=global_decomposition.evaluate(out/label,direct,root,fig,run_xfitter,global_decomposition.table_parameters(direct,name));tot=global_decomposition.totals(like,nuisances)
   high=x>=.05; metrics={}
   for pid in (21,2,-2,1,-1):
    i=flavours.index(pid);rel=residual[i,high]/np.maximum(abs(target[i,high]),1e-12);metrics[str(pid)]={'rms_relative':float(np.sqrt(np.mean(rel*rel))),'max_absolute_relative':float(np.max(abs(rel)))}
   record={'highx_valence_weight':valence_weight,'highx_f2_weight':f2_weight,'chi2':tot,'delta_chi2':{k:tot[k]-base_totals[k] for k in tot},'metrics':metrics,'coeff_norm':float(np.linalg.norm(coeff))};records.append(record);print(label,record['delta_chi2'])
 finally:
  if old:os.environ['LHAPDF_DATA_PATH']=old
  else:os.environ.pop('LHAPDF_DATA_PATH',None)
 (out/'results.yaml').write_text(yaml.safe_dump({'description':__doc__,'baseline':base_totals,'records':records},sort_keys=False))
if __name__=='__main__':main()
