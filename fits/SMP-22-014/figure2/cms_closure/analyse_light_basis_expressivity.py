#!/usr/bin/env python3
"""Diagnose whether high-x light combinations are missing from the POD span.

Compares the configured projection with the best possible least-squares
reconstruction of high-x uv, dv, and F2 combinations using all 100 modes.
The latter has no gluon or global-PDF protection, so its collateral gluon
error distinguishes lack of span from an inter-flavour trade-off.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import yaml

def rms(x): return float(np.sqrt(np.mean(np.asarray(x)**2)))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path);a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 here=Path(__file__).resolve().parent;fig,root=here.parent,here.parents[3];sys.path[:0]=[str(fig/'likelihood_scans'),str(root)]
 from projection_metrics import Figure2ProjectionOperator
 from pod_projection.pod_projection import LHAPDF_XGRID
 cfg=yaml.safe_load((fig/'likelihood_scans/scan_config.yaml').read_text())['projection'];start,stop=cfg['x_slice'];x=np.asarray(LHAPDF_XGRID[start:stop]);flavours=tuple(int(v) for v in cfg['flavors'])
 with np.load(a.source/'runs/_reference/full_pod_projection.npz',allow_pickle=False) as d: target,configured=d['target_grid'],d['projected_grid']
 op=Figure2ProjectionOperator.build(cfg['basis_set'],100,flavours,x,1.65,cfg['metric'],relative_weight=cfg['relative_weight'],relative_x_range=tuple(cfg['relative_x_range']),relative_floor=cfg['relative_floor'],relative_valence_weight=cfg['relative_valence_weight'],relative_valence_x_range=tuple(cfg['relative_valence_x_range']),relative_valence_floor=cfg['relative_valence_floor'],relative_f2_weight=cfg['relative_f2_weight'],relative_f2_x_range=tuple(cfg['relative_f2_x_range']),relative_f2_floor=cfg['relative_f2_floor'])
 reference=op.base.reference_grid; shifts=op.base.matrix.reshape(len(flavours),len(x),-1); disp=target-reference
 ix={pid:flavours.index(pid) for pid in flavours}; mask=(x>=.05)&(x<=.7); charges={2:4/9,-2:4/9,1:1/9,-1:1/9,3:1/9,-3:1/9,4:4/9,-4:4/9,5:1/9,-5:1/9}
 def combo(values): return np.stack((values[ix[2]]-values[ix[-2]],values[ix[1]]-values[ix[-1]],sum(q*values[ix[p]] for p,q in charges.items())))
 base_combo=combo(reference); target_combo=combo(target)
 shift_combo=np.stack((shifts[ix[2]]-shifts[ix[-2]], shifts[ix[1]]-shifts[ix[-1]], sum(q*shifts[ix[p]] for p,q in charges.items()))) # combo,x,mode
 denom=np.maximum(abs(target_combo[:,mask]),1e-10); A=(shift_combo[:,mask,:]/denom[:,:,None]).reshape(-1,100); b=((target_combo-base_combo)[:,mask]/denom).reshape(-1)
 coeff,*_=np.linalg.lstsq(A,b,rcond=None); light_only=reference+np.einsum('m,fxm->fx',coeff,shifts)
 def metrics(grid):
  comb=combo(grid);out={name:rms((comb[i,mask]-target_combo[i,mask])/np.maximum(abs(target_combo[i,mask]),1e-10)) for i,name in enumerate(('uv','dv','F2'))};g=ix[21];out['gluon']=rms((grid[g,mask]-target[g,mask])/np.maximum(abs(target[g,mask]),1e-10));return out
 s=np.linalg.svd(A,compute_uv=False); result={'configured':metrics(configured),'light_only_best_possible':metrics(light_only),'light_combination_matrix_rank':int(np.linalg.matrix_rank(A)),'condition_number':float(s[0]/s[-1]),'coefficient_norm':float(np.linalg.norm(coeff))}
 (a.output_dir/'expressivity.yaml').write_text(yaml.safe_dump(result,sort_keys=False))
 fig_,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True);axes[0].semilogy(s/s[0],'o-');axes[0].set(title='High-x light-combination singular spectrum',xlabel='mode index',ylabel=r'$s_i/s_0$');axes[0].grid(alpha=.25)
 names=('uv','dv','F2','gluon');y=np.arange(4);axes[1].barh(y-.18,[result['configured'][n] for n in names],.36,label='configured');axes[1].barh(y+.18,[result['light_only_best_possible'][n] for n in names],.36,label='light-only optimum');axes[1].set_yticks(y,names);axes[1].set_xscale('log');axes[1].set(title=r'RMS relative closure, $0.05\leq x\leq0.7$');axes[1].legend(frameon=False);axes[1].grid(axis='x',alpha=.25);fig_.savefig(a.output_dir/'light_basis_expressivity.png',dpi=180)
 print(yaml.safe_dump(result,sort_keys=False))
if __name__=='__main__':main()
