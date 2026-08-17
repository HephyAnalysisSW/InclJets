#!/usr/bin/env python3
"""Insert the shared scale-aligned reference into each scan direction and plot."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import yaml

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--scan',required=True,type=Path);p.add_argument('--source',required=True,type=Path);p.add_argument('--reference',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 here=Path(__file__).resolve().parent;fig=here.parent;sys.path.insert(0,str(fig/'likelihood_scans'));from scan_tools import parse_likelihood,likelihood_groups
 nuis=yaml.safe_load((fig/'reference_fit/nuisances.yaml').read_text())['nuisances']; direct=likelihood_groups(parse_likelihood(a.reference/'matched_direct/output/likelihood.txt'),nuis);pod=likelihood_groups(parse_likelihood(a.reference/'full_pod/output/likelihood.txt'),nuis)
 def total(x):return {k:float(x[k]['total_chi2']) for k in ('HERA','CMS')}|{'sum':float(x['HERA']['total_chi2']+x['CMS']['total_chi2'])}
 direct,pod=total(direct),total(pod)
 with np.load(a.scan,allow_pickle=False) as d: arrays={k:np.array(d[k]) for k in d.files}
 names=list(dict.fromkeys(arrays['parameter_name'].astype(str))); ref=yaml.safe_load((a.source/'runs/_reference/direct/evaluation.yaml').read_text())['parameter_values']
 arrays['parameter_name']=np.concatenate((arrays['parameter_name'],np.asarray(names)))
 arrays['scan_coordinate_sigma']=np.concatenate((arrays['scan_coordinate_sigma'],np.zeros(len(names))))
 arrays['parameter_value']=np.concatenate((arrays['parameter_value'],np.asarray([float(ref[n]) for n in names])))
 for route,values in (('direct',direct),('full_pod',pod)):
  for group,value in values.items(): arrays[f'{route}_chi2_{group}']=np.concatenate((arrays[f'{route}_chi2_{group}'],np.full(len(names),value)))
 order=np.lexsort((arrays['scan_coordinate_sigma'],arrays['parameter_name']));arrays={k:v[order] for k,v in arrays.items()};a.output.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(a.output,**arrays);print(f'Wrote {a.output}')
if __name__=='__main__':main()
