#!/usr/bin/env python3
"""Plot scale-aligned global flavour attribution and CMS jet-bin residuals."""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import yaml

LABEL={21:'g',2:'u',-2:r'$\bar u$',1:'d',-1:r'$\bar d$',3:'s',-3:r'$\bar s$',4:'c',-4:r'$\bar c$',5:'b',-5:r'$\bar b$'}
RAPIDITY=(r'$|y|<0.5$',r'$0.5<|y|<1.0$',r'$1.0<|y|<1.5$',r'$1.5<|y|<2.0$')
def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--summary',required=True,type=Path); p.add_argument('--bins',required=True,type=Path); p.add_argument('--output-dir',required=True,type=Path); a=p.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
 summary=yaml.safe_load(a.summary.read_text()); shifts=summary['shifts']; pids=[int(x) for x in summary['flavours']]
 labels=[LABEL[x] for x in pids]+['all']; cases=[f'flavour_{x:+d}' for x in pids]+['all_flavours']; groups=('HERA','CMS','sum')
 fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True); y=np.arange(len(cases)); width=.24
 for i,g in enumerate(groups): ax.barh(y+(i-1)*width,[shifts[c][g] for c in cases],height=width,label=g)
 ax.axvline(0,color='.35',lw=.8); ax.set_yticks(y,labels); ax.set_xlabel(r'hybrid minus matched-direct $\chi^2$'); ax.set_title('Scale-aligned full-POD residual attributed by input flavour'); ax.legend(frameon=False); ax.grid(axis='x',alpha=.2); fig.savefig(a.output_dir/'flavour_likelihood_attribution.png',dpi=180); plt.close(fig)
 with np.load(a.bins,allow_pickle=False) as d:
  direct={k.removeprefix('direct_input_'):d[k] for k in d.files if k.startswith('direct_input_')}; allf={k.removeprefix('all_flavours_'):d[k] for k in d.files if k.startswith('all_flavours_')}
  rows=[]
  for pid in pids:
   theory=d[f'flavour_{pid:+d}_theory']; rows.append(100*(theory/direct['theory']-1))
 fig,axes=plt.subplots(2,4,figsize=(15,6.4),sharex='col',gridspec_kw={'height_ratios':[2.2,1]},constrained_layout=True)
 for i,(top,bot) in enumerate(zip(axes[0],axes[1]),1):
  m=direct['rapidity_bin']==i; pt=direct['pt'][m]; top.plot(pt,direct['theory'][m],'o-',ms=3,label='matched direct'); top.plot(pt,allf['theory'][m],'s--',ms=3,label='all-flavour POD hybrid'); top.set_xscale('log'); top.set_yscale('log'); top.set_title(RAPIDITY[i-1]); top.grid(alpha=.2); bot.plot(pt,100*(allf['theory'][m]/direct['theory'][m]-1),'D-',color='#6a3d9a',ms=3); bot.axhline(0,color='.35',lw=.8); bot.set_xscale('log'); bot.set_xlabel(r'$p_T$ [GeV]'); bot.set_ylabel('POD/direct−1 [%]'); bot.grid(alpha=.2)
 axes[0,0].set_ylabel('theory'); axes[0,0].legend(frameon=False,fontsize=8); fig.suptitle('CMS inclusive jets: scale-aligned all-flavour POD bin closure'); fig.savefig(a.output_dir/'cms_bin_closure.png',dpi=180); plt.close(fig)
 fig,ax=plt.subplots(figsize=(15,5),constrained_layout=True); image=ax.imshow(np.asarray(rows),aspect='auto',cmap='coolwarm',vmin=-.35,vmax=.35); ax.set_yticks(np.arange(len(pids)),[LABEL[x] for x in pids]); ax.set_xlabel('CMS bin (four rapidity blocks, increasing $p_T$)'); ax.set_ylabel('one replaced input flavour'); ax.set_title('Individual flavour-induced CMS theory shifts [%]'); [ax.axvline(v-.5,color='k',lw=.6,alpha=.35) for v in np.where(np.diff(direct['rapidity_bin'])!=0)[0]+1]; fig.colorbar(image,ax=ax,label='hybrid/direct−1 [%]'); fig.savefig(a.output_dir/'cms_flavour_bin_response.png',dpi=180); plt.close(fig)
 print(f'Wrote figures in {a.output_dir}')
if __name__=='__main__': main()
