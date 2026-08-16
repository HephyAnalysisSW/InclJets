#!/usr/bin/env python3
"""Compare original and scale-aligned HERA/CMS ADbar likelihood decompositions."""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

GROUPS=(("HERA","HERA"),("CMS","CMS inclusive jets"),("sum","HERA + CMS"))

def original(path: Path):
    with np.load(path,allow_pickle=False) as d:
        mask=d["parameter_name"].astype(str)=="ADbar"; order=np.argsort(d["scan_coordinate_sigma"][mask])
        return (d["scan_coordinate_sigma"][mask][order], {f"direct_{g}":d[f"direct_chi2_{g}"][mask][order] for g,_ in GROUPS}|{f"pod_{g}":d[f"full_pod_chi2_{g}"][mask][order] for g,_ in GROUPS})

def aligned(path: Path):
    with np.load(path,allow_pickle=False) as d:
        x=d["ADbar"]; order=np.argsort(x); x=x[order]; zero=np.argmin(abs(x-0.267003)); sigma=(x-x[zero])/((x[-1]-x[0])/4)
        return sigma, {f"direct_{g}":d[f"matched_direct_{g}"][order] for g,_ in GROUPS}|{f"pod_{g}":d[f"full_pod_{g}"][order] for g,_ in GROUPS}

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--before",required=True,type=Path); p.add_argument("--after",required=True,type=Path); p.add_argument("--output",required=True,type=Path); a=p.parse_args()
    routes=("before: original direct/POD route",original(a.before)),("after: common table route at $Q_0=1.65$ GeV",aligned(a.after))
    fig,axes=plt.subplots(4,3,figsize=(13.2,10),sharex="col",gridspec_kw={"height_ratios":[2.1,1,2.1,1]},constrained_layout=True)
    for ir,(label,(x,data)) in enumerate(routes):
        zero=np.argmin(abs(x)); top_row=2*ir
        for column,(key,title) in enumerate(GROUPS):
            direct,pod=data[f"direct_{key}"],data[f"pod_{key}"]; base=direct[zero]
            top,bot=axes[top_row,column],axes[top_row+1,column]
            top.plot(x,direct-base,"o-",lw=1.8,ms=4.5,label="direct")
            top.plot(x,pod-base,"s--",lw=1.7,ms=4.2,label="full 100-mode POD")
            top.axhline(0,color=".72",lw=.8); top.axvline(0,color=".82",lw=.8); top.grid(alpha=.2)
            top.set_title(title if ir==0 else ""); top.set_ylabel(r"$\chi^2-\chi^2_{\rm direct}(0)$")
            bot.plot(x,pod-direct,"D-",color="#7a3e9d",lw=1.6,ms=4)
            bot.axhline(0,color=".35",lw=.9); bot.axvline(0,color=".82",lw=.8); bot.grid(alpha=.2)
            bot.set_ylabel(r"$\chi^2_{\rm POD}-\chi^2_{\rm direct}$")
            bot.set_xlabel(r"ADbar scan coordinate [$\sigma_{\rm HESSE}$]")
        axes[top_row,0].legend(frameon=False,loc="best")
        axes[top_row,0].annotate(label,xy=(0,0.5),xycoords="axes fraction",xytext=(-54,0),textcoords="offset points",rotation=90,ha="center",va="center",fontsize=11,fontweight="bold")
    fig.suptitle("Unprofiled ADbar likelihood decomposition: effect of scale-aligned POD closure",fontsize=15)
    a.output.parent.mkdir(parents=True,exist_ok=True); fig.savefig(a.output,dpi=180); print(f"Wrote {a.output}")
if __name__=="__main__": main()
