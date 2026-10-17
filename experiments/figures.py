from pathlib import Path
import sys,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'src'))
from plane_shadow import *
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
blue='#195A87';orange='#BA4924'
fig,axs=plt.subplots(1,2,figsize=(8.2,2.7),gridspec_kw={'width_ratios':[1.3,1]})
x=np.arange(1,10)
axs[0].plot(x,[0 if m%3 else 405 for m in x],'-o',label=r'$X_a$',color=blue)
axs[0].plot(x,[27 if m%3 else 405 for m in x],'--s',label=r'$Y_a$',color=orange,markersize=4)
axs[0].set(xlabel='Extension degree m',ylabel='Rational planes',xticks=x,yticks=[0,27,405]);axs[0].legend(frameon=False)
axs[1].bar([0,1],[135,126],color=[blue,orange],width=.5,label='Length-3 orbits')
axs[1].bar([0,1],[0,27],bottom=[135,126],color='#ADBCC6',width=.5,label='Fixed planes')
axs[1].set(xticks=[0,1],xticklabels=[r'$X_a$',r'$Y_a$'],ylabel='Number of Frobenius orbits',ylim=(0,180))
axs[1].legend(frameon=False,fontsize=8,loc='upper center')
fig.tight_layout();fig.savefig(R/'paper/figures/plane_descent.pdf');plt.close(fig)
ns=list(class_representatives()); sig=[signature(tuple(i for i,n in enumerate(z) for _ in range(n))) for z in ns]
fig,axs=plt.subplots(2,1,figsize=(8.2,4),sharex=True,gridspec_kw={'height_ratios':[1.15,1]})
x=np.arange(10);width=.25
for j,color in enumerate([blue,'#7199AF',orange]):
 axs[0].bar(x+(j-1)*width,[s[1][j] for s in sig],width,color=color,label=f'b{j}')
axs[0].set(ylabel='Balanced-character multiplicity',ylim=(0,23));axs[0].legend(ncol=3,frameon=False)
for i,z in enumerate(ns):
 cs=tuple(j for j,n in enumerate(z) for _ in range(n))
 axs[1].bar(i,plane_count(cs),color=orange if z in [(1,1,4),(2,2,2)] else blue,width=.6)
axs[1].set(ylabel='Rational planes',xticks=x,xticklabels=[str(z).replace(' ','') for z in ns],xlabel='Cubic-class multiplicities (n0, n1, n2)')
for i in (6,9):
 axs[0].axvspan(i-.44,i+.44,color=orange,alpha=.08);axs[1].axvspan(i-.44,i+.44,color=orange,alpha=.08)
axs[0].text(7.4,18,'Collision: same S = 0\nand (b0,b1,b2) = (8,6,6)',ha='center',fontsize=9)
fig.tight_layout();fig.savefig(R/'paper/figures/twist_census.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(8.2,2.4));ax.axis('off')
boxes=[(.01,.26,.28,.55,'Diagonal equations\n(nonzero coefficients)'),(.37,.58,.28,.35,'Jacobi character sums\nSame zeta function'),(.37,.08,.28,.35,'Exact plane equations\n0 versus 27 over the base field'),(.75,.26,.24,.55,'Two regression fixtures\nfor arithmetic databases')]
from matplotlib.patches import FancyBboxPatch
for x,y,w,h,text in boxes:
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.015',facecolor='#EFF4F7',edgecolor=blue,lw=1))
 ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=9)
for start,end in [((.29,.6),(.36,.74)),((.29,.42),(.36,.26)),((.66,.74),(.74,.6)),((.66,.26),(.74,.42))]:
 ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','color':blue})
fig.tight_layout();fig.savefig(R/'paper/figures/workflow.pdf');plt.close(fig)
