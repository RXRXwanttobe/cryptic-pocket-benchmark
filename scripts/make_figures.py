"""Publication-style English figures from the frozen, audited results.

No prediction or selection is rerun. All pointwise intervals are Wilson 95% CIs.
The paired difference intervals are the frozen target bootstrap results.
"""
from pathlib import Path
import json,math,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches,colors,ticker

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'docs/figures';OUT.mkdir(exist_ok=True)
S=json.loads((BASE/'data/published/statistics.json').read_text())
DATA=json.loads((BASE/'data/published/target_metrics.json').read_text());ROWS=DATA['targets'];N=len(ROWS)
MODELS=['p2rank','learned','learned_plm']
NAMES={'p2rank':'P2Rank','learned':'Lacuna (learned)','learned_plm':'Lacuna (learned + PLM)'}
COLORS={'p2rank':'#D17C32','learned':'#008C95','learned_plm':'#7752A0'}
MARKERS=dict(zip(MODELS,['o','s','^']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9.5,'axes.titlesize':10.5,
 'legend.fontsize':8.3,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'lines.linewidth':1.6,
 'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
manifest={'n':N,'source_results_sha256':hashlib.sha256((BASE/'data/published/target_metrics.json').read_bytes()).hexdigest(),
 'source_statistics_sha256':hashlib.sha256((BASE/'data/published/statistics.json').read_bytes()).hexdigest(),
 'intervals':'Wilson 95% pointwise binomial intervals; paired differences: existing 10000-resample target bootstrap. CIs are not adjusted for simultaneous coverage.',
 'threshold':.25,'files':[],'plotted_points':{}}
def wilson(h,n=N):
 z=1.959963984540054;p=np.asarray(h,dtype=float)/n;den=1+z*z/n
 c=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return np.stack([c-half,c+half],axis=-1)*100
def finish(fig,name,caption):
 for ext in ['png','svg','pdf']:
  p=OUT/(name+'.'+ext);fig.savefig(p,dpi=600,bbox_inches='tight',pad_inches=.10)
  manifest['files'].append({'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 plt.close(fig);(OUT/(name+'.caption.txt')).write_text(caption+'\n')
def style(ax):
 ax.set_axisbelow(True);ax.grid(axis='y',color='#E6E6E6',linewidth=.55);ax.yaxis.set_major_formatter(ticker.PercentFormatter(100,decimals=0))
def hit(m,k):return S['totals'][m][str(k)]['hits']

# Figure 1: asymmetric intervals, count labels above upper caps.
fig,ax=plt.subplots(figsize=(7.2,4.15));fig.subplots_adjust(top=.79,bottom=.17,left=.10,right=.98)
x=np.arange(3);pts={}
for j,m in enumerate(MODELS):
 y=np.array([hit(m,k) for k in [1,3,5]]);ci=wilson(y);xx=x+(j-1)*.235
 assert np.allclose(ci,np.array([S['totals'][m][str(k)]['ci'] for k in [1,3,5]])*100)
 ax.bar(xx,y,.215,color=COLORS[m],zorder=3,label=NAMES[m])
 ax.errorbar(xx,y,yerr=np.vstack((y-ci[:,0],ci[:,1]-y)),fmt='none',ecolor='#242424',elinewidth=.85,capsize=2.5,capthick=.85,zorder=4)
 for xx0,y0,hi in zip(xx,y,ci[:,1]):ax.text(xx0,hi+2,f'{y0}/100',ha='center',va='bottom',fontsize=8)
 pts[m]={'k':[1,3,5],'hits':y.tolist(),'ci_percent':ci.tolist()}
style(ax);ax.set(xticks=x,xticklabels=['Top-1','Top-3','Top-5'],ylim=(0,100),ylabel='Target recovery (%)',xlabel='Number of top-ranked pockets evaluated')
fig.suptitle('Figure 1  Recovery within the top-ranked pockets',x=.10,y=.98,ha='left',fontsize=11)
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.54,.93),ncol=3,frameon=False,handlelength=1.4,columnspacing=1.4)
fig.text(.10,.015,'n = 100; success: residue IoU >= 0.25. Error bars: Wilson 95% CI.',fontsize=8,color='#454545')
manifest['plotted_points']['01_top135']=pts
finish(fig,'01_top135','Figure 1. Recovery of the annotated site within the top 1, 3 or 5 predictions for the same 100 targets. Error bars show asymmetric Wilson 95% confidence intervals for each recovery proportion, not uncertainty in binding affinity or a paired test of model differences. The primary endpoint is Top-5 at residue IoU >= 0.25.')

# Figures 2 and 9: pointwise (not simultaneous) uncertainty; identical denominators.
for kind in ['k','threshold']:
 fig,ax=plt.subplots(figsize=(7.2,4.0));fig.subplots_adjust(top=.79,bottom=.20,left=.10,right=.98)
 pts={}
 for m in MODELS:
  xx=np.arange(1,11) if kind=='k' else np.array(sorted(map(float,S['thresholds'])))
  yy=np.array([hit(m,int(k)) for k in xx]) if kind=='k' else np.array([S['thresholds'][f'{t:.2f}'][m] for t in xx])
  if kind=='threshold':
   independently=[sum(r['metrics'][m]['top5_best_iou']>=round(float(t),2) for r in ROWS) for t in xx]
   assert independently==yy.tolist(),(m,independently,yy)
  ci=wilson(yy)
  ax.fill_between(xx,ci[:,0],ci[:,1],color=COLORS[m],alpha=.075,linewidth=0,zorder=1)
  ax.plot(xx,yy,marker=MARKERS[m],ms=4,color=COLORS[m],label=NAMES[m],zorder=3)
  pts[m]={'x':xx.tolist(),'hits':yy.tolist(),'ci_percent':ci.tolist()}
 style(ax);ax.set(ylim=(0,100),ylabel='Target recovery (%)' if kind=='k' else 'Top-5 target recovery (%)')
 if kind=='k':
  name='02_topk_curve';title='Figure 2  Recovery as more candidates are reviewed'
  ax.set(xlim=(.8,10.65),xticks=range(1,11),xlabel='Number of top-ranked pockets evaluated (k)')
  ax.axvline(5,color='#777777',ls=(0,(3,3)),lw=.8,zorder=2);ax.text(5.12,93,'Primary k = 5',fontsize=8,color='#555555')
  for m,offset in [('p2rank',0),('learned',-9),('learned_plm',8)]:ax.annotate(str(hit(m,10))+'%',(10,hit(m,10)),xytext=(6,offset),textcoords='offset points',color=COLORS[m],fontsize=8,va='center')
  note='Shading: pointwise Wilson 95% CI; n = 100 at every k.'
 else:
  name='09_iou_sensitivity';title='Figure 9  Sensitivity to the required residue overlap'
  ax.set(xlim=(.09,.51),xticks=xx,xlabel='Minimum residue IoU required for a Top-5 hit')
  ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
  ax.axvline(.25,color='#777777',ls=(0,(3,3)),lw=.8,zorder=2);ax.text(.257,94,'Primary IoU = 0.25',fontsize=8,color='#555555')
  note='Shading: pointwise Wilson 95% CI. Threshold sensitivity is exploratory.'
 fig.suptitle(title,x=.10,y=.98,ha='left',fontsize=11)
 fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.54,.93),ncol=3,frameon=False,handlelength=1.7,columnspacing=1.3)
 fig.text(.10,.025,note,fontsize=8,color='#454545')
 manifest['plotted_points'][name]=pts
 finish(fig,name,note+' Shaded intervals are not simultaneous confidence bands. The same targets are reused at each point and the points are correlated. Model superiority must be assessed using paired differences, not overlap of marginal intervals. At a given k, only actual returned predictions are scored if fewer than k are returned.')

# Figure 3: threshold-centered blue-white-red overlap map.
matrix=np.array([[r['metrics'][m]['top5_best_iou'] for m in MODELS] for r in ROWS])
from matplotlib.colors import TwoSlopeNorm
cmap=plt.get_cmap('RdBu_r');norm=TwoSlopeNorm(vmin=0,vcenter=.25,vmax=1)
fig,axs=plt.subplots(1,2,figsize=(7.4,8.8));fig.subplots_adjust(top=.88,bottom=.085,left=.12,right=.98,wspace=.52)
for block,ax in enumerate(axs):
 sub=matrix[block*50:(block+1)*50]
 im=ax.imshow(sub,aspect='auto',cmap=cmap,norm=norm,interpolation='none')
 ax.set_xticks(range(3),['P2Rank','Lacuna\nlearned','Lacuna\n+ PLM'],fontsize=8);ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=4)
 ax.set_yticks(range(50),[r['id']+' '+r['pdb'] for r in ROWS[block*50:(block+1)*50]],fontsize=7)
 for i in range(50):
  for j in range(3):
   v=sub[i,j];rgb=cmap(norm(v))[:3]
   linear=[q/12.92 if q<=.04045 else ((q+.055)/1.055)**2.4 for q in rgb]
   lum=sum(q*w for q,w in zip(linear,[.2126,.7152,.0722]))
   white_contrast=1.05/(lum+.05);black_contrast=(lum+.05)/.05
   ink='white' if white_contrast>black_contrast else '#111111'
   ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=7,color=ink)
   if v>=.25:ax.add_patch(patches.Rectangle((j-.475,i-.455),.95,.91,fill=False,edgecolor=ink,linewidth=.6))
 for spine in ax.spines.values():spine.set_visible(False)
cax=fig.add_axes([.25,.035,.55,.012]);bar=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,.125,.25,.5,.75,1]);bar.ax.set_xticklabels(['0','0.125','0.25','0.50','0.75','1.00']);bar.ax.tick_params(labelsize=7,length=2,pad=2);bar.set_label('Best residue IoU among Top-5 predictions',fontsize=8,labelpad=2)
fig.suptitle('Figure 3  Target-level pocket overlap',x=.12,y=.982,ha='left',fontsize=11)
fig.text(.12,.946,'Blue: below 0.25; red: above 0.25; outlines: IoU >= 0.25.',fontsize=8)
fig.text(.12,.923,'Piecewise color scale centered on 0.25; displayed values are rounded only.',fontsize=7.5,color='#454545')

manifest['plotted_points']['03_target_heatmap']={'ids':[r['id'] for r in ROWS],'pdbs':[r['pdb'] for r in ROWS],'models':MODELS,'iou':matrix.tolist(),'hits':(matrix>=.25).sum(axis=0).tolist()}
manifest['figure3_color_scale']={'colormap':'RdBu_r','vmin':0,'vcenter':.25,'vmax':1,'normalization':'piecewise linear centered on primary threshold'}
finish(fig,'03_target_heatmap','Best residue IoU among Top-5 predictions. Blue below 0.25, white at the threshold, red above it. Contrasting outlines mark IoU >= 0.25. The color scale is piecewise linear with different slopes either side of 0.25. Values are rounded for display only. No single-target error bars are inferred.')

# Recommended supplemental figures: paired uncertainty and ranking ceiling.
fig,ax=plt.subplots(figsize=(7.2,3.0));fig.subplots_adjust(top=.77,bottom=.28,left=.23,right=.98)
comparisons=[]
for i,m in enumerate(['learned','p2rank']):
 c=next(c for c in S['comparisons'] if c['a']=='learned_plm' and c['b']==m and c['k']==5);comparisons.append(c)
 d=100*c['difference'];lo,hi=np.array(c['difference_ci'])*100
 ax.errorbar(d,i,xerr=[[d-lo],[hi-d]],fmt='o',ms=6,color=COLORS['learned_plm'],elinewidth=1.5,capsize=4,capthick=1)
 ax.text(20,i,f'{d:+.0f} pp  [{lo:+.0f}, {hi:+.0f}]\nHolm p = {c["holm_adjusted_p"]:.3f}',va='center',fontsize=8)
ax.axvline(0,color='#666666',ls=(0,(3,3)),lw=.8);ax.set(yticks=[0,1],yticklabels=['PLM vs learned','PLM vs P2Rank'],xlim=(-7,35),ylim=(1.55,-.55),xlabel='Paired Top-5 recovery difference (percentage points)')
fig.suptitle('Figure 10  Uncertainty in the observed improvement',x=.08,y=.98,ha='left',fontsize=11)
fig.text(.08,.04,'Intervals: target-paired bootstrap 95% CI (10,000 resamples).\nTwo planned tests: exact McNemar with Holm adjustment. CIs are per-comparison.',fontsize=8,color='#454545')
finish(fig,'10_paired_uncertainty','Recommended supplement. Paired target bootstrap confidence intervals for the Top-5 difference. Exact McNemar tests have Holm adjustment for two planned comparisons; displayed intervals are per-comparison, not family-wise confidence intervals. A marginal CI overlap in Figure 1 is not a test of these differences.')

fig,ax=plt.subplots(figsize=(7.2,3.8));fig.subplots_adjust(top=.81,bottom=.22,left=.10,right=.98)
ceiling=np.array([sum(r['metrics'][m]['oracle_any_hit'] for r in ROWS) for m in MODELS]);top=np.array([hit(m,5) for m in MODELS]);x=np.arange(3)
ax.bar(x,ceiling,.53,color='#E3E6E8',zorder=2)
for i,m in enumerate(MODELS):
 ax.bar(i,top[i],.53,color=COLORS[m],zorder=3)
 ax.text(i,ceiling[i]+2,f'Any candidate: {ceiling[i]}/100',ha='center',fontsize=8)
 ax.text(i,top[i]-8,f'Top-5: {top[i]}/100',ha='center',va='center',color='white',fontsize=9)
style(ax);ax.set(xticks=x,xticklabels=['P2Rank','Lacuna\nlearned','Lacuna\nlearned + PLM'],ylim=(0,100),ylabel='Target recovery (%)')
fig.suptitle('Figure 12  Detection coverage and ranking conversion',x=.10,y=.98,ha='left',fontsize=11)
fig.text(.10,.88,'Gray: retrospective upper bound over all returned candidates.',fontsize=8,color='#454545')
fig.text(.10,.015,'The two Lacuna rankers share the same candidate pool.\nCandidate counts differ between P2Rank and Lacuna; the gray bars are not budget-matched.',fontsize=8,color='#454545')
finish(fig,'12_candidate_ceiling','Recommended mechanistic supplement. Gray bars show the retrospective oracle coverage using every returned candidate. Colored bars show actual Top-5 recovery. Both Lacuna configurations share exactly the same complete candidate pool. All-candidate counts differ between tools; no equal-budget superiority claim is made for the gray bars. This deterministic oracle decomposition is not an independent experimental estimate of attainable deployment performance.')

(OUT/'figure-manifest.json').write_text(json.dumps(manifest,indent=2))
(OUT/'README.txt').write_text('English publication-style figures\nPNG: 600 dpi; SVG: editable text; PDF: vector.\nModels: P2Rank 2.5.1 default, Lacuna 1.2.0 learned and learned-plm.\nFrozen target selection and numerical results are unchanged.\nFigures 1, 2, 3 and 9: requested report figures.\nFigures 10 and 12: recommended supplements, not automatically added to the report.\nSee individual captions and manifest for uncertainty definitions and source hashes.\nNMA conformers are structural samples, not a molecular-dynamics trajectory.\n')
print('Created six English figures, 18 image/vector files; all plotted values verified.',flush=True)
