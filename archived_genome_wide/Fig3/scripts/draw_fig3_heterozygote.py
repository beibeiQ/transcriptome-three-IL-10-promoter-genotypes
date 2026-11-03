"""Plot comprehensive heterozygote-associated patterns from archived contrasts."""
from pathlib import Path
import sys,argparse,json
sys.path.insert(0,str(Path(__file__).resolve().parent/'python_packages'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox

ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
out=ap.parse_args().output;tab=out/'source_tables';fd=out/'figures';fd.mkdir(exist_ok=True)
def read(name):return pd.read_csv(tab/name,dtype={'gene_id':str})
w=read('A_all_genes_paired_comparisons.csv.gz');counts=read('B_pattern_counts.csv')
totals=read('B_totals_and_denominators.csv');selected=read('C_selected_response_examples.csv')
points=read('C_individual_expression_points.csv');effects=read('D_selected_effects_and_intervals.csv')
all_effects=read('supplement_selected_examples_all_stages.csv')
T=['Spleen','Intestine'];stages=['Baseline','Response_d1','Response_d5']
P=['Both_same','Both_opposite','Del_only','Ins_only']
pc={'Both_same':'#7a659e','Both_opposite':'#c96d92','Del_only':'#287d8e','Ins_only':'#d79533'}
short={'Both_same':'Both, same signs','Both_opposite':'Both, opposite signs','Del_only':'H–Del only','Ins_only':'H–Ins only'}
long={'Both_same':'Both detected, same signs','Both_opposite':'Both detected, opposite signs',
      'Del_only':'Only Hetero–HomoDel detected','Ins_only':'Only Hetero–HomoIns detected'}
G=['HomoDel','Hetero','HomoIns'];gc={'HomoDel':'#397e85','Hetero':'#825ca6','HomoIns':'#d58a33'}
gm={'HomoDel':'o','Hetero':'^','HomoIns':'s'}
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':7.5,'axes.titlesize':7.5,'axes.labelsize':7,'xtick.labelsize':6,
    'ytick.labelsize':6,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','legend.frameon':False,
    'text.color':'#222a30','axes.labelcolor':'#222a30','xtick.color':'#30383f','ytick.color':'#30383f'})
W,H=210/25.4,340/25.4
fig=plt.figure(figsize=(W,H),facecolor='white')
def text(x,y,s,size=7.5,**kw):return fig.text(x,y,s,fontsize=size,**kw)
def heading(letter,title,y):
    text(.025,y,letter,12,fontweight='bold',va='top')
    text(.075,y-.001,title,9.2,fontweight='bold',va='top')
def stage_label(s):return {'Baseline':'Baseline','Response_d1':'Day 1 response','Response_d5':'Day 5 response'}[s]
def day_label(s):return 'day 1' if s=='Response_d1' else 'day 5'
ph=[Line2D([],[],ls='',marker='D' if p=='Both_opposite' else 'o',color=pc[p],ms=4,label=long[p]) for p in P]
gh=[Line2D([],[],color=gc[g],marker=gm[g],ms=3,lw=1,label=g) for g in G]

text(.5,.990,'Heterozygote-associated expression and response patterns',10.8,ha='center',va='top',fontweight='bold')
heading('A','Complete comparisons with both homozygous genotypes',.967)
fig.legend(handles=ph,loc='upper center',bbox_to_anchor=(.54,.950),ncol=2,
           fontsize=6.5,handletextpad=.5,columnspacing=1.6,handlelength=1)
for ri,tissue in enumerate(T):
    for ci,stage in enumerate(stages):
        ax=fig.add_axes([.100+ci*.300,.798-ri*.146,.238,.106])
        d=w[(w.tissue==tissue)&(w.stage==stage)]
        finite=np.isfinite(d.effect_D)&np.isfinite(d.effect_I)
        gray=d[finite&d.pattern.eq('Neither_detected')]
        ax.scatter(gray.effect_D,gray.effect_I,s=1.4,c='#d7dcdf',alpha=.58,lw=0,rasterized=True,zorder=1)
        missing=d[finite&d.pattern.eq('Incomplete')]
        ax.scatter(missing.effect_D,missing.effect_I,s=7,c='#8b9399',marker='x',lw=.45,rasterized=True,zorder=2)
        for p in P:
            z=d[finite&d.pattern.eq(p)]
            ax.scatter(z.effect_D,z.effect_I,s=5.2 if p!='Both_opposite' else 10,
                       c=pc[p],marker='D' if p=='Both_opposite' else 'o',alpha=.82,
                       lw=0,rasterized=True,zorder=3 if p!='Both_opposite' else 4)
        ax.set_xscale('asinh',linear_width=1);ax.set_yscale('asinh',linear_width=1)
        extent=max(2,float(np.nanmax(np.abs(d.loc[finite,['effect_D','effect_I']].to_numpy()))))*1.04
        ax.set_xlim(-extent,extent);ax.set_ylim(-extent,extent)
        ticks=[v for v in [-64,-16,-4,-1,0,1,4,16,64] if abs(v)<extent]
        if len(ticks)>7:ticks=[v for v in [-64,-16,-4,0,4,16,64] if abs(v)<extent]
        ax.set_xticks(ticks,ticks);ax.set_yticks(ticks,ticks)
        ax.axhline(0,color='#aeb5b9',lw=.5,ls='--',zorder=0)
        ax.axvline(0,color='#aeb5b9',lw=.5,ls='--',zorder=0)
        ax.tick_params(length=2,pad=2)
        n=int(((d.pattern.isin(P))&d.fully_testable).sum())
        ax.text(.03,.96,f'Any comparison detected: {n:,}',transform=ax.transAxes,
                fontsize=5.9,ha='left',va='top',bbox={'facecolor':'white','edgecolor':'none','alpha':.83,'pad':1.3})
        if ri==0:ax.set_title(stage_label(stage),pad=5,fontsize=8,fontweight='bold')
        if ci==0:ax.set_ylabel(f'{tissue}\nHetero − HomoIns',fontsize=7,labelpad=3)
        if ri==1:ax.set_xlabel('Hetero − HomoDel',labelpad=2)
text(.5,.623,'Effects in log2 units; gray = neither detected, × = incomplete tests; signed asinh spacing',6.3,ha='center')

heading('B','Distribution of detected comparison patterns',.603)
for ti,tissue in enumerate(T):
    ax=fig.add_axes([.105+ti*.49,.518,.355,.057])
    bottom=np.zeros(3)
    for p in P:
        values=np.array([int(counts[(counts.tissue==tissue)&(counts.stage==s)&(counts.pattern==p)].genes.sum()) for s in stages])
        ax.bar(range(3),values,bottom=bottom,width=.60,color=pc[p],edgecolor='white',linewidth=.35)
        bottom+=values
    for j,value in enumerate(bottom):ax.text(j,value+max(bottom)*.025,f'{int(value):,}',ha='center',va='bottom',fontsize=6.2)
    ax.set_xticks(range(3),['Baseline','Day 1','Day 5']);ax.set_ylim(0,max(bottom)*1.25)
    ax.yaxis.set_major_locator(MaxNLocator(3,integer=True));ax.set_title(tissue,fontsize=8,fontweight='bold',pad=3)
    ax.set_ylabel('Genes',fontsize=6.7,labelpad=2);ax.tick_params(length=2,pad=2)
text(.5,.494,'Bars count genes detected in at least one comparison; the four display classes do not overlap.',6.3,ha='center')

heading('C','Representative response trajectories across comparison patterns',.477)
fig.legend(handles=gh,loc='upper right',bbox_to_anchor=(.973,.490),ncol=3,fontsize=6.2,
           columnspacing=.8,handlelength=1.2,handletextpad=.4)
for ri,tissue in enumerate(T):
    for ci,p in enumerate(P):
        s=selected[(selected.tissue==tissue)&(selected.pattern==p)].iloc[0]
        x=.085+ci*.235;y=.373-ri*.117
        ax=fig.add_axes([x,y,.177,.069])
        dp=points[(points.tissue==tissue)&(points.gene_id==s.gene_id)]
        for gi,g in enumerate(G):
            means=[];sd=[];offset=(gi-1)*.13
            for day,xx in zip(['d0','d1','d5'],[0,1,5]):
                vals=dp[(dp.genotype==g)&(dp.time==day)].sort_values('rep').log2_normalized_count_plus1.to_numpy()
                ax.scatter(xx+offset+np.array([-.022,0,.022]),vals,s=7,color=gc[g],marker=gm[g],alpha=.6,lw=0,zorder=3)
                means.append(vals.mean());sd.append(vals.std(ddof=1))
            ax.errorbar(np.array([0,1,5])+offset,means,yerr=sd,color=gc[g],marker=gm[g],
                        lw=.85,ms=2.4,capsize=1.2,elinewidth=.6,zorder=4)
        ax.set_xlim(-.5,5.5);ax.set_xticks([0,1,5]);ax.set_xlabel('Day',fontsize=6.6,labelpad=1.5)
        ax.yaxis.set_major_locator(MaxNLocator(3));ax.tick_params(length=2,pad=1.5)
        title=f'{tissue} | {s.gene_id}\n{short[p]} | {day_label(s.stage)}'
        ax.set_title(title,fontsize=6.5,pad=4,color=pc[p])
text(.026,.350,'log2(normalized count + 1)',7,rotation=90,ha='center',va='center')
text(.5,.231,'Individual fish and mean ± SD; n = 3 per genotype × day. The stage in each title defines the selected comparison.',6.1,ha='center')

heading('D','Effect estimates for the selected stage in each example',.214)
dh=[Line2D([],[],color=gc[g],marker='o' if g=='HomoDel' else 's',lw=.8,ms=3,
           label='Hetero − '+g) for g in ['HomoDel','HomoIns']]
fig.legend(handles=dh,loc='upper right',bbox_to_anchor=(.97,.225),ncol=2,fontsize=6.3,
           columnspacing=1,handlelength=1.6)
for ti,tissue in enumerate(T):
    ax=fig.add_axes([.210+ti*.49,.061,.245,.116])
    labels=[]
    for row,p in enumerate(P):
        s=selected[(selected.tissue==tissue)&(selected.pattern==p)].iloc[0]
        labels.append(f'{s.gene_id}\n{day_label(s.stage)}')
        ed=effects[(effects.tissue==tissue)&(effects.gene_id==s.gene_id)]
        for g,off,marker in [('HomoDel',-.13,'o'),('HomoIns',.13,'s')]:
            e=ed[ed.comparator==g].iloc[0]
            sig=e.q_pair_family<.05
            ax.errorbar(e.effect,row+off,xerr=np.array([[e.effect-e.ci95_low],[e.ci95_high-e.effect]]),
                        fmt=marker,color=gc[g],mfc=gc[g] if sig else 'white',mec=gc[g],
                        ms=3.4,elinewidth=.8,capsize=1.6,mew=.7)
    ax.axvline(0,color='#aab1b5',lw=.6,ls='--');ax.set_yticks(range(4),labels,fontsize=6.2)
    ax.set_ylim(3.6,-.6);ax.tick_params(length=2,pad=2)
    ax.xaxis.set_major_locator(MaxNLocator(4));ax.set_xlabel('Response difference (log2 units)',fontsize=6.8,labelpad=2)
    ax.set_title(tissue,fontsize=8,fontweight='bold',pad=5)
text(.5,.027,'Filled / open: q < 0.05 / q ≥ 0.05 for each contrast. Bars: pointwise 95% Wald intervals.',6.4,ha='center')

def export(f,name):
    f.savefig(fd/(name+'.png'),dpi=300,facecolor='white')
    f.savefig(fd/(name+'.pdf'),facecolor='white',metadata={'Title':name,'Creator':'Matplotlib; archived DESeq2 contrasts'})
export(fig,'Figure3_heterozygote_patterns')
fig.savefig(fd/'Figure3_heterozygote_patterns.svg',facecolor='white')
regions={'Figure3A_complete_comparisons':(.617,.977),'Figure3B_pattern_counts':(.488,.611),
         'Figure3C_response_examples':(.225,.487),'Figure3D_effect_intervals':(0,.224)}
for name,(lo,hi) in regions.items():
    box=Bbox.from_bounds(0,lo*H,W,(hi-lo)*H)
    fig.savefig(fd/(name+'.png'),dpi=300,bbox_inches=box,facecolor='white')
    fig.savefig(fd/(name+'.pdf'),bbox_inches=box,facecolor='white')
plt.close(fig)

# Full-stage effects for the same selected genes, retaining all archived q values.
fs,axes=plt.subplots(1,2,figsize=(W,175/25.4))
fs.subplots_adjust(left=.18,right=.98,bottom=.10,top=.88,wspace=.8)
for ax,tissue in zip(axes,T):
    labels=[]
    for pi,p in enumerate(P):
        s=selected[(selected.tissue==tissue)&(selected.pattern==p)].iloc[0]
        for si,stage in enumerate(stages):
            row=pi*3+si
            labels.append(f'{s.gene_id}\n'+{'Baseline':'Baseline','Response_d1':'Day 1 response','Response_d5':'Day 5 response'}[stage])
            dd=all_effects[(all_effects.tissue==tissue)&(all_effects.gene_id==s.gene_id)&(all_effects.stage==stage)]
            for g,off,marker in [('HomoDel',-.13,'o'),('HomoIns',.13,'s')]:
                e=dd[dd.comparator==g].iloc[0]
                if not np.isfinite(e.ci95_low):continue
                ax.errorbar(e.effect,row+off,xerr=np.array([[e.effect-e.ci95_low],[e.ci95_high-e.effect]]),
                            fmt=marker,color=gc[g],mfc=gc[g] if e.q_pair_family<.05 else 'white',
                            ms=3,elinewidth=.7,capsize=1.3,mew=.65)
    for yy in [2.5,5.5,8.5]:ax.axhline(yy,color='#e3e6e8',lw=.6)
    ax.axvline(0,color='#b4bbc0',ls='--',lw=.6)
    ax.set_yticks(range(12),labels,fontsize=6.0);ax.set_ylim(11.7,-.7)
    ax.set_title(tissue,fontweight='bold',fontsize=9)
    ax.set_xlabel('Effect (log2 units)',fontsize=7);ax.xaxis.set_major_locator(MaxNLocator(4))
fs.legend(handles=dh,loc='upper center',bbox_to_anchor=(.5,.98),ncol=2,fontsize=7)
export(fs,'Figure3_supplement_all_stage_effects');plt.close(fs)
(out/'figure_QA.json').write_text(json.dumps({'main_size_mm':[210,340],'dpi':300,
    'panels':['all baseline/day1/day5 paired comparisons','four evidence-pattern counts',
              'eight representative response trajectories','two contrast estimates per example'],
    'outside_both_filter':False,'extra_fold_cutoff':None,'IUT_required':False,
    'selected_example_count':len(selected),'matplotlib':matplotlib.__version__},indent=2),encoding='utf-8')
print('Figure 3 and separate panels exported.')
