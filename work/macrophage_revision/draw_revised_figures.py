# Project paths.
import os as _il10_os
import sys as _il10_sys
from pathlib import Path as _Il10Path
_il10_repo = next((p for p in _Il10Path(__file__).resolve().parents if (p / ".il10-project").is_file()), None)
if _il10_repo is None: raise RuntimeError("Cannot locate the analysis checkout")
_il10_sys.path.insert(0, str(_il10_repo))
from config.paths import REPO_ROOT, data_path, result_path, annotation_path

from pathlib import Path
import sys,json,shutil,re
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python_libs'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.text import Text

R=Path(__file__).resolve().parents[2];W=R/'work/macrophage_revision';O=R/'outputs/figure_revision_v2';O.mkdir(exist_ok=True)
S=O/'source_tables';S.mkdir(exist_ok=True)
P=data_path("figure_source_tables")
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
G=['HomoDel','Hetero','HomoIns'];GC={'HomoDel':'#397E85','Hetero':'#825CA6','HomoIns':'#D58A33'};GM={'HomoDel':'o','Hetero':'^','HomoIns':'s'}
PAT=['Both_same','Both_opposite','Del_only','Ins_only'];PC={'Both_same':'#7A659E','Both_opposite':'#C96D92','Del_only':'#287D8E','Ins_only':'#D79533'}
PL={'Both_same':'Both, same signs','Both_opposite':'Both, opposite signs','Del_only':'H−Del only','Ins_only':'H−Ins only'}
ST=['Baseline','Response_d1','Response_d5'];SL={'Baseline':'Baseline','Response_d1':'Day 1 response','Response_d5':'Day 5 response'}
T=['Spleen','Intestine']
def rd(path):return pd.read_csv(path,dtype={'gene_id':str},low_memory=False)
def export(fig,name,main=False):
    # Figure labels.
    for txt in fig.texts:
        if txt.get_position()[1] >= .985 or txt.get_text().startswith(('Log2 effects;', 'Joint candidates across', 'A, C:', 'A: individual', 'C: Hetero', 'All 40 distinct')):
            txt.set_visible(False)
    for txt in fig.findobj(Text):
        value=txt.get_text()
        if value not in list('ABCD'): txt.set_fontweight('normal')
        value=re.sub(r'(?<=\d),(?=\d{3}(?:\D|$))','\u202f',value)
        value=value.replace('normalized count + 1','normalized count+1')
        value=value.replace('Log2 effects','log$_2$ effects')
        value=re.sub(r'\b([PQn])\s*([=<>≥])\s*',r'$\1$\2',value)
        value=value.replace('joint Q','joint $Q$').replace(' P (',' $P$ (')
        txt.set_text(value)
    for ext in ['pdf','svg','png']:
        fig.savefig(O/f'{name}.{ext}',dpi=600 if main else 300)
    if main:fig.savefig(O/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)
def heading(fig,letter,title,y):
    fig.text(.025,y,letter,fontsize=12,fontweight='bold',va='top')
    fig.text(.08,y-.002,title,fontsize=8.8,fontweight='bold',va='top')
def gh():return [Line2D([],[],color=GC[g],marker=GM[g],ms=3,lw=1,label=g) for g in G]
def ch():return [Line2D([],[],color=GC[g],marker=GM[g],ms=3,lw=1,label='Hetero − '+g) for g in ['HomoDel','HomoIns']]
def trajectory(ax,z,value,title=None,ylabel=None):
    for gi,g in enumerate(G):
        zz=z[z.genotype==g];av=[];sd=[]
        for day in [0,1,5]:
            vals=zz[zz.time==f'd{day}'][value].to_numpy();assert len(vals)==3,(g,day,len(vals))
            ax.scatter(np.full(3,day)+(gi-1)*.065+np.array([-.012,0,.012]),vals,s=7,color=GC[g],marker=GM[g],alpha=.45,lw=0,zorder=3)
            av.append(vals.mean());sd.append(vals.std(ddof=1))
        ax.errorbar([0,1,5],av,yerr=sd,color=GC[g],marker=GM[g],ms=3,lw=.85,capsize=1.4,elinewidth=.6,label=g,zorder=4)
    ax.set_xticks([0,1,5]);ax.set_xlim(-.4,5.4);ax.set_xlabel('Days post-challenge',labelpad=2)
    ax.yaxis.set_major_locator(MaxNLocator(4));ax.tick_params(length=2,pad=2)
    if title:ax.set_title(title,pad=6)
    if ylabel:ax.set_ylabel(ylabel,labelpad=3)

# New Figure 3 combines former 3A/B and 4A without changing tests or classifications.
w=rd(P/'Fig3/source_tables/A_all_genes_paired_comparisons.csv.gz')
counts=rd(P/'Fig3/source_tables/B_pattern_counts.csv')
joint=rd(P/'Fig4/source_tables/A_all_gene_test_points.csv.gz')
fig=plt.figure(figsize=(180/25.4,235/25.4))
fig.text(.5,.995,'Heterozygote-associated patterns and joint candidate evidence',ha='center',va='top',fontsize=10,fontweight='bold')
heading(fig,'A','Direction and timing of differences from both homozygotes',.966)
fig.legend(handles=[Line2D([],[],color=PC[p],marker='D' if p=='Both_opposite' else 'o',ls='',ms=3,label=PL[p]) for p in PAT],loc='upper center',bbox_to_anchor=(.55,.938),ncol=2,fontsize=6.3,columnspacing=1.8,handletextpad=.35)
for ri,t in enumerate(T):
 for ci,st in enumerate(ST):
    ax=fig.add_axes([.10+ci*.303,.695-ri*.195,.235,.155]);d=w[(w.tissue==t)&(w.stage==st)]
    finite=np.isfinite(d.effect_D)&np.isfinite(d.effect_I)
    for pat in ['Neither_detected','Incomplete']+PAT:
        z=d[finite&(d.pattern==pat)];color=PC.get(pat,'#CED5D9' if pat=='Neither_detected' else '#8B9399')
        ax.scatter(z.effect_D,z.effect_I,s=1.1 if pat=='Neither_detected' else (7 if pat=='Both_opposite' else 3),c=color,marker='D' if pat=='Both_opposite' else ('x' if pat=='Incomplete' else 'o'),alpha=.5 if pat=='Neither_detected' else .8,linewidths=.35 if pat=='Incomplete' else 0,rasterized=True)
    ax.set_xscale('asinh',linear_width=1);ax.set_yscale('asinh',linear_width=1)
    extent=max(2,float(np.nanmax(abs(d.loc[finite,['effect_D','effect_I']].to_numpy()))))*1.04
    ax.set_xlim(-extent,extent);ax.set_ylim(-extent,extent)
    ticks=[v for v in [-16,-4,-1,0,1,4,16] if abs(v)<extent];ax.set_xticks(ticks,ticks);ax.set_yticks(ticks,ticks)
    ax.axhline(0,color='#ADB5B9',lw=.5,ls='--');ax.axvline(0,color='#ADB5B9',lw=.5,ls='--')
    n=int((d.pattern.isin(PAT)&d.fully_testable).sum())
    ax.text(.03,.96,f'{n:,} detected in ≥1 contrast',transform=ax.transAxes,va='top',fontsize=5.7,bbox={'facecolor':'white','edgecolor':'none','alpha':.83,'pad':1})
    if ri==0:ax.set_title(SL[st],fontsize=7.5,fontweight='bold')
    if ci==0:ax.set_ylabel(t+'\nHetero − HomoIns')
    if ri==1:ax.set_xlabel('Hetero − HomoDel',labelpad=2)
    ax.tick_params(length=2,pad=1.5)
fig.text(.54,.456,'Log2 effects; signed asinh spacing. Gray: neither detected; ×: incomplete tests.',ha='center',fontsize=6)
heading(fig,'B','Distribution of detected comparison patterns',.433)
for ti,t in enumerate(T):
    ax=fig.add_axes([.11+ti*.49,.318,.345,.075]);bottom=np.zeros(3)
    for pat in PAT:
        v=np.array([counts.loc[(counts.tissue==t)&(counts.stage==st)&(counts.pattern==pat),'genes'].sum() for st in ST])
        ax.bar(range(3),v,bottom=bottom,color=PC[pat],width=.62,edgecolor='white',linewidth=.3);bottom+=v
    for x,v in enumerate(bottom):ax.text(x,v+max(bottom)*.025,f'{int(v):,}',ha='center',fontsize=6)
    ax.set_xticks(range(3),['Baseline','Day 1','Day 5']);ax.set_ylim(0,max(bottom)*1.23);ax.set_ylabel('Genes');ax.set_title(t,fontweight='bold',pad=5);ax.yaxis.set_major_locator(MaxNLocator(3))
heading(fig,'C','Joint evidence for differences from both homozygotes',.271)
for ci,(t,fam) in enumerate([(t,f) for t in T for f in ['Baseline','Response']]):
    ax=fig.add_axes([.085+ci*.238,.099,.195,.128]);z=joint[(joint.tissue==t)&(joint.family==fam)&joint.p_D.notna()&joint.p_I.notna()].copy()
    z['xx']=-np.log10(z.p_D.clip(lower=1e-300));z['yy']=-np.log10(z.p_I.clip(lower=1e-300))
    priority=z['finite_effect_priority_pass'].fillna(False).astype(bool) if 'finite_effect_priority_pass' in z else pd.Series(False,index=z.index)
    core=z.q_joint<.05
    if 'finite_effect_priority_pass' not in z:
        cand=rd(P/'Fig4/source_tables/B_complete_candidate_order.csv')
        pr=set(cand.loc[(cand.tissue==t)&(cand.family==fam)&cand.finite_effect_priority_pass,'gene_id']);priority=z.gene_id.isin(pr)
    for mask,color,size in [(~core,'#C8CFD3',2),(core&~priority,'#D79A3D',5),(core&priority,'#825CA6',6)]:
        ax.scatter(z.loc[mask,'xx'].clip(upper=25),z.loc[mask,'yy'].clip(upper=25),s=size,c=color,alpha=.6 if color=='#C8CFD3' else .9,lw=0,rasterized=True)
    if core.any():
        # BH was performed across tissues within each family: use that same boundary.
        cutoff=-np.log10(joint.loc[(joint.family==fam)&(joint.q_joint<.05),'p_joint'].max());ax.axvline(cutoff,color='#7E8894',ls='--',lw=.6);ax.axhline(cutoff,color='#7E8894',ls='--',lw=.6)
    ax.set_xlim(-.6,25.6);ax.set_ylim(-.6,25.6);ax.set_xticks([0,10,20]);ax.set_yticks([0,10,20]);ax.set_title(t+' | '+fam,fontsize=6.8,pad=6)
    ax.text(.03,.94,f'n = {len(z):,}',transform=ax.transAxes,va='top',fontsize=5.8,color='#6F7A82')
    ax.text(.98,.04,f'{core.sum()} joint / {int((core&priority).sum())} priority',transform=ax.transAxes,ha='right',fontsize=5.8)
    ax.set_xlabel('−log$_{10}$ P (H−Ins)' if False else '−log$_{10}$ P (H−Del)',fontsize=6.1,labelpad=2)
    if ci==0:ax.set_ylabel('−log$_{10}$ P (H−Ins)',fontsize=6.1,labelpad=2)
fig.legend(handles=[Line2D([],[],ls='',marker='o',color=c,ms=3,label=l) for c,l in [('#C8CFD3','Other tested genes'),('#D79A3D','Other joint candidates'),('#825CA6','Passed finite-effect checks')]],loc='lower center',bbox_to_anchor=(.53,.039),ncol=3,fontsize=5.9,columnspacing=1.2)
fig.text(.5,.019,'Joint candidates across baseline and response: 24 distinct splenic and 247 intestinal genes.',ha='center',fontsize=6.3)
export(fig,'Figure3_merged',True)

# New Figure 4: transparent signature trajectories, response contrasts and ME1.
scores=rd(W/'signature_individual_scores.csv');tests=rd(W/'signature_contrasts.csv');omni=rd(W/'signature_interaction_tests.csv')
me=rd(W/'ME1_individual_counts.csv');me=me[me.tissue=='Spleen'].copy();me['log_count']=np.log2(me.normalized_count+1)
meeff=rd(W/'ME1_existing_contrasts.csv');meeff=meeff[meeff.tissue=='Spleen']
fig=plt.figure(figsize=(180/25.4,205/25.4))
fig.text(.5,.995,'Splenic macrophage-associated signatures and ME1 dynamics',ha='center',va='top',fontsize=10,fontweight='bold')
heading(fig,'A','Macrophage-associated expression programs',.964)
fig.legend(handles=gh(),loc='upper center',bbox_to_anchor=(.57,.944),ncol=3,fontsize=6.5)
siglabels={'Macrophage':'Macrophage-associated','M1':'M1-associated','M2':'M2-associated'}
for i,sig in enumerate(['Macrophage','M1','M2']):
    ax=fig.add_axes([.105+i*.303,.713,.225,.166]);z=scores[scores.signature==sig];q=float(omni.loc[omni.signature==sig,'Q'].iloc[0]);size={'Macrophage':9,'M1':17,'M2':19}[sig]
    trajectory(ax,z,'score',siglabels[sig]+f' ({size} genes)','GSVA score' if i==0 else None)
    ax.text(.5,-.36,f'Interaction Q = {q:.3g}',ha='center',transform=ax.transAxes,fontsize=6.2)
heading(fig,'B','Differences in the response from baseline',.627)
fig.legend(handles=ch(),loc='upper right',bbox_to_anchor=(.98,.621),ncol=2,fontsize=6.2)
for ti,time in enumerate(['d1','d5']):
    ax=fig.add_axes([.225+ti*.46,.423,.275,.147])
    for si,sig in enumerate(['Macrophage','M1','M2']):
      for comp,off in [('HomoDel',-.14),('HomoIns',.14)]:
        e=tests[(tests.signature==sig)&(tests.question=='Response')&(tests.time==time)&(tests.comparator==comp)].iloc[0];hit=e.Q<.05
        ax.errorbar(e.estimate,si+off,xerr=[[e.estimate-e.CI_low],[e.CI_high-e.estimate]],fmt=GM[comp],ms=3.5,color=GC[comp],mfc=GC[comp] if hit else 'white',elinewidth=.75,capsize=1.5)
        if hit:ax.text(e.CI_low,si+off-.2,f'Q={e.Q:.3f}',fontsize=5.8,color=GC[comp])
    ax.axvline(0,color='#A7B0B6',ls='--',lw=.6);ax.set_yticks(range(3),['Macrophage','M1-associated','M2-associated']);ax.set_ylim(2.6,-.7);ax.set_xlim(-1.2,1.2);ax.set_xticks([-1,0,1]);ax.set_xlabel('Response difference in GSVA score');ax.set_title('Day '+time[1],fontweight='bold')
heading(fig,'C','ME1 expression trajectory',.354)
heading(fig,'D','ME1 response contrasts',.354)
# Reposition the D heading into the second half.
fig.texts[-2].set_x(.54);fig.texts[-1].set_x(.60)
ax=fig.add_axes([.11,.113,.33,.185]);trajectory(ax,me,'log_count','127497534 | malic enzyme 1','log$_2$(normalized count + 1)')
ax=fig.add_axes([.65,.113,.29,.185])
for i,stage in enumerate(['Response_d1','Response_d5']):
 for comp,off in [('HomoDel',-.14),('HomoIns',.14)]:
    e=meeff[(meeff.stage==stage)&(meeff.comparator==comp)].iloc[0];hit=e.q_pair_family<.05
    ax.errorbar(e.effect,i+off,xerr=[[e.effect-e.ci95_low],[e.ci95_high-e.effect]],fmt=GM[comp],ms=3.5,color=GC[comp],mfc=GC[comp] if hit else 'white',elinewidth=.75,capsize=1.5)
    if hit:ax.text(e.ci95_low,i+off-.18,f'Q={e.q_pair_family:.1e}',fontsize=5.7,color=GC[comp])
ax.axvline(0,color='#A7B0B6',ls='--',lw=.6);ax.set_yticks([0,1],['Day 1','Day 5']);ax.set_ylim(1.5,-.6);ax.set_xlabel('Response difference (log$_2$ units)');ax.xaxis.set_major_locator(MaxNLocator(4))
fig.text(.5,.044,'A, C: individual fish and mean ± SD; n = 3 fish per genotype and day.\nB, D: pointwise 95% intervals; filled/open symbols indicate Q < 0.05 / Q ≥ 0.05 in the stated test family.',ha='center',fontsize=6.1)
export(fig,'Figure4_macrophage_ME1',True)

# Supplement S2: old pattern examples and their effect intervals.
sel=rd(P/'Fig3/source_tables/C_selected_response_examples.csv');points=rd(P/'Fig3/source_tables/C_individual_expression_points.csv');eff=rd(P/'Fig3/source_tables/D_selected_effects_and_intervals.csv')
fig=plt.figure(figsize=(180/25.4,190/25.4));fig.text(.5,.99,'Supplementary Figure S2. Examples of comparison patterns',ha='center',va='top',fontweight='bold',fontsize=10)
heading(fig,'A','Representative trajectories for the four detection patterns',.954)
fig.legend(handles=gh(),loc='upper center',bbox_to_anchor=(.58,.933),ncol=3,fontsize=6.4)
for ri,t in enumerate(T):
 for ci,pat in enumerate(PAT):
    e=sel[(sel.tissue==t)&(sel.pattern==pat)].iloc[0];ax=fig.add_axes([.09+ci*.24,.686-ri*.247,.18,.145])
    trajectory(ax,points[(points.tissue==t)&(points.gene_id==e.gene_id)],'log2_normalized_count_plus1',f'{t} | {e.gene_id}\n{PL[pat]} | day {e.stage[-1]}')
    ax.title.set_color(PC[pat]);ax.title.set_fontsize(6.6)
    if ci==0:ax.set_ylabel('log$_2$(normalized count + 1)',fontsize=6)
heading(fig,'B','The two response contrasts at each selected stage',.358)
fig.legend(handles=ch(),loc='upper right',bbox_to_anchor=(.97,.354),ncol=2,fontsize=6)
for ti,t in enumerate(T):
 ax=fig.add_axes([.23+ti*.47,.09,.25,.2]);labs=[]
 for i,pat in enumerate(PAT):
    se=sel[(sel.tissue==t)&(sel.pattern==pat)].iloc[0];labs.append(f'{se.gene_id}\nday {se.stage[-1]}')
    for comp,off in [('HomoDel',-.14),('HomoIns',.14)]:
      e=eff[(eff.tissue==t)&(eff.gene_id==se.gene_id)&(eff.comparator==comp)].iloc[0]
      ax.errorbar(e.effect,i+off,xerr=[[e.effect-e.ci95_low],[e.ci95_high-e.effect]],fmt=GM[comp],ms=3,color=GC[comp],mfc=GC[comp] if e.q_pair_family<.05 else 'white',elinewidth=.7,capsize=1.4)
 ax.set_yticks(range(4),labs,fontsize=6);ax.set_ylim(3.5,-.5);ax.axvline(0,color='#A7B0B6',ls='--',lw=.6);ax.set_title(t,fontweight='bold');ax.set_xlabel('Response difference (log$_2$ units)')
fig.text(.5,.024,'A: individual fish and mean ± SD. B: pointwise 95% Wald intervals; filled/open = Q < 0.05 / Q ≥ 0.05.',ha='center',fontsize=6)
export(fig,'Supplementary_Figure_S2_pattern_examples')

# Supplement S3 preserves the candidate ranking and all ten original examples.
cand=rd(P/'Fig4/source_tables/B_complete_candidate_order.csv');cs=rd(P/'Fig4/source_tables/C_display_genes_fixed_rule.csv');cp=rd(P/'Fig4/source_tables/C_individual_fish_expression.csv');ce=rd(P/'Fig4/source_tables/D_effects_Wald_intervals_and_deletion_ranges.csv')
fig=plt.figure(figsize=(180/25.4,255/25.4));fig.text(.5,.995,'Supplementary Figure S3. Joint candidates and sensitivity',ha='center',va='top',fontsize=10,fontweight='bold')
heading(fig,'A','Ordering of all jointly detected candidate entries',.964)
for ti,t in enumerate(T):
    ax=fig.add_axes([.10+ti*.49,.79,.35,.12]);z=cand[cand.tissue==t].sort_values('entry_order').reset_index(drop=True)
    for _,e in z.iterrows():
        color='#397E85' if e.family=='Baseline' else '#D58A33';marker='o' if e.family=='Baseline' else 's'
        ax.scatter(e.entry_order,-np.log10(e.q_joint),s=13,marker=marker,edgecolor=color,facecolor=color if e.finite_effect_priority_pass else 'white',linewidth=.7)
    ax.axhline(-np.log10(.05),color='#A7B0B6',ls='--',lw=.5);ax.set_title(f'{t} | {len(z)} entries; {z.gene_id.nunique()} genes',fontsize=7);ax.set_xlabel('Order by joint Q within tissue');ax.set_ylabel('−log$_{10}$ joint Q');ax.tick_params(length=2)
fig.legend(handles=[Line2D([],[],ls='',marker='o',color='#397E85',ms=3,label='Baseline'),Line2D([],[],ls='',marker='s',color='#D58A33',ms=3,label='Response'),Line2D([],[],ls='',marker='o',color='#809099',mfc='white',ms=3,label='Other joint candidates')],loc='upper center',bbox_to_anchor=(.55,.754),ncol=3,fontsize=6)
heading(fig,'B','Highest-ranked distinct candidates passing effect checks',.719)
fig.legend(handles=gh(),loc='upper right',bbox_to_anchor=(.99,.709),ncol=3,fontsize=6)
for ri,t in enumerate(T):
 for ci,(_,e) in enumerate(cs[cs.tissue==t].sort_values('display_index').iterrows()):
    ax=fig.add_axes([.075+ci*.187,.514-ri*.191,.145,.12]);z=cp[(cp.tissue==t)&(cp.gene_id==e.gene_id)]
    short_labels={'127510871':'plexin-C1-like','127510740':'FGD4b-like','127511380':'DEF6-like','127498082':'NLRP12-like','127497534':'malic enzyme 1','127502293':'fos homolog Ab','127505992':'Uncharacterized','127497599':'Uncharacterized','127494715':'Uncharacterized homolog','127502635':'Jun dimerization\nprotein 2b'}
    label=short_labels[e.gene_id]
    trajectory(ax,z,'log2_normalized_count_plus1',f'{t} | [{int(e.display_index)}] {e.gene_id}\n{label}')
    ax.title.set_fontsize(6);ax.set_xlabel('Day',fontsize=6);ax.tick_params(labelsize=5.5)
    if ci==0:ax.set_ylabel('log$_2$(normalized count + 1)',fontsize=6)
heading(fig,'C','Effect estimates and leave-one-fish ranges',.265)
for ti,t in enumerate(T):
 genes=cs[cs.tissue==t].sort_values('display_index').gene_id.tolist()
 for si,stage in enumerate(ST):
    ax=fig.add_axes([.14+ti*.48+si*.13,.07,.105,.135])
    for gi,gid in enumerate(genes):
     for comp,off in [('HomoDel',-.13),('HomoIns',.13)]:
      e=ce[(ce.tissue==t)&(ce.gene_id==gid)&(ce.stage==stage)&(ce.comparator==comp)].iloc[0]
      ax.plot([e.deletion_min,e.deletion_max],[gi+off]*2,color='#BDC4C8',lw=3,solid_capstyle='butt')
      ax.errorbar(e.effect,gi+off,xerr=[[e.effect-e.ci95_low],[e.ci95_high-e.effect]],fmt=GM[comp],ms=2.5,color=GC[comp],elinewidth=.6,capsize=1)
    ax.set_yticks(range(5),[f'[{i+1}] {gid}' for i,gid in enumerate(genes)] if si==0 else ['']*5,fontsize=5.2);ax.set_ylim(4.6,-.6);ax.axvline(0,color='#A7B0B6',ls='--',lw=.5);ax.set_title(('Baseline' if stage=='Baseline' else 'Day '+stage[-1]+' response'),fontsize=6)
    ax.xaxis.set_major_locator(MaxNLocator(3));ax.tick_params(length=2,labelsize=5.5)
 fig.text(.35+ti*.48,.227,t,ha='center',fontweight='bold',fontsize=7)
fig.text(.5,.027,'C: Hetero−HomoDel (teal) and Hetero−HomoIns (orange); pointwise 95% Wald intervals.\nGray segments: ranges over 27 fixed-parameter leave-one-fish estimates. Baseline and response use separate BH families.',ha='center',fontsize=6)
export(fig,'Supplementary_Figure_S3_joint_candidates')

# Supplement S4: sensitivity to overlapping markers and all individual-gene profiles.
members=rd(W/'marker_membership_and_annotation.csv');meta=rd(W/'sample_metadata.csv');vst=rd(P/'Fig2/inputs/Spleen_VST.csv.gz').set_index('gene_id')
fig=plt.figure(figsize=(180/25.4,255/25.4));fig.text(.5,.995,'Supplementary Figure S4. Marker membership and sensitivity',ha='center',va='top',fontsize=10,fontweight='bold')
heading(fig,'A','Scores after removing the five shared M1/M2 markers',.962)
fig.legend(handles=gh(),loc='upper right',bbox_to_anchor=(.98,.941),ncol=3,fontsize=6.3)
for i,sig in enumerate(['M1_without_shared','M2_without_shared']):
 ax=fig.add_axes([.11+i*.49,.76,.34,.13]);trajectory(ax,scores[scores.signature==sig],'score',sig.replace('_',' '),'GSVA score' if i==0 else None)
heading(fig,'B','Complete mapped marker profiles across genotype and day',.706)
order=[];labels=[];membership=[]
for label in ['Macrophage','M1','M2']:
 for _,e in members[members.signature==label].sort_values('gene_symbol').iterrows():
    if e.gene_id in order:continue
    order.append(e.gene_id);labels.append(e.gene_symbol.replace('LOC127518037','properdin-like').replace('LOC127498996','npc2-like'));membership.append('Shared M1/M2' if e.in_shared_M1_M2 else label)
cols=[];clabels=[]
for g in G:
 for t in ['d0','d1','d5']:
    cols.append(meta.loc[(meta.tissue=='Spleen')&(meta.genotype==g)&(meta.time==t),'sample'].tolist());clabels.append(g+'\n'+t[1])
spmeta=meta[meta.tissue=='Spleen'];vv=vst.loc[order,spmeta['sample']].copy();z=(vv.sub(vv.mean(axis=1),axis=0)).div(vv.std(axis=1,ddof=1),axis=0)
heat=np.array([z[c].mean(axis=1).to_numpy() for c in cols]).T
ax=fig.add_axes([.29,.105,.49,.54]);im=ax.imshow(heat,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2,interpolation='nearest');ax.set_yticks(range(len(order)),labels,fontsize=5.6);ax.set_xticks(range(9),clabels,fontsize=5.6);ax.tick_params(length=0,pad=3)
for v in [2.5,5.5]:ax.axvline(v,color='white',lw=1)
for i in range(1,len(order)):
 if membership[i]!=membership[i-1]:ax.axhline(i-.5,color='white',lw=.7)
for i,lab in enumerate(membership):ax.text(9.0,i,lab,ha='left',va='center',fontsize=5.1)
cax=fig.add_axes([.88,.66,.07,.008]);fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2]);cax.set_title('Gene-wise z',fontsize=6,pad=3)
fig.text(.5,.046,'All 40 distinct mapped markers are shown. Values are group means of gene-wise standardized spleen VST expression.\nShared M1/M2 membership is retained in the primary analysis and removed only for the sensitivity analysis.',ha='center',fontsize=6)
export(fig,'Supplementary_Figure_S4_marker_sensitivity')

# Explicit source-table copies keep figures and their evidence together.
copies={
 'Figure3_all_paired_comparisons.csv.gz':P/'Fig3/source_tables/A_all_genes_paired_comparisons.csv.gz',
 'Figure3_pattern_counts.csv':P/'Fig3/source_tables/B_pattern_counts.csv',
 'Figure3_joint_test_points.csv.gz':P/'Fig4/source_tables/A_all_gene_test_points.csv.gz',
 'Figure3_joint_candidate_order.csv':P/'Fig4/source_tables/B_complete_candidate_order.csv',
 'S2_selected_examples.csv':P/'Fig3/source_tables/C_selected_response_examples.csv',
 'S2_individual_points.csv':P/'Fig3/source_tables/C_individual_expression_points.csv',
 'S2_effects.csv':P/'Fig3/source_tables/D_selected_effects_and_intervals.csv',
 'S3_selected_candidates.csv':P/'Fig4/source_tables/C_display_genes_fixed_rule.csv',
 'S3_individual_points.csv':P/'Fig4/source_tables/C_individual_fish_expression.csv',
 'S3_effects_and_sensitivity.csv':P/'Fig4/source_tables/D_effects_Wald_intervals_and_deletion_ranges.csv'}
for dest,src in copies.items():shutil.copy2(src,S/dest)
for name in ['signature_individual_scores.csv','signature_group_summary.csv','signature_contrasts.csv','signature_interaction_tests.csv','ME1_existing_contrasts.csv','ME1_individual_counts.csv','marker_membership_and_annotation.csv','ME1_adjusted_score_correlations.csv','archived_score_reproduction.csv','marker_provenance.json','R_sessionInfo.txt','marker_sets.gmt']:
    shutil.copy2(W/name,S/name)
pd.DataFrame(heat,index=order,columns=[g+'_'+t for g in G for t in ['d0','d1','d5']]).rename_axis('gene_id').to_csv(S/'S4_marker_heatmap_group_means.csv')
shutil.copy2(__file__,O/'draw_revised_figures.py');shutil.copy2(W/'analyse_signatures.R',O/'analyse_signatures.R')
print('Generated figures:',O)
