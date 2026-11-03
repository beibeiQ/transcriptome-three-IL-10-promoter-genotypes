from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/python_libs'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
W=ROOT/'work/tissue_function_revision';A=W/'analysis';I=W/'inputs'
O=ROOT/'outputs/figure_revision_v3';O.mkdir(exist_ok=True)
OLD=ROOT/'outputs/figure_revision_v2/source_tables'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':8.4,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white',
 'mathtext.fontset':'dejavusans','xtick.direction':'in','ytick.direction':'in'})
G=['HomoDel','Hetero','HomoIns'];GC={'HomoDel':'#397E85','Hetero':'#825CA6','HomoIns':'#D58A33'};GM={'HomoDel':'o','Hetero':'^','HomoIns':'s'}
NEG='#397E85';POS='#825CA6';MUTED='#BBC4C9'
CON=[('HomoDel','d1'),('HomoIns','d1'),('HomoDel','d5'),('HomoIns','d5')]
CT=['Day 1\nH−Del','Day 1\nH−Ins','Day 5\nH−Del','Day 5\nH−Ins']
cm=LinearSegmentedColormap.from_list('effect',['#397E85','#F8F8F6','#B56753'])
def rd(f):return pd.read_csv(f,dtype={'gene_id':str},low_memory=False)
ps=rd(I/'pathway_summary.csv');PIDS=ps.pathway_id.tolist()
PL=dict(zip(PIDS,['RIG-I-like receptor','Cytokine–receptor','Tight junction','Adherens junction','Apoptosis','Lysosome','Glutathione','Oxidative phosphorylation']))
scores=rd(A/'individual_scores.csv');ef=rd(A/'program_contrasts.csv');om=rd(A/'program_interactions.csv')
direct=rd(A/'direct_tissue_program_contrasts.csv');oldscore=rd(OLD/'signature_individual_scores.csv');oldcontr=rd(OLD/'signature_contrasts.csv');oldom=rd(OLD/'signature_interaction_tests.csv')
genes=rd(A/'focused_gene_individual_counts.csv');genes['log_count']=np.log2(genes.normalized_count+1)
me=rd(OLD/'ME1_existing_contrasts.csv')
def export(fig,name,main=False):
    fig.savefig(O/f'{name}.pdf');fig.savefig(O/f'{name}.svg');fig.savefig(O/f'{name}.png',dpi=300)
    if main:fig.savefig(O/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)
def head(fig,letter,title,y):
    fig.text(.025,y,letter,fontweight='bold',fontsize=13,va='top')
    fig.text(.08,y-.001,title,fontsize=9.3,va='top')
def gh():return [Line2D([],[],color=GC[g],marker=GM[g],ms=4,lw=1,label=g) for g in G]
def trajectory(ax,z,value,title,q=None,ylabel=None):
    for gi,g in enumerate(G):
        av=[];sd=[]
        for day in [0,1,5]:
            vals=z[(z.genotype==g)&(z.time==f'd{day}')][value].to_numpy();assert len(vals)==3,(g,day,title,len(vals))
            ax.scatter(np.full(3,day)+(gi-1)*.05+np.array([-.012,0,.012]),vals,s=9,marker=GM[g],c=GC[g],alpha=.45,lw=0,zorder=3)
            av.append(vals.mean());sd.append(vals.std(ddof=1))
        ax.errorbar([0,1,5],av,yerr=sd,color=GC[g],marker=GM[g],ms=3.3,lw=1,capsize=1.5,elinewidth=.6,zorder=4)
    ax.set_xticks([0,1,5]);ax.set_xlim(-.35,5.35);ax.set_xlabel('Days post-challenge',labelpad=2)
    ax.yaxis.set_major_locator(MaxNLocator(4));ax.tick_params(length=2,pad=2)
    ax.set_title(title+('\nGenotype×time: '+('*' if q<.05 else 'ns') if q is not None else ''),pad=7)
    if ylabel:ax.set_ylabel(ylabel,labelpad=3)
def point_ci(ax,estimate,lo,hi,y,q,color,marker='o',size=4):
    ax.plot([lo,hi],[y,y],color=color,lw=.85,alpha=.75)
    ax.plot([lo,hi],[y,y],ls='',marker='|',ms=3,color=color,alpha=.75)
    ax.plot(estimate,y,marker=marker,ms=size,mec=color,mfc=color if q<.05 else 'white',mew=.8,zorder=4)

# Figure 4: same fixed program screen in both organs, then the focused trajectories.
fig=plt.figure(figsize=(180/25.4,250/25.4))
head(fig,'A','Program response contrasts in intestine and spleen',.985)
extent=float(np.ceil(ef[ef.question=='Response'].estimate.abs().max()*10)/10)
for j,tis in enumerate(['Intestine','Spleen']):
    ax=fig.add_axes([.25+j*.36,.750,.28,.160]);m=np.zeros((8,4));q=np.zeros((8,4))
    for r,pid in enumerate(PIDS):
        for c,(comp,t) in enumerate(CON):
            row=ef[(ef.tissue==tis)&(ef.question=='Response')&(ef.pathway_id==pid)&(ef.comparator==comp)&(ef.time==t)].iloc[0]
            m[r,c]=row.estimate;q[r,c]=row.Q
    ax.pcolormesh(np.arange(5),np.arange(9),m,cmap=cm,vmin=-extent,vmax=extent,edgecolors='white',linewidth=.7)
    ax.set_ylim(8,0);ax.set_xticks(np.arange(4)+.5,CT);ax.set_yticks(np.arange(8)+.5,[PL[p] for p in PIDS] if j==0 else [])
    ax.set_title(tis,pad=7);ax.tick_params(length=0,pad=3)
    for r in range(8):
        for c in range(4):ax.text(c+.5,r+.5,f'{m[r,c]:+.2f}'+('*' if q[r,c]<.05 else ''),ha='center',va='center',fontsize=6.8,color='white' if abs(m[r,c])>extent*.72 else '#263238')
    for spine in ax.spines.values():spine.set_visible(False)
cax=fig.add_axes([.72,.950,.19,.008]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-extent,extent),cmap=cm),cax=cax,orientation='horizontal');cb.ax.xaxis.set_label_position('top');cb.set_label('Response difference (GSVA units)',fontsize=6.3,labelpad=3);cb.ax.tick_params(labelsize=6,length=2,pad=1)
fig.text(.25,.708,'* $Q$<0.05 across all 64 program response contrasts',fontsize=6.8)
head(fig,'B','Intestinal programs and the splenic metabolic candidate',.683)
fig.legend(handles=gh(),loc='upper center',bbox_to_anchor=(.57,.452),ncol=3,fontsize=7.3,columnspacing=1.3,handlelength=1.8)
for j,(tis,pid,title) in enumerate([('Intestine','cide04142','Intestine | Lysosome'),('Intestine','cide04210','Intestine | Apoptosis'),('Spleen','ME1','Spleen | ME1')]):
    ax=fig.add_axes([.10+j*.30,.480,.225,.140])
    if pid=='ME1':
        z=genes[(genes.tissue==tis)&(genes.gene_id=='127497534')];trajectory(ax,z,'log_count',title,ylabel='log$_2$(normalized count+1)')
    else:
        z=scores[(scores.tissue==tis)&(scores.pathway_id==pid)];q=om[(om.tissue==tis)&(om.pathway_id==pid)].iloc[0].Q
        trajectory(ax,z,'score',title,q,ylabel='GSVA score' if j==0 else None)
head(fig,'C','Splenic macrophage-associated expression programs',.416)
for j,sig in enumerate(['Macrophage','M1','M2']):
    ax=fig.add_axes([.10+j*.30,.224,.225,.130]);z=oldscore[oldscore.signature==sig]
    q=oldom[oldom.signature==sig].iloc[0].Q
    trajectory(ax,z,'score',f'{sig}-associated',q,ylabel='GSVA score' if j==0 else None)
head(fig,'D','Splenic response estimates',.181)
ax=fig.add_axes([.22,.060,.33,.078])
for yi,sig in enumerate(['Macrophage','M1','M2']):
    for j,(comp,t) in enumerate(CON):
        row=oldcontr[(oldcontr.signature==sig)&(oldcontr.question=='Response')&(oldcontr.comparator==comp)&(oldcontr.time==t)].iloc[0]
        point_ci(ax,row.estimate,row.CI_low,row.CI_high,yi+(j-1.5)*.15,row.Q,GC[comp],'o' if t=='d1' else 's',3)
ax.axvline(0,color='#89969D',lw=.6,ls='--');ax.set_yticks(range(3),['Macrophage','M1','M2']);ax.set_ylim(2.5,-.5);ax.set_xlim(-1.2,1.1);ax.set_xticks([-1,0,1]);ax.set_xlabel('Response difference (GSVA units)',fontsize=7);ax.set_title('Signature contrasts',fontsize=8,pad=6);ax.tick_params(length=2,pad=2)
ax=fig.add_axes([.73,.060,.21,.078]);focus=rd(A/'focused_gene_pairwise_effects.csv')
for yi,t in enumerate(['d1','d5']):
    for j,comp in enumerate(['HomoDel','HomoIns']):
        row=focus[(focus.gene_id=='127497534')&(focus.tissue=='Spleen')&(focus.stage=='Response_'+t)&(focus.comparator==comp)].iloc[0]
        point_ci(ax,row.effect,row.CI_low,row.CI_high,yi+(j-.5)*.17,row.Q,GC[comp],'o' if t=='d1' else 's',3.3)
ax.axvline(0,color='#89969D',lw=.6,ls='--');ax.set_yticks([0,1],['Day 1','Day 5']);ax.set_ylim(1.5,-.5);ax.set_xlim(-1.5,.5);ax.set_xticks([-1,0]);ax.set_xlabel('Response difference (log$_2$)',fontsize=7);ax.set_title('ME1 contrasts',fontsize=8,pad=6);ax.tick_params(length=2,pad=2)
fig.legend(handles=[Line2D([],[],color=GC['HomoDel'],lw=1,label='H−Del'),Line2D([],[],color=GC['HomoIns'],lw=1,label='H−Ins'),Line2D([],[],color='#515B60',marker='o',ls='',ms=3,label='Day 1'),Line2D([],[],color='#515B60',marker='s',ls='',ms=3,label='Day 5')],loc='lower center',bbox_to_anchor=(.55,.002),ncol=4,fontsize=6.6,columnspacing=1.4)
export(fig,'Figure4_intestine_spleen_functions',True)

# Figure 5: genome-wide direct tests, common program tests and linked gene examples.
gd=rd(I/'gene_direct_tissue_effects.csv.gz');fg=rd(A/'focused_gene_direct_effects.csv')
fig=plt.figure(figsize=(180/25.4,235/25.4));head(fig,'A','Genome-wide direct tissue differences',.985)
mx=float(np.nanmax(np.abs(gd[['spleen_effect','intestine_effect']].to_numpy())))*1.03
for j,(comp,t) in enumerate(CON):
    ax=fig.add_axes([.085+j*.235,.76,.18,.16]);z=gd[(gd.comparator==comp)&(gd.stage=='Response_'+t)]
    for detected in [False,True]:
        zz=z[(z.q_family<.05)==detected];col=np.where(zz.effect>0,POS,NEG) if detected else '#D4DBDE'
        ax.scatter(zz.intestine_effect,zz.spleen_effect,c=col,s=1 if not detected else 3,alpha=.55 if not detected else .8,lw=0,rasterized=True)
    ax.set_xscale('asinh',linear_width=1);ax.set_yscale('asinh',linear_width=1);ax.set_xlim(-mx,mx);ax.set_ylim(-mx,mx)
    ticks=[v for v in [-16,-4,0,4,16] if abs(v)<mx];ax.set_xticks(ticks,[str(v) for v in ticks]);ax.set_yticks(ticks,[str(v) for v in ticks])
    ax.plot([-mx,mx],[-mx,mx],ls='--',lw=.6,color='#7D888E');ax.set_title(CT[j].replace('\n',' | '),pad=6)
    ax.text(.04,.97,f'{int((z.q_family<.05).sum()):,} genes',transform=ax.transAxes,va='top',fontsize=6.8)
    ax.set_xlabel('Intestinal response',fontsize=7,labelpad=2)
    if j==0:ax.set_ylabel('Splenic response',fontsize=7,labelpad=2)
    ax.tick_params(length=2,pad=2,labelsize=6)
fig.text(.54,.717,'Hetero–homozygote response effects in log$_2$ units; signed asinh spacing',ha='center',fontsize=6.6)
head(fig,'B','Direct tissue contrasts for the eight functional programs',.684)
dx=direct[direct.question=='Response'];lim=max(abs(dx[['CI_low','CI_high']].to_numpy().ravel()))*1.03
for j,(comp,t) in enumerate(CON):
    ax=fig.add_axes([.30+j*.169,.379,.138,.245])
    for yi,pid in enumerate(PIDS):
        row=dx[(dx.pathway_id==pid)&(dx.comparator==comp)&(dx.time==t)].iloc[0]
        point_ci(ax,row.estimate,row.CI_low,row.CI_high,yi,row.Q,POS if row.estimate>0 else NEG, size=3.4)
    ax.axvline(0,color='#89969D',lw=.6,ls='--');ax.set_yticks(range(8),[PL[p] for p in PIDS] if j==0 else []);ax.set_ylim(7.6,-.6);ax.set_xlim(-lim,lim);ax.set_xticks([-1,0,1]);ax.set_title(CT[j],pad=7);ax.tick_params(length=2,pad=2)
fig.text(.63,.337,'Spleen − intestine response contrast (GSVA units)',ha='center',fontsize=7.5)
head(fig,'C','Gene-level context for lysosome, apoptosis and ME1',.299)
gids=['127512542','127509098','127497534'];gn={'127512542':'HEXB\n127512542','127509098':'PDPK1a\n127509098','127497534':'ME1\n127497534'}
lim=max(abs(fg[['ci95_low','ci95_high']].to_numpy().ravel()))*1.03
for j,(comp,t) in enumerate(CON):
    ax=fig.add_axes([.30+j*.169,.092,.138,.139])
    for yi,gid in enumerate(gids):
        row=fg[(fg.gene_id==gid)&(fg.comparator==comp)&(fg.stage=='Response_'+t)].iloc[0]
        point_ci(ax,row.effect,row.ci95_low,row.ci95_high,yi,row.q_family,POS if row.effect>0 else NEG,size=3.4)
    ax.axvline(0,color='#89969D',lw=.6,ls='--');ax.set_yticks(range(3),[gn[g] for g in gids] if j==0 else []);ax.set_ylim(2.6,-.6);ax.set_xlim(-lim,lim);ax.set_xticks([-4,0,4]);ax.set_title(CT[j].replace('\n',' | '),fontsize=7,pad=7);ax.tick_params(length=2,pad=2)
fig.text(.63,.052,'Spleen − intestine response contrast (log$_2$ units)',ha='center',fontsize=7.5)
fig.legend(handles=[Line2D([],[],color=POS,marker='o',ms=3,lw=0,label='Positive direct effect'),Line2D([],[],color=NEG,marker='o',ms=3,lw=0,label='Negative direct effect'),Line2D([],[],color='#515B60',marker='o',mfc='white',ms=3,lw=0,label='$Q$≥0.05')],loc='lower center',bbox_to_anchor=(.55,.003),ncol=3,fontsize=6.5)
export(fig,'Figure5_functional_tissue_integration',True)

# Supplement: all program trajectories in both organs.
fig=plt.figure(figsize=(180/25.4,250/25.4))
fig.legend(handles=gh(),loc='upper center',bbox_to_anchor=(.53,.996),ncol=3,fontsize=8)
for pi,pid in enumerate(PIDS):
    row=pi//2; pair=pi%2
    for ti,tis in enumerate(['Intestine','Spleen']):
        col=pair*2+ti;ax=fig.add_axes([.08+col*.24,.785-row*.234,.175,.132])
        z=scores[(scores.pathway_id==pid)&(scores.tissue==tis)];q=om[(om.pathway_id==pid)&(om.tissue==tis)].iloc[0].Q
        trajectory(ax,z,'score',f'{PL[pid]}\n{tis}',ylabel='GSVA score' if col==0 else None)
        ax.set_title(f'{PL[pid]}\n{tis} | '+('*' if q<.05 else 'ns'),fontsize=6.8,pad=6);ax.tick_params(labelsize=6);ax.set_xlabel('Day',fontsize=6.8)
        lo=scores[scores.pathway_id==pid].score.min();hi=scores[scores.pathway_id==pid].score.max();pad=.15*(hi-lo)
        ax.set_ylim(lo-pad,hi+pad)
export(fig,'Supplementary_Figure_S5_all_program_trajectories')

# Supplement: gene trajectories and predefined scoring sensitivity.
fig=plt.figure(figsize=(180/25.4,230/25.4));head(fig,'A','Expression context of the three linked genes',.985)
fig.legend(handles=gh(),loc='upper center',bbox_to_anchor=(.56,.96),ncol=3,fontsize=7.5)
for ti,tis in enumerate(['Intestine','Spleen']):
    for gi,gid in enumerate(gids):
        ax=fig.add_axes([.10+gi*.30,.73-ti*.235,.225,.13]);z=genes[(genes.tissue==tis)&(genes.gene_id==gid)]
        trajectory(ax,z,'log_count',f'{tis} | '+gn[gid].replace('\n',' '),ylabel='log$_2$(normalized count+1)' if gi==0 else None)
head(fig,'B','Score-definition sensitivity',.425)
sens=rd(A/'scoring_sensitivity_comparison.csv');sene=rd(A/'within_tissue_scoring_sensitivity.csv')
ax=fig.add_axes([.09,.105,.29,.245])
ax.scatter(sens.pooled_effect,sens.separate_effect,c=np.where(sens.pooled_Q<.05,POS,'#AFBBC1'),s=16,alpha=.85,edgecolors='none')
lim=max(np.abs(sens[['pooled_effect','separate_effect']].to_numpy()).ravel())*1.15
ax.plot([-lim,lim],[-lim,lim],color='#79868E',lw=.7,ls='--');ax.axhline(0,color='#BBC4C9',lw=.5);ax.axvline(0,color='#BBC4C9',lw=.5);ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
ax.set_xlabel('Pooled scoring effect');ax.set_ylabel('Separate-tissue scoring effect');ax.set_title('All 64 response contrasts\n60 share effect direction',fontsize=8)
ax=fig.add_axes([.69,.105,.24,.245])
chosen=[('Intestine','cide04142','HomoDel','d1'),('Intestine','cide04210','HomoIns','d5'),('Spleen','cide04520','HomoIns','d5')]
labels=[]
for yi,(tis,pid,comp,t) in enumerate(chosen):
    labels.append(f'{tis}: {PL[pid]}\nH−'+('Del' if comp=='HomoDel' else 'Ins')+f' | day {t[1]}')
    for si,(table,color) in enumerate([(ef,POS),(sene,NEG)]):
        row=table[(table.tissue==tis)&(table.pathway_id==pid)&(table.comparator==comp)&(table.time==t)&(table.question=='Response')].iloc[0]
        point_ci(ax,row.estimate,row.CI_low,row.CI_high,yi+(si-.5)*.22,row.Q,color,size=3.8)
        ax.text(.98,yi+(si-.5)*.22,'*' if row.Q<.05 else 'ns',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6.5,color=color)
ax.axvline(0,color='#89969D',lw=.6,ls='--');ax.set_ylim(2.5,-.5);ax.set_yticks(range(3),labels,fontsize=6.8);ax.set_xlim(-.7,1.65);ax.set_xticks([-.5,0,.5,1]);ax.set_xlabel('Response difference (GSVA units)');ax.set_title('Selected contrasts:\nscore-definition sensitivity',fontsize=7.7)
fig.legend(handles=[Line2D([],[],color=POS,marker='o',lw=1,ms=3,label='Pooled across tissues'),Line2D([],[],color=NEG,marker='o',mfc='white',lw=1,ms=3,label='Within each tissue')],loc='lower center',bbox_to_anchor=(.58,.015),ncol=2,fontsize=7)
export(fig,'Supplementary_Figure_S6_gene_context_sensitivity')

# Diagnostics remain separate from scientific main figures.
fig,axs=plt.subplots(4,4,figsize=(10,9));diag=rd(A/'score_model_diagnostics.csv')
for ax,(tis,pid) in zip(axs.ravel(),[(t,p) for p in PIDS for t in ['Intestine','Spleen']]):
    z=diag[(diag.tissue==tis)&(diag.pathway_id==pid)]
    ax.scatter(z.fitted,z.residual,s=9,c=[GC[g] for g in z.genotype]);ax.axhline(0,color='#666',ls='--',lw=.6);ax.set_title(f'{tis} | {PL[pid]}',fontsize=7);ax.set_xlabel('Fitted score',fontsize=7);ax.set_ylabel('Residual',fontsize=7)
fig.tight_layout();export(fig,'Analysis_diagnostics_score_residuals')
print('Rendered Figures4/5, supplementaryS5/S6 and diagnostics.')
