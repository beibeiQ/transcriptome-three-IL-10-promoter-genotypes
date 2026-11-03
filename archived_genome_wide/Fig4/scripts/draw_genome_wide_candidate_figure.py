"""Figure 4: all-gene discovery, complete evidence ordering, and declared display subset."""
from pathlib import Path
import sys,json,shutil,hashlib,textwrap
sys.path.insert(0,str(Path(__file__).resolve().parent/'python_packages'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.transforms import Bbox

root=Path('outputs/joint_candidates').resolve()
out=Path('outputs/joint_candidate_figures').resolve()
fd=out/'figures';td=out/'source_tables';sd=out/'scripts'
for p in [fd,td,sd]:p.mkdir(parents=True,exist_ok=True)
def rd(p):return pd.read_csv(p,dtype={'gene_id':str},low_memory=False)
def save(d,n):d.to_csv(td/n,index=False,encoding='utf-8-sig')
allg=rd(root/'source_tables/all_genome_wide_tests.csv.gz')
core=rd(root/'source_tables/all_core_candidates_with_sensitivity.csv')
points=rd(root/'source_tables/all_core_candidate_expression.csv.gz')
sens=rd(root/'source_tables/all_candidate_contrast_sensitivity.csv')
pair=rd(Path('outputs/Figure3_杂合型功能_20261022/Heterozygote_function_20261022/inputs/DESeq2_pairwise.csv.gz'))
T=['Spleen','Intestine'];F=['Baseline','Response'];G=['HomoDel','Hetero','HomoIns'];days=['d0','d1','d5']
gc={'HomoDel':'#397e85','Hetero':'#825ca6','HomoIns':'#d58a33'};gm={'HomoDel':'o','Hetero':'^','HomoIns':'s'}
fc={'Baseline':'#397e85','Response':'#c78632'};fm={'Baseline':'o','Response':'s'}
priority_color='#825ca6';other_color='#c99949';background='#bdc5cb';ink='#252d33'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7.2,'axes.labelsize':6.8,
 'axes.titlesize':7.5,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,'xtick.labelsize':6,'ytick.labelsize':6,
 'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','text.color':ink,'axes.labelcolor':ink})

# A display ordering across two separately corrected questions; no new combined FDR claim.
ranked=core.sort_values(['tissue','q_joint','p_joint','gene_id','family']).copy()
ranked['entry_order']=ranked.groupby('tissue').cumcount()+1
sel=ranked[ranked.finite_effect_priority_pass].drop_duplicates(['tissue','gene_id']).groupby('tissue',sort=False).head(5).copy()
sel['display_index']=sel.groupby('tissue').cumcount()+1
assert len(sel)==10 and sel.groupby('tissue').size().eq(5).all()
sel=sel.rename(columns={'family':'selection_question'})
qwide=allg.pivot(index=['tissue','gene_id'],columns='family',values='q_joint').rename(columns={'Baseline':'q_baseline','Response':'q_response'}).reset_index()
sel=sel.merge(qwide,on=['tissue','gene_id'],validate='one_to_one')
save(sel,'C_display_genes_fixed_rule.csv')
save(ranked,'B_complete_candidate_order.csv')
aa=allg[['tissue','family','gene_id','p_D','p_I','p_joint','q_joint','core']].merge(
 core[['tissue','family','gene_id','finite_effect_priority_pass']],on=['tissue','family','gene_id'],how='left',validate='one_to_one')
aa['priority']=aa.finite_effect_priority_pass.eq(True)
cap=25.
for ref in ['D','I']:
 aa['minus_log10_p_'+ref]=-np.log10(np.maximum(aa['p_'+ref],1e-300))
 aa['display_'+ref]=aa['minus_log10_p_'+ref].clip(upper=cap)
save(aa,'A_all_gene_test_points.csv.gz')
selected_points=points.merge(sel[['tissue','gene_id','display_index','selection_question']],on=['tissue','gene_id'],validate='many_to_one')
selected_points['log2_normalized_count_plus1']=np.log2(selected_points.normalized_count+1)
save(selected_points,'C_individual_fish_expression.csv')
ef=pair.merge(sel[['tissue','gene_id','display_index']],on=['tissue','gene_id'],validate='many_to_one')
ef=ef.merge(sens,on=['tissue','gene_id','stage','comparator'],validate='one_to_one')
save(ef,'D_effects_Wald_intervals_and_deletion_ranges.csv')
cutoffs={family:float(allg[(allg.family==family)&allg.core].p_joint.max()) for family in F}
shortnames={'127510871':'plexin-C1-like','127510740':'FGD4b-like','127511380':'DEF6-like','127498082':'NLRP12-like annotation',
 '127497534':'malic enzyme 1','127502293':'fos homolog Ab','127505992':'Uncharacterized','127497599':'Uncharacterized',
 '127494715':'Uncharacterized homolog','127502635':'Jun dimerization protein 2b'}
def qfmt(q):
 if not np.isfinite(q):return 'NA'
 return f'{q:.1e}' if q<.001 else f'{q:.2g}'
W,H=210/25.4,340/25.4
fig=plt.figure(figsize=(W,H),facecolor='white')
def txt(x,y,s,size=7,**kw):return fig.text(x,y,s,fontsize=size,**kw)
def heading(letter,title,y):
 txt(.025,y,letter,12,fontweight='bold',va='top');txt(.075,y-.001,title,9,fontweight='bold',va='top')
txt(.5,.992,'Genome-wide prioritization of heterozygote-associated genes',10.6,ha='center',va='top',fontweight='bold')

heading('A','Joint evidence from comparisons with both homozygotes',.964)
scatter_axes=[]
contexts=[('Spleen','Baseline'),('Spleen','Response'),('Intestine','Baseline'),('Intestine','Response')]
for i,(tis,fam) in enumerate(contexts):
 x=.085+i*.231;ax=fig.add_axes([x,.807,.177,.117]);scatter_axes.append(ax)
 d=aa[(aa.tissue==tis)&(aa.family==fam)&np.isfinite(aa.p_joint)]
 for mask,col,sz,alpha,z in [(~d.core,background,2.0,.42,1),(d.core&~d.priority,other_color,7,.85,2),(d.priority,priority_color,9,.95,3)]:
  zdata=d[mask];ax.scatter(zdata.display_D,zdata.display_I,s=sz,c=col,alpha=alpha,lw=0,zorder=z,rasterized=(z==1))
 threshold=-np.log10(cutoffs[fam])
 ax.plot([threshold,threshold,cap],[cap,threshold,threshold],color='#6e7881',ls='--',lw=.6,zorder=0)
 ax.set(xlim=(-.6,cap+.6),ylim=(-.6,cap+.6),xticks=[0,10,20],yticks=[0,10,20])
 ax.tick_params(length=2,pad=2)
 ax.set_title(f'{tis} | {fam}',fontweight='bold',fontsize=7.4,pad=5)
 ax.set_xlabel('−log10 P (H−Del)',labelpad=2,fontsize=6.3)
 if i==0:ax.set_ylabel('−log10 P (H−Ins)',labelpad=2,fontsize=6.3)
 ax.text(.98,.04,f'{int(d.core.sum())} joint / {int(d.priority.sum())} priority',transform=ax.transAxes,ha='right',va='bottom',fontsize=5.8,
         bbox=dict(facecolor='white',edgecolor='none',alpha=.87,pad=1.3))
 ax.text(.03,.97,f'n = {len(d):,}',transform=ax.transAxes,ha='left',va='top',fontsize=5.8,color='#68747b',
         bbox=dict(facecolor='white',edgecolor='none',alpha=.87,pad=1))
handles_a=[Line2D([],[],marker='o',ls='',color=background,ms=3,label='Other tested genes'),
 Line2D([],[],marker='o',ls='',color=other_color,ms=3,label='Other joint candidates'),
 Line2D([],[],marker='o',ls='',color=priority_color,ms=3,label='Passed finite-effect checks')]
fig.legend(handles=handles_a,loc='upper center',bbox_to_anchor=(.51,.778),ncol=3,fontsize=6.05,columnspacing=1.1,handletextpad=.4)
txt(.5,.753,'Dashed corner: joint BH discovery boundary. P-value axes capped at 25; all values are retained in source tables.',6.0,ha='center')

heading('B','Ordering of all jointly detected candidates',.730)
for tis,x in zip(T,[.10,.59]):
 ax=fig.add_axes([x,.562,.345,.125]);d=ranked[ranked.tissue==tis]
 for fam in F:
  z=d[d.family==fam]
  for status in [False,True]:
   a=z[z.finite_effect_priority_pass==status]
   ax.scatter(a.entry_order,-np.log10(a.q_joint),s=12 if status else 14,marker=fm[fam],
      facecolors=fc[fam] if status else 'white',edgecolors=fc[fam],lw=.65,zorder=3 if status else 2)
 for _,r in sel[sel.tissue==tis].iterrows():
  xx=r.entry_order;yy=-np.log10(r.q_joint)
  ax.scatter(xx,yy,s=27,facecolors='none',edgecolors=ink,lw=.6,zorder=5)
  offsets={'Spleen':[(4,8),(-7,15),(7,-13),(14,17),(15,-12)],
           'Intestine':[(9,8),(14,14),(29,12),(33,-9),(48,-18)]}
  dx,dy=offsets[tis][int(r.display_index)-1]
  ax.annotate(str(r.display_index),(xx,yy),xytext=(dx,dy),textcoords='offset points',fontsize=6.8,fontweight='bold',
     color=ink,ha='center',bbox=dict(boxstyle='round,pad=.05',facecolor='white',edgecolor='none',alpha=.92),
     arrowprops=dict(arrowstyle='-',color='#8d979e',lw=.45,shrinkA=1,shrinkB=3))
 ax.axhline(-np.log10(.05),color='#aab4ba',ls='--',lw=.5,zorder=0)
 ax.set_xlim(.4,len(d)+max(1,len(d)*.025));ax.set_ylim(1.1,16.2)
 ax.xaxis.set_major_locator(MaxNLocator(4,integer=True));ax.set_yticks([2,6,10,14]);ax.tick_params(length=2,pad=2)
 ax.set_xlabel('Order by joint q within tissue',fontsize=6.7,labelpad=2)
 if tis=='Spleen':ax.set_ylabel('−log10 joint q',labelpad=2)
 ax.set_title(f'{tis} | {len(d)} entries; {d.gene_id.nunique()} genes',fontsize=7.7,fontweight='bold',pad=5)
handles_b=[Line2D([],[],marker=fm[f],ls='',color=fc[f],ms=3,label=f) for f in F]+[
 Line2D([],[],marker='o',ls='',mec='#68747b',mfc='white',ms=3,label='Other joint candidates')]
fig.legend(handles=handles_b,loc='upper center',bbox_to_anchor=(.52,.539),ncol=3,fontsize=6.2,columnspacing=1.4,handletextpad=.5)
txt(.5,.515,'Numbers identify the five displayed distinct genes per tissue. Baseline and response retain separate BH families.',6.0,ha='center')

heading('C','Expression of the top-ranked distinct genes passing checks',.495)
gh=[Line2D([],[],color=gc[g],marker=gm[g],lw=.8,ms=3,label=g) for g in G]
fig.legend(handles=gh,loc='upper right',bbox_to_anchor=(.983,.500),ncol=3,fontsize=6.05,columnspacing=.7,handlelength=1)
expression_axes=[]
def expression(ax,tis,gene,legend=False):
 p=selected_points[(selected_points.tissue==tis)&(selected_points.gene_id==gene)]
 if p.empty:
  p=points[(points.tissue==tis)&(points.gene_id==gene)].copy();p['log2_normalized_count_plus1']=np.log2(p.normalized_count+1)
 for gi,g in enumerate(G):
  means=[];sd=[];offset=(gi-1)*.095
  for day,xx in zip(days,[0,1,5]):
   vals=p[(p.genotype==g)&(p.time==day)].sort_values('fish_id').log2_normalized_count_plus1.to_numpy()
   assert len(vals)==3
   ax.scatter(xx+offset+np.array([-.018,0,.018]),vals,s=7,c=gc[g],marker=gm[g],alpha=.58,lw=0,zorder=2)
   means.append(vals.mean());sd.append(vals.std(ddof=1))
  ax.errorbar(np.array([0,1,5])+offset,means,yerr=sd,color=gc[g],marker=gm[g],lw=.85,ms=2.5,capsize=1.1,elinewidth=.6,zorder=3)
 ax.set_xlim(-.45,5.45);ax.set_xticks([0,1,5]);ax.yaxis.set_major_locator(MaxNLocator(3));ax.tick_params(length=2,pad=2)
for ri,tis in enumerate(T):
 rows=sel[sel.tissue==tis].sort_values('display_index')
 ya=.348-ri*.148
 for ci,(_,r) in enumerate(rows.iterrows()):
  x=.10+ci*.176;ax=fig.add_axes([x,ya,.132,.089]);expression_axes.append(ax)
  expression(ax,tis,r.gene_id)
  if ri==1:ax.set_xlabel('Day',fontsize=6.2,labelpad=1)
  if ci==0:ax.set_ylabel('log2(normalized count + 1)',fontsize=6.0,labelpad=3)
  txt(x,ya+.130,f'[{int(r.display_index)}] {r.gene_id}',7.0,fontweight='bold')
  txt(x,ya+.116,shortnames.get(r.gene_id,'Uncharacterized'),6.0)
  txt(x,ya+.102,f'B q={qfmt(r.q_baseline)}; R q={qfmt(r.q_response)}',5.75)
 txt(.023,ya+.044,tis,8,rotation=90,ha='center',va='center',fontweight='bold')
txt(.5,.166,'Individual fish and mean ± SD; n = 3 per genotype × day. B: baseline; R: joint response across days 1 and 5.',6.1,ha='center')

heading('D','Effects and individual-fish sensitivity',.148)
dh=[Line2D([],[],color=gc[r],marker='o' if r=='HomoDel' else 's',lw=.7,ms=3,label='H−'+('Del' if r=='HomoDel' else 'Ins')) for r in ['HomoDel','HomoIns']]
dh.append(Line2D([],[],color='#b9c1c5',lw=3,label='Leave-one-fish range'))
fig.legend(handles=dh,loc='upper right',bbox_to_anchor=(.985,.153),ncol=3,fontsize=6.0,columnspacing=.7,handlelength=1.2)
stages=['Baseline','Response_d1','Response_d5']
effect_limits={}
for stage in stages:
 z=ef[ef.stage==stage]
 lo=min(z.ci95_low.min(),z.deletion_min.min());hi=max(z.ci95_high.max(),z.deletion_max.max());pad=.08*(hi-lo)
 effect_limits[stage]=(float(lo-pad),float(hi+pad))
for tis,x0 in zip(T,[.14,.65]):
 rows=sel[sel.tissue==tis].sort_values('display_index')
 for j,stage in enumerate(stages):
  ax=fig.add_axes([x0+j*.113,.037,.094,.071])
  for k,(_,r) in enumerate(rows.iterrows()):
   for ref,off,marker in [('HomoDel',-.13,'o'),('HomoIns',.13,'s')]:
    z=ef[(ef.tissue==tis)&(ef.gene_id==r.gene_id)&(ef.stage==stage)&(ef.comparator==ref)].iloc[0]
    ax.plot([z.deletion_min,z.deletion_max],[k+off,k+off],color='#b9c1c5',lw=3.4,alpha=.65,zorder=1)
    ax.errorbar(z.effect,k+off,xerr=np.array([[z.effect-z.ci95_low],[z.ci95_high-z.effect]]),fmt=marker,color=gc[ref],mec=gc[ref],mfc=gc[ref],
      ms=2.6,mew=.55,capsize=1,elinewidth=.65,zorder=3)
  ax.axvline(0,c='#adb6bc',ls='--',lw=.5,zorder=0)
  ax.set_ylim(4.55,-.55);ax.set_xlim(effect_limits[stage])
  ticks=MaxNLocator(3).tick_values(*effect_limits[stage]);ax.set_xticks(ticks[(ticks>=effect_limits[stage][0])&(ticks<=effect_limits[stage][1])])
  ax.tick_params(length=2,pad=2,labelsize=5.7)
  ax.set_yticks(range(5),[f'[{int(r.display_index)}] {r.gene_id}' for r in rows.itertuples()] if j==0 else ['']*5,fontsize=5.75)
  ax.set_title(['Baseline','Day 1 response','Day 5 response'][j],fontsize=6.2,pad=2)
 txt(x0+.16,.126,tis,7.5,fontweight='bold',ha='center')
 txt(x0+.16,.020,'Hetero−homozygote effect (log2 units)',6.0,ha='center')
txt(.5,.007,'Pointwise 95% Wald intervals; gray segments: range of 27 fixed-parameter leave-one-fish estimates.',6.0,ha='center')

fig.canvas.draw()
renderer=fig.canvas.get_renderer();outside=[]
for artist in fig.findobj(matplotlib.text.Text):
 if not artist.get_visible() or not artist.get_text():continue
 bb=artist.get_window_extent(renderer)
 if bb.x0<-1 or bb.y0<-1 or bb.x1>fig.bbox.width+1 or bb.y1>fig.bbox.height+1:outside.append(artist.get_text())
def export(figure,name):
 figure.savefig(fd/(name+'.png'),dpi=300,facecolor='white');figure.savefig(fd/(name+'.pdf'),facecolor='white');figure.savefig(fd/(name+'.svg'),facecolor='white')
export(fig,'Figure4_genome_wide_candidates')
fig.savefig(fd/'Figure4_preview.png',dpi=150,facecolor='white')
regions={'Figure4A_joint_discovery':(.745,.975),'Figure4B_complete_ordering':(.507,.741),'Figure4C_top_gene_expression':(.159,.505),'Figure4D_effect_sensitivity':(0,.157)}
for name,(lo,hi) in regions.items():
 box=Bbox.from_bounds(0,lo*H,W,(hi-lo)*H)
 fig.savefig(fd/(name+'.png'),dpi=300,bbox_inches=box,facecolor='white')
 fig.savefig(fd/(name+'.pdf'),bbox_inches=box,facecolor='white')
plt.close(fig)

# Full candidate catalogue: every jointly significant gene once per tissue.
catalogue=ranked.drop_duplicates(['tissue','gene_id']).copy()
save(catalogue,'Supplement_all_unique_candidate_order.csv')
page_rows=[]
with PdfPages(fd/'Supplement_all_candidate_expression.pdf') as pdf:
 page=0
 for tis in T:
  c=catalogue[catalogue.tissue==tis]
  for start in range(0,len(c),16):
   page+=1;chunk=c.iloc[start:start+16]
   ff,axes=plt.subplots(4,4,figsize=(210/25.4,265/25.4));ff.subplots_adjust(left=.08,right=.975,bottom=.065,top=.88,hspace=.95,wspace=.50)
   ff.suptitle(f'All jointly detected candidates | {tis} | {start+1}–{start+len(chunk)} of {len(c)}',fontsize=10,fontweight='bold',y=.976)
   ff.legend(handles=gh,loc='upper center',bbox_to_anchor=(.52,.952),ncol=3,fontsize=7)
   for ax,(_,r) in zip(axes.flat,chunk.iterrows()):
    expression(ax,tis,r.gene_id);ax.set_xlabel('Day',fontsize=6.3,labelpad=2)
    qq=qwide[(qwide.tissue==tis)&(qwide.gene_id==r.gene_id)].iloc[0]
    ann=r.annotation if isinstance(r.annotation,str) and r.annotation!='-' else 'Uncharacterized'
    ann=textwrap.shorten(ann,width=30,placeholder='…')
    status='Best-q context: checks passed' if r.finite_effect_priority_pass else 'Best-q context: see flags'
    ax.set_title(f'{r.gene_id}\n{ann}\nB q={qfmt(qq.q_baseline)}; R q={qfmt(qq.q_response)}\n{status}',fontsize=6.1,loc='left',pad=5)
    page_rows.append({'page':page,'tissue':tis,'gene_id':r.gene_id,'q_baseline':qq.q_baseline,'q_response':qq.q_response})
   for ax in axes.flat[len(chunk):]:ax.axis('off')
   ff.text(.017,.49,'log2(normalized count + 1)',fontsize=7,rotation=90,ha='center',va='center')
   ff.text(.5,.020,'Individual fish and mean ± SD; n = 3 per genotype × day. Different fish were sampled at each day.',fontsize=6.4,ha='center')
   pdf.savefig(ff);plt.close(ff)
save(pd.DataFrame(page_rows),'Supplement_expression_page_index.csv')

qa={'main_size_mm':[210,340],'main_dpi':300,'tested_points':len(aa),'finite_joint_tests':int(np.isfinite(aa.p_joint).sum()),
 'core_gene_question_entries':len(core),'priority_gene_question_entries':int(core.finite_effect_priority_pass.sum()),
 'selected_distinct_genes_per_tissue':sel.groupby('tissue').size().to_dict(),
 'selected_questions':sel.groupby(['tissue','selection_question']).size().rename('genes').reset_index().to_dict('records'),
 'main_individual_expression_points':len(selected_points),'main_effect_estimates':len(ef),
 'supplement_gene_tissue_entries':len(catalogue),'supplement_pages':page,'P_axis_cap':cap,'joint_P_discovery_cutoffs':cutoffs,
 'text_outside_canvas':outside,'matplotlib':matplotlib.__version__,'no_new_significance_tests':True,
 'selection_rule':'All finite-effect-priority entries sorted within tissue by joint q, joint P, gene ID, question; first five distinct genes. Display cap only, no stage or annotation quota.',
 'ordering_is_descriptive_across_two_separate_FDR_families':True}
(out/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(Path(__file__),sd/Path(__file__).name)
shutil.copy2(Path(__file__),fd/'gen_fig_4.py')
print(json.dumps(qa,ensure_ascii=False,indent=2))
