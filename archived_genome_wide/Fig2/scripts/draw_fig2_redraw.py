"""Draw revised Figure 2 exclusively from exported source tables."""
from pathlib import Path
import sys, argparse, json, textwrap
sys.path.insert(0, str(Path(__file__).resolve().parent / 'python_packages'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox

ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();out=args.output;tab=out/'source_tables'
figdir=out/'figures';figdir.mkdir(exist_ok=True)
def read(name): return pd.read_csv(tab/name,dtype={'gene_id':str})
order=read('B_gene_order_and_cluster.csv.gz')
heat=read('B_group_mean_VST_z.csv.gz')
profiles=read('B_cluster_response_profiles.csv')
selection=read('C_gene_selection.csv')
points=read('C_individual_points.csv')
ora=read('D_KEGG_ORA_all.csv')
background=read('D_background_and_detection_summary.csv')
counts=read('A_venn_counts.csv')
G=['HomoDel','Hetero','HomoIns'];T=['Spleen','Intestine'];days=['d0','d1','d5']
colors={'HomoDel':'#397e85','Hetero':'#825ca6','HomoIns':'#d58a33'}
markers={'HomoDel':'o','Hetero':'^','HomoIns':'s'}
cluster_colors=['#7c9cb8','#aac6ac','#c0a8ca','#dec19e']
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':7.8,'axes.titlesize':8,'axes.labelsize':7.5,'xtick.labelsize':6.5,
    'ytick.labelsize':6.5,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,
    'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'text.color':'#20262c','axes.labelcolor':'#20262c','xtick.color':'#30383f','ytick.color':'#30383f'})
cmap=LinearSegmentedColormap.from_list('expression',['#34648b','#f7f7f6','#bb534c'])
enrich_cmap=LinearSegmentedColormap.from_list('enrichment',['#c4dfdc','#23727b'])
enrich_norm=Normalize(1.3,14)
W,H=210/25.4,365/25.4
fig=plt.figure(figsize=(W,H),facecolor='white')

def text(x,y,s,size=8,**kw):
    return fig.text(x,y,s,fontsize=size,**kw)
def heading(letter,title,y):
    text(.025,y,letter,12,fontweight='bold',va='top')
    text(.068,y-.001,title,9.3,fontweight='bold',va='top')
def qfmt(q):
    if q>=.01: return f'{q:.3f}'
    if q>=.001: return f'{q:.4f}'
    return f'{q:.1e}'

text(.5,.990,'Baseline differences and genotype-dependent responses',11,
     ha='center',va='top',fontweight='bold')
heading('A','Overlap of baseline and response differences',.964)
for tissue,x in zip(T,[.065,.555]):
    ax=fig.add_axes([x,.850,.39,.098]);ax.set_aspect('equal');ax.axis('off')
    ax.set(xlim=(-2.15,2.15),ylim=(-1.7,1.65))
    c=counts[counts.tissue.eq(tissue)].set_index('category').genes
    b,r,both=int(c['Baseline only']),int(c['Response only']),int(c['Both'])
    left=Circle((-.65,0),1.05,facecolor='#deefed',edgecolor='#579c9d',lw=1)
    right=Circle((.65,0),1.05,facecolor='#fff0d8',edgecolor='#cf9b54',lw=1)
    ax.add_patch(left);ax.add_patch(right)
    lens=Circle((-.65,0),1.05,facecolor='#d9cae8',edgecolor='none')
    lens.set_clip_path(right);ax.add_patch(lens)
    for center,col in [((-.65,0),'#579c9d'),((.65,0),'#cf9b54')]:
        ax.add_patch(Circle(center,1.05,fill=False,edgecolor=col,lw=1))
    ax.text(0,1.48,tissue,ha='center',va='center',fontsize=9,fontweight='bold')
    ax.text(-1.05,1.17,'Baseline',ha='center',fontsize=7.8,color='#357f80')
    ax.text(1.05,1.17,'Response',ha='center',fontsize=7.8,color='#a87231')
    for xx,n,lab in [(-1.02,b,'Baseline\nonly'),(0,both,'Both'),(1.02,r,'Response\nonly')]:
        ax.text(xx,.16,f'{n:,}',ha='center',fontsize=10,fontweight='bold')
        ax.text(xx,-.30,lab,ha='center',va='center',fontsize=6.6,linespacing=1.05)
    ax.text(0,-1.25,f'Baseline total: {b+both:,}   |   Response total: {r+both:,}',ha='center',fontsize=6.8)
    ax.text(0,-1.58,f"Neither: {int(c['Neither']):,}   |   Tested: {b+r+both+int(c['Neither']):,}",ha='center',fontsize=6.8)

heading('B','Expression patterns of response-associated genes',.827)
handles=[Line2D([],[],color=colors[g],marker=markers[g],lw=1,ms=3,label=g) for g in G]
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.974,.841),ncol=3,
           fontsize=6.7,handlelength=1.4,columnspacing=.9)
heat_images=[]
for tissue,x in zip(T,[.07,.56]):
    rows=order[order.tissue.eq(tissue)].sort_values('row_index')
    hh=heat[heat.tissue.eq(tissue)].copy();hh['key']=hh.time+'_'+hh.genotype
    cols=[d+'_'+g for d in days for g in G]
    matrix=hh.pivot(index='gene_id',columns='key',values='mean_VST_z').loc[rows.gene_id,cols].to_numpy()
    ax=fig.add_axes([x,.591,.222,.201])
    im=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=-2,vmax=2,interpolation='nearest',rasterized=True)
    heat_images.append(im);ax.set_yticks([])
    ax.set_xticks(range(9),['Del','Het','Ins']*3,fontsize=6.1)
    ax.tick_params(axis='x',length=0,pad=3)
    for line in [2.5,5.5]: ax.axvline(line,color='white',lw=1.5)
    for xx,day in zip([1,4,7],[0,1,5]):
        ax.text(xx,1.012,f'Day {day}',transform=ax.get_xaxis_transform(),ha='center',fontsize=7)
    for spine in ax.spines.values():spine.set_visible(False)
    for cluster,rr in rows.groupby('cluster',sort=True):
        start=int(rr.row_index.min())-1;end=int(rr.row_index.max())
        ax.axhline(end-.5,color='white',lw=.8)
        ax.add_patch(Rectangle((-.078,start-.5),.047,len(rr),transform=ax.get_yaxis_transform(),
                               clip_on=False,color=cluster_colors[int(cluster)-1],lw=0))
    text(x+.11,.808,f'{tissue} | {len(rows):,} genes',8.5,ha='center',fontweight='bold')
    # Equal-height panels make small clusters visible even when heatmap bands are thin.
    for j,cluster in enumerate([1,2,3,4]):
        yy=.591+(3-j)*.0507
        pa=fig.add_axes([x+.283,yy+.005,.120,.0335])
        pp=profiles[(profiles.tissue==tissue)&(profiles.cluster==cluster)]
        for g in G:
            p=pp[pp.genotype.eq(g)].set_index('time').loc[days]
            pa.plot([0,1,5],p.mean_response_log2FC,color=colors[g],marker=markers[g],ms=2.4,lw=.85)
        pa.axhline(0,color='#acb1b4',lw=.45,ls='--',zorder=0)
        pa.set_xticks([0,1,5]);pa.tick_params(length=2,pad=1,labelsize=5.8)
        pa.yaxis.set_major_locator(MaxNLocator(nbins=2))
        pa.set_xlim(-.25,5.25)
        n=int(pp.genes.iloc[0])
        pa.set_title(f'C{cluster}  |  n = {n:,}',fontsize=6.6,pad=2,loc='left',
                     color='#30383f')
        if j==3:pa.set_xlabel('Day',fontsize=6.4,labelpad=1)
    text(x+.343,.806,'Mean response (log2FC)',6.4,ha='center')
cax=fig.add_axes([.285,.563,.19,.007])
cb=fig.colorbar(heat_images[0],cax=cax,orientation='horizontal',ticks=[-2,0,2])
cb.ax.tick_params(labelsize=6.3,length=2,pad=1)
cb.set_label('Mean VST row z-score',fontsize=6.8,labelpad=2)
text(.53,.565,'Heatmap: relative expression\nCurves: change from own day 0',6.5,va='center',linespacing=1.4)

heading('C','IL10 and representative response-associated genes',.523)
main=selection[selection.main_panel].copy()
mainlist=[]
for tissue in T:
    for role in ['Prespecified IL10','Response representative']:
        mainlist.append(main[(main.tissue==tissue)&(main.role==role)].iloc[0])
for k,s in enumerate(mainlist):
    x=.080+k*.235
    ax=fig.add_axes([x,.391,.180,.095])
    p=points[(points.tissue==s.tissue)&(points.gene_id==s.gene_id)]
    for gi,g in enumerate(G):
        pp=p[p.genotype==g]
        means=[];sds=[]
        offset=(gi-1)*.13
        for day,xx in zip(days,[0,1,5]):
            vals=pp[pp.time==day].sort_values('rep').log2_normalized_count_plus1.to_numpy()
            ax.scatter(xx+offset+np.array([-.022,0,.022]),vals,s=9,marker=markers[g],
                       color=colors[g],alpha=.65,lw=0,zorder=3)
            means.append(vals.mean());sds.append(vals.std(ddof=1))
        ax.errorbar(np.array([0,1,5])+offset,means,yerr=sds,color=colors[g],
                    lw=.95,marker=markers[g],ms=2.8,capsize=1.5,elinewidth=.7,zorder=4)
    ax.set_xlim(-.5,5.5);ax.set_xticks([0,1,5]);ax.set_xlabel('Day',labelpad=2)
    ax.yaxis.set_major_locator(MaxNLocator(4));ax.tick_params(length=2,pad=2)
    gene='IL10' if s.role=='Prespecified IL10' else f'Gene {s.gene_id} | C{int(s.cluster)}'
    text(x,.504,f'{s.tissue} | {gene}',7.5,fontweight='bold')
    text(x,.493,f'Baseline q = {qfmt(s.baseline_q)}; response q = {qfmt(s.response_q)}',6.3)
text(.026,.437,'log2(normalized count + 1)',7.2,rotation=90,ha='center',va='center')
text(.5,.361,'Individual fish and mean ± SD; n = 3 per genotype × day',6.7,ha='center')

heading('D','KEGG enrichment of baseline and response gene sets',.343)
# Display ALL significant pathways, including all 13 intestine response pathways.
for tissue,x,colstart in zip(T,[.260,.755],[.068,.562]):
    for test,y,h,title_y in [('Baseline',.239,.072,.326),('Response',.053,.151,.222)]:
        d=ora[(ora.tissue==tissue)&(ora.test==test)&(ora.q_within<.05)].sort_values(['q_within','pvalue','pathway_id'])
        bg=background[(background.tissue==tissue)&(background.test==test)].iloc[0]
        ax=fig.add_axes([x,y,.205,h])
        yy=np.arange(len(d))
        ax.scatter(d.fold_enrichment,yy,s=8+.34*d.hits,c=-np.log10(d.q_within),
                   cmap=enrich_cmap,norm=enrich_norm,edgecolors='white',linewidths=.3,zorder=3)
        ax.set_yticks(yy,d.pathway.tolist(),fontsize=6.4)
        ax.set_ylim(len(d)-.45,-.55);ax.set_xlim(.95,4.05)
        ax.set_xticks([1,2,3,4]);ax.axvline(1,color='#bfc3c5',lw=.6,ls='--')
        ax.tick_params(axis='y',length=0,pad=4);ax.tick_params(axis='x',length=2,pad=2)
        ax.set_xlabel('Fold enrichment',fontsize=6.8,labelpad=2)
        for sp in ['left','top','right']:ax.spines[sp].set_visible(False)
        text(colstart,title_y,f'{tissue} | {test} ({len(d)} pathways)',8,fontweight='bold')
        text(colstart,title_y-.011,f'{int(bg.annotated_targets):,} annotated targets / {int(bg.annotated_background):,} background genes',6.4)
cbax=fig.add_axes([.135,.015,.18,.006])
cb=fig.colorbar(plt.cm.ScalarMappable(norm=enrich_norm,cmap=enrich_cmap),cax=cbax,orientation='horizontal',ticks=[2,6,10,14])
cb.ax.tick_params(labelsize=6.1,length=2,pad=1)
text(.225,.027,'−log10(q)',6.8,ha='center')
sizehandles=[plt.scatter([],[],s=8+.34*n,color='#6d9e9e',edgecolors='white',linewidths=.3,label=str(n)) for n in [25,75,150]]
fig.legend(handles=sizehandles,title='Hit genes',loc='lower center',bbox_to_anchor=(.60,.006),ncol=3,
           fontsize=6.5,title_fontsize=6.7,handletextpad=.6,columnspacing=1)
text(.98,.008,'All q < 0.05 pathways shown',6.3,ha='right')

def export(figure,name):
    figure.savefig(figdir/(name+'.png'),dpi=300,facecolor='white')
    figure.savefig(figdir/(name+'.pdf'),facecolor='white',metadata={'Title':name,'Creator':'Matplotlib; archived model outputs'})

export(fig,'Figure2_revised')
fig.savefig(figdir/'Figure2_revised.svg',facecolor='white')
regions={'Figure2A_venn':(.845,.975),'Figure2B_patterns':(.536,.840),
         'Figure2C_IL10_and_examples':(.353,.535),'Figure2D_all_KEGG':(0,.350)}
for name,(lo,hi) in regions.items():
    bbox=Bbox.from_bounds(0,lo*H,W,(hi-lo)*H)
    fig.savefig(figdir/(name+'.png'),dpi=300,bbox_inches=bbox,facecolor='white')
    fig.savefig(figdir/(name+'.pdf'),bbox_inches=bbox,facecolor='white')
plt.close(fig)

# Full-sample expression heatmaps preserve individual biological replicates.
ind=read('B_individual_VST_z.csv.gz')
md=pd.read_csv(out/'inputs'/'sample_metadata.csv')
fs=plt.figure(figsize=(W,185/25.4))
for i,tissue in enumerate(T):
    rows=order[order.tissue==tissue].sort_values('row_index')
    ss=[]
    for day in days:
        for g in G:
            ss+=md[(md.tissue==tissue)&(md.time==day)&(md.genotype==g)].sort_values('rep')['sample'].tolist()
    matrix=ind[ind.tissue==tissue].pivot(index='gene_id',columns='sample',values='row_z').loc[rows.gene_id,ss].to_numpy()
    ax=fs.add_axes([.12,.56-i*.46,.83,.33])
    im=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=-2,vmax=2,interpolation='nearest',rasterized=True)
    ax.set_xticks(range(27),ss,rotation=90,fontsize=6);ax.set_yticks([])
    ax.set_title(f'{tissue} | {len(rows):,} response genes | all 27 samples',fontsize=9,fontweight='bold')
    for c,rr in rows.groupby('cluster'):
        ax.axhline(rr.row_index.max()-.5,color='white',lw=.7)
        ax.text(-.02,rr.row_index.mean()-1,f'C{c}',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6.3)
    for x in [8.5,17.5]:ax.axvline(x,color='white',lw=1.5)
    for sp in ax.spines.values():sp.set_visible(False)
fs.text(.5,.98,'Individual-sample VST row z-scores',ha='center',fontsize=11,fontweight='bold')
cax=fs.add_axes([.36,.025,.28,.012]);fs.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2])
export(fs,'Figure2_supplement_individual_heatmaps');plt.close(fs)

# Every cluster has a representative in the supplemental panel.
rs=selection[selection.role=='Response representative'].sort_values(['tissue','cluster'])
fc,axes=plt.subplots(2,4,figsize=(W,120/25.4))
fc.subplots_adjust(left=.09,right=.98,bottom=.13,top=.86,hspace=.70,wspace=.48)
for ax,(_,s) in zip(axes.flat,rs.iterrows()):
    pp=points[(points.tissue==s.tissue)&(points.gene_id==s.gene_id)]
    for gi,g in enumerate(G):
        mm=[];sd=[];offset=(gi-1)*.12
        for day,xx in zip(days,[0,1,5]):
            v=pp[(pp.genotype==g)&(pp.time==day)].sort_values('rep').log2_normalized_count_plus1.to_numpy()
            ax.scatter(xx+offset+np.array([-.022,0,.022]),v,s=8,color=colors[g],marker=markers[g],alpha=.6,lw=0)
            mm.append(v.mean());sd.append(v.std(ddof=1))
        ax.errorbar(np.array([0,1,5])+offset,mm,yerr=sd,color=colors[g],lw=.8,ms=2.2,marker=markers[g],capsize=1)
    ax.set_title(f'{s.tissue} C{int(s.cluster)}\nGene {s.gene_id}\nResponse q = {qfmt(s.response_q)}',fontsize=7)
    ax.set_xticks([0,1,5]);ax.set_xlim(-.5,5.5);ax.set_xlabel('Day',fontsize=7)
fc.text(.018,.5,'log2(normalized count + 1)',rotation=90,ha='center',va='center',fontsize=8)
fc.legend(handles=handles,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.0),fontsize=7)
export(fc,'Figure2_supplement_cluster_examples');plt.close(fc)

qa={'main_size_mm':[210,365],'raster_dpi':300,
    'all_significant_KEGG_pathways_displayed':True,'KEGG_display_rows':int((ora.q_within<.05).sum()),
    'heatmap_group_columns':9,'heatmap_individual_samples_per_tissue':27,
    'IL10_tissues':['Spleen','Intestine'],'fold_cutoff':None,
    'matplotlib':matplotlib.__version__}
(out/'figure_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa,indent=2))
