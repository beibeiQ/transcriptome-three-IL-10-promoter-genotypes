from pathlib import Path
import pandas as pd,numpy as np,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[2];W=ROOT/'work/tissue_function_revision';I=W/'inputs';A=W/'analysis'
member=pd.read_csv(I/'pathway_membership.csv',dtype={'gene_id':str})
member=member[member.in_pooled_universe].copy()
ann=pd.read_csv(I/'gene_annotation.csv.gz',dtype={'gene_id':str}).fillna('')
pair=pd.read_csv(I/'gene_pairwise_comparisons.csv.gz',dtype={'gene_id':str},low_memory=False)
direct=pd.read_csv(I/'gene_direct_tissue_effects.csv.gz',dtype={'gene_id':str},low_memory=False)
assert not pair.duplicated(['gene_id','tissue','stage']).any()
assert not direct.duplicated(['gene_id','comparator','stage']).any()
long=[]
for comp,code in [('HomoDel','D'),('HomoIns','I')]:
    x=pair[['tissue','stage','gene_id','strict','baseMean','gene_name','primary_annotation',f'effect_{code}',f'lfcSE_{code}',f'pvalue_{code}',f'q_pair_family_{code}',f'ci95_low_{code}',f'ci95_high_{code}']].copy()
    x.columns=['tissue','stage','gene_id','strict','baseMean','gene_name','primary_annotation','effect','SE','P','Q','CI_low','CI_high'];x['comparator']=comp;long.append(x)
pl=pd.concat(long,ignore_index=True)
joined=member[['pathway_id','pathway','domain','gene_id']].merge(pl,on='gene_id',how='left')
joined.to_csv(A/'complete_pathway_gene_contrasts.csv.gz',index=False)
stats=[]
for keys,d in joined.groupby(['pathway_id','pathway','tissue','stage','comparator'],dropna=False):
    good=d[np.isfinite(d.effect)]
    stats.append(dict(zip(['pathway_id','pathway','tissue','stage','comparator'],keys))|dict(mapped_genes=d.gene_id.nunique(),finite_genes=len(good),median_effect=good.effect.median(),fraction_positive=(good.effect>0).mean(),fraction_negative=(good.effect<0).mean(),positive_detected=((good.Q<.05)&(good.effect>0)).sum(),negative_detected=((good.Q<.05)&(good.effect<0)).sum()))
pd.DataFrame(stats).to_csv(A/'pathway_gene_direction_summary.csv',index=False)
dj=member[['pathway_id','pathway','domain','gene_id']].merge(direct,on='gene_id',how='left').merge(ann[['gene_id','gene_name','primary_description']],on='gene_id',how='left')
dj.to_csv(A/'complete_pathway_gene_direct_effects.csv.gz',index=False)
dstat=[]
for keys,d in dj.groupby(['pathway_id','pathway','stage','comparator'],dropna=True):
    good=d[np.isfinite(d.effect)]
    dstat.append(dict(zip(['pathway_id','pathway','stage','comparator'],keys))|dict(tested_genes=len(good),median_direct_effect=good.effect.median(),positive_detected=((good.q_family<.05)&(good.effect>0)).sum(),negative_detected=((good.q_family<.05)&(good.effect<0)).sum()))
pd.DataFrame(dstat).to_csv(A/'pathway_gene_direct_summary.csv',index=False)
rule={'two_pathway_examples':[{'pathway_id':'cide04142','stage':'Response_d1','comparator':'HomoDel','expected_direct_sign':-1},{'pathway_id':'cide04210','stage':'Response_d5','comparator':'HomoIns','expected_direct_sign':1}],
      'fixed_reference_gene':'127497534'}
selected=[];used=set()
strict=pair.groupby(['gene_id','tissue']).strict.first().unstack().fillna(False)
both_strict=set(strict.index[strict.all(axis=1)])
for s in rule['two_pathway_examples']:
    z=dj[(dj.pathway_id==s['pathway_id'])&(dj.stage==s['stage'])&(dj.comparator==s['comparator'])].copy()
    z=z[z.gene_id.isin(both_strict)&(~z.gene_id.isin(used))&(np.sign(z.effect)==s['expected_direct_sign'])]
    z=z[z.primary_description.fillna('').str.len().gt(3)&~z.primary_description.fillna('').str.contains('uncharacterized|hypothetical',case=False)]
    z=z.sort_values(['q_family','gene_id']);assert len(z)>0
    row=z.iloc[0];selected.append(row);used.add(row.gene_id)
sel=pd.DataFrame(selected);sel.to_csv(A/'illustrative_pathway_genes.csv',index=False)
focus=list(sel.gene_id)+['127497534']
direct[direct.gene_id.isin(focus)].merge(ann,on='gene_id',how='left').to_csv(A/'focused_gene_direct_effects.csv',index=False)
pl[pl.gene_id.isin(focus)].to_csv(A/'focused_gene_pairwise_effects.csv',index=False)
norm=pd.read_csv(I/'normalized_counts_long.csv.gz',dtype={'gene_id':str})
norm[norm.gene_id.isin(focus)].to_csv(A/'focused_gene_individual_counts.csv',index=False)
e=pd.read_csv(A/'program_contrasts.csv');om=pd.read_csv(A/'program_interactions.csv')
e=e[(e.tissue=='Intestine')&(e.question=='Response')];qual=[]
for (pid,time),d in e.groupby(['pathway_id','time']):
    if len(d)==2 and (d.Q<.05).all() and np.sign(d.estimate).nunique()==1:
        qual.append({'pathway_id':pid,'time':time,'max_Q':d.Q.max(),'median_abs_effect':d.estimate.abs().median()})
if qual:
    selected_programs=pd.DataFrame(qual).sort_values(['max_Q','median_abs_effect'],ascending=[True,False]).drop_duplicates('pathway_id').head(2);mode='Bilateral same-stage same-direction support'
else:
    selected_programs=om[om.tissue=='Intestine'].sort_values(['Q','pathway_id']).head(2).copy();mode='No bilateral score-response detections; two smallest intestinal omnibus Q values'
selected_programs['selection_rule']=mode;selected_programs.to_csv(A/'illustrative_program_selection.csv',index=False)
print('Selected pathway genes:');print(sel[['pathway_id','gene_id','primary_description','stage','comparator','effect','q_family']].to_string(index=False))
print('Selected intestinal program trajectories:');print(selected_programs.to_string(index=False))
ss=pd.DataFrame(stats);print('Intestinal selected-pathway gene directions:');print(ss[(ss.tissue=='Intestine')&(((ss.pathway_id=='cide04142')&(ss.stage=='Response_d1'))|((ss.pathway_id=='cide04210')&(ss.stage=='Response_d5')))].to_string(index=False))
print('Focused direct summaries:');print(pd.DataFrame(dstat).query("pathway_id in ['cide04142','cide04210']").to_string(index=False))
