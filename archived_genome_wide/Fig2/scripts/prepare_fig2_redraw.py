"""Prepare a new Figure 2 from archived fitted models; no model refitting."""
from pathlib import Path
import sys, argparse, json, shutil, hashlib
sys.path.insert(0, str(Path(__file__).resolve().parent / 'python_packages'))
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ap = argparse.ArgumentParser()
ap.add_argument('--source', type=Path)
ap.add_argument('--output', type=Path, required=True)
args = ap.parse_args()
out = args.output
inp, tables = out/'inputs', out/'source_tables'
inp.mkdir(parents=True, exist_ok=True)
tables.mkdir(exist_ok=True)
files = {
    'omnibus_tests.csv.gz': 'Figure2_20261008/omnibus_tests.csv.gz',
    'normalized_counts_long.csv.gz': 'Figure2_20261008/normalized_counts_long.csv.gz',
    'sample_metadata.csv': 'Figure2_20261008/sample_metadata.csv',
    'KEGG_annotation_input.tsv': 'Figure2_20261008/KEGG_annotation_input.tsv',
    'original_response_ORA.csv': 'Figure2_20261008/KEGG_ORA_all.csv',
    'row_order.csv.gz': 'Figure2_revision_20261009/Figure2B_row_order.csv.gz',
    'ordering_features.csv.gz': 'Figure2_revision_20261009/Figure2B_ordering_features.csv.gz',
    'own_baseline_effects.csv.gz': 'Figure2_revision_20261009/Figure2B_own_baseline_source.csv.gz',
    'Spleen_VST.csv.gz': 'Figure3_20261009/Spleen_VST.csv.gz',
    'Intestine_VST.csv.gz': 'Figure3_20261009/Intestine_VST.csv.gz',
}
if args.source:
    for name, relative in files.items():
        shutil.copy2(args.source/relative, inp/name)

def read(name, **kw):
    return pd.read_csv(inp/name, dtype={'gene_id': str}, **kw)

def save(df, name):
    df.to_csv(tables/name, index=False, encoding='utf-8-sig')

def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    ans = np.empty(len(p)); ans[order] = np.minimum(q,1)
    return ans

om, md = read('omnibus_tests.csv.gz'), read('sample_metadata.csv')
order, effects = read('row_order.csv.gz'), read('own_baseline_effects.csv.gz')
features = read('ordering_features.csv.gz')
assert not om.duplicated(['tissue','gene_id','test']).any()
assert md['sample'].is_unique and len(md)==54
assert (md.groupby(['tissue','genotype','time']).size()==3).all()
assert not order.duplicated(['tissue','gene_id']).any()
baseline = om[om.test.eq('Baseline')].set_index(['tissue','gene_id'])
response = om[om.test.eq('Response')].set_index(['tissue','gene_id'])
classes = baseline[['pvalue','q_within']].join(response[['pvalue','q_within']], lsuffix='_baseline', rsuffix='_response').reset_index()
classes['eligible'] = np.isfinite(classes.pvalue_baseline)&np.isfinite(classes.pvalue_response)
classes['category'] = np.select([
    ~classes.eligible,
    (classes.q_within_baseline<.05)&(classes.q_within_response<.05),
    classes.q_within_baseline<.05, classes.q_within_response<.05],
    ['Not testable','Both','Baseline only','Response only'], default='Neither')
save(classes, 'A_gene_membership.csv.gz')
counts = classes.groupby(['tissue','category'], sort=False).size().rename('genes').reset_index()
save(counts, 'A_venn_counts.csv')
expected = {'Spleen': {'Both':925,'Baseline only':1225,'Response only':906,'Neither':18502},
            'Intestine': {'Both':4273,'Baseline only':2669,'Response only':3272,'Neither':13001}}
for tissue, cc in expected.items():
    actual = counts[counts.tissue.eq(tissue)].set_index('category').genes.to_dict()
    for category, n in cc.items(): assert actual[category]==n, (tissue,category,actual)

genotypes, days = ['HomoDel','Hetero','HomoIns'], ['d0','d1','d5']
heat, individual, profile = [], [], []
for tissue in ['Spleen','Intestine']:
    rows = order[order.tissue.eq(tissue)].sort_values('row_index')
    ids = rows.gene_id.tolist()
    expected_ids = set(response.loc[tissue].query('q_within < .05').index)
    assert set(ids)==expected_ids
    samples = md[md.tissue.eq(tissue)]['sample'].tolist()
    vst = read(tissue+'_VST.csv.gz').set_index('gene_id').loc[ids, samples]
    E = vst.to_numpy()
    sd = E.std(axis=1, ddof=1)
    assert (sd>0).all() and np.isfinite(E).all()
    z = (E-E.mean(axis=1,keepdims=True))/sd[:,None]
    assert np.max(np.abs(z.mean(axis=1)))<1e-12
    for j,sample in enumerate(samples):
        m = md.set_index('sample').loc[sample]
        zz = rows.copy()
        zz['sample'], zz['genotype'], zz['time'] = sample, m.genotype, m.time
        zz['VST'],zz['row_z'] = E[:,j],z[:,j]
        individual.append(zz)
    for day in days:
        for genotype in genotypes:
            ss=md[(md.tissue==tissue)&(md.time==day)&(md.genotype==genotype)]['sample']
            indices=[samples.index(s) for s in ss]
            zz=rows.copy();zz['genotype']=genotype;zz['time']=day
            zz['mean_VST_z']=z[:,indices].mean(axis=1)
            heat.append(zz)
    ef=effects[effects.tissue.eq(tissue)]
    assert ef.groupby('gene_id').size().eq(6).all()
    assert ef.display_available.all() and np.isfinite(ef.effect_display).all()
    for (cluster,g,day), d in ef.groupby(['cluster','g1','time']):
        profile.append({'tissue':tissue,'cluster':cluster,'genotype':g,'time':day,
                        'mean_response_log2FC':d.effect_display.mean(),'genes':d.gene_id.nunique()})
    for cluster,d in rows.groupby('cluster'):
        for g in genotypes:
            profile.append({'tissue':tissue,'cluster':cluster,'genotype':g,'time':'d0',
                            'mean_response_log2FC':0.,'genes':len(d)})
heat, individual, profile = pd.concat(heat), pd.concat(individual), pd.DataFrame(profile)
save(heat,'B_group_mean_VST_z.csv.gz')
save(individual,'B_individual_VST_z.csv.gz')
save(profile,'B_cluster_response_profiles.csv')
save(order,'B_gene_order_and_cluster.csv.gz')

# Representatives: closest to the response-feature centroid in each cluster,
# among genes passing the archived nine-cell expression criterion and baseMean >=20.
selections=[]
feature_cols=[c for c in features.columns if c not in ['tissue','gene_id']]
for tissue in ['Spleen','Intestine']:
    oo=order[order.tissue.eq(tissue)]
    largest=int(oo.groupby('cluster').size().sort_values(ascending=False,kind='stable').index[0])
    for cluster,rows in oo.groupby('cluster'):
        f=features[features.tissue.eq(tissue)].set_index('gene_id').loc[rows.gene_id,feature_cols]
        centroid=f.mean(axis=0).to_numpy()
        dist=np.sum((f.to_numpy()-centroid)**2,axis=1)
        cand=response.loc[tissue].loc[f.index].copy()
        cand['centroid_distance_squared']=dist
        cand=cand[cand.expressed_all_groups & (cand.baseMean>=20)]
        cand=cand.reset_index().sort_values(['centroid_distance_squared','q_within','gene_id'])
        assert len(cand)>0
        s=cand.iloc[0].to_dict()
        s.update(tissue=tissue,cluster=int(cluster),role='Response representative',
                 main_panel=(int(cluster)==largest))
        selections.append(s)
    s=response.loc[(tissue,'127522978')].to_dict()
    s.update(tissue=tissue,gene_id='127522978',cluster=None,role='Prespecified IL10',main_panel=True)
    selections.append(s)
selected=pd.DataFrame(selections)
selected=selected.merge(baseline[['q_within']].reset_index().rename(columns={'q_within':'baseline_q'}),on=['tissue','gene_id'])
selected=selected.rename(columns={'q_within':'response_q'})
save(selected,'C_gene_selection.csv')
nc=read('normalized_counts_long.csv.gz')
points=nc.merge(selected[['tissue','gene_id','role','cluster','main_panel']],on=['tissue','gene_id'])
points['log2_normalized_count_plus1']=np.log2(points.normalized_count+1)
assert points.groupby(['tissue','gene_id','genotype','time']).size().eq(3).all()
save(points,'C_individual_points.csv')
save(points.groupby(['tissue','gene_id','genotype','time']).log2_normalized_count_plus1.agg(['mean','std','size']).reset_index(),'C_group_summary.csv')

# Omnibus-set ORA; separate hypothesis family for every tissue and test.
# Both categories are included in each corresponding full target list, once.
kegg=read('KEGG_annotation_input.tsv',sep='\t')
sets={str(row.Entry):{s.strip() for s in str(row.Gene_id).split(',') if s.strip()}
      for row in kegg.itertuples() if pd.notna(row.Gene_id)}
annotated=set().union(*sets.values())
names=kegg.set_index('Entry').Name.to_dict()
ora,backgrounds=[],[]
for tissue in ['Spleen','Intestine']:
    for test in ['Baseline','Response']:
        tested=om[(om.tissue==tissue)&(om.test==test)&np.isfinite(om.pvalue)]
        universe=set(tested.gene_id)&annotated
        targets=set(tested[tested.q_within<.05].gene_id)&universe
        records=[]
        for pathway, members in sets.items():
            ss=members&universe
            if not 10<=len(ss)<=500: continue
            hits=ss&targets
            N,M,K,x=len(universe),len(targets),len(ss),len(hits)
            records.append({'tissue':tissue,'test':test,'pathway_id':pathway,'pathway':names[pathway],
                'background_genes':N,'target_genes':M,'pathway_genes':K,'hits':x,
                'fold_enrichment':(x/M)/(K/N), 'pvalue':hypergeom.sf(x-1,N,K,M),
                'hit_gene_ids':';'.join(sorted(hits))})
        df=pd.DataFrame(records)
        df['q_within']=bh(df.pvalue)
        df=df.sort_values(['q_within','pvalue','pathway_id'])
        ora.append(df)
        backgrounds.append({'tissue':tissue,'test':test,'all_detected_genes':int((tested.q_within<.05).sum()),
            'annotated_background':len(universe),'annotated_targets':len(targets),
            'pathways_tested':len(df),'significant_pathways':int((df.q_within<.05).sum())})
ora=pd.concat(ora,ignore_index=True)
save(ora,'D_KEGG_ORA_all.csv')
save(pd.DataFrame(backgrounds),'D_background_and_detection_summary.csv')
shown=ora[ora.q_within<.05]
save(shown,'D_main_displayed_pathways.csv')
old=read('original_response_ORA.csv')
check=ora[ora.test=='Response'].merge(old,on=['tissue','pathway_id'],suffixes=('_new','_old'))
assert len(check)==len(old)==336
for field in ['fold_enrichment','pvalue','q_within']:
    assert np.allclose(check[field+'_new'],check[field+'_old'],rtol=1e-9,atol=1e-15),field
checksums=[]
for f in sorted(inp.iterdir()):
    checksums.append({'file':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size})
save(pd.DataFrame(checksums),'input_checksums.csv')
qa={'sample_count':len(md),'counts':counts.to_dict('records'),
    'clusters':order.groupby(['tissue','cluster']).size().rename('genes').reset_index().to_dict('records'),
    'main_selected_genes':json.loads(selected[selected.main_panel][['tissue','gene_id','role','cluster','baseline_q','response_q']].to_json(orient='records',double_precision=15)),
    'ORA_summary':backgrounds,'original_response_ORA_reproduced':True,
    'new_DESeq2_or_limma_fits':False,'extra_fold_change_cutoff':None,
    'heatmap':'Per-gene VST z over 27 samples; mean across 3 fish for each of 9 groups.',
    'profiles':'Mean archived own-baseline model log2FC across cluster genes; no inferential error bars.',
    'ORA':'Unoriented omnibus target sets; separate BH family per tissue and test.',
    'versions':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}}
(out/'numerical_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'main_genes':qa['main_selected_genes'],'clusters':qa['clusters'],'ORA':backgrounds},ensure_ascii=True,indent=2))
