"""Describe all heterozygote contrasts without an outside-both filter."""
from pathlib import Path
import sys, argparse, shutil, json, hashlib
from urllib.parse import unquote
sys.path.insert(0,str(Path(__file__).resolve().parent/'python_packages'))
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();out=args.output;inp=out/'inputs';tab=out/'source_tables'
inp.mkdir(parents=True,exist_ok=True);tab.mkdir(exist_ok=True)
files={'DESeq2_pairwise.csv.gz':'Figure4_5_20261009/DESeq2_pairwise.csv.gz',
       'normalized_counts_long.csv.gz':'Figure4_5_20261009/normalized_counts_long.csv.gz',
       'sample_metadata.csv':'Figure4_5_20261009/sample_metadata.csv',
       'contrast_definitions.csv':'Figure4_5_20261009/DESeq2_contrast_definitions.csv',
       'omnibus_tests.csv.gz':'Figure2_20261008/omnibus_tests.csv.gz'}
if args.source:
    for name,rel in files.items():shutil.copy2(args.source/rel,inp/name)
def read(name):return pd.read_csv(inp/name,dtype={'gene_id':str})
def save(df,name):df.to_csv(tab/name,index=False,encoding='utf-8-sig')

d=read('DESeq2_pairwise.csv.gz')
assert not d.duplicated(['tissue','stage','gene_id','comparator']).any()
assert set(d.comparator)=={'HomoDel','HomoIns'}
assert set(d.stage)=={'Baseline','Response_d1','Response_d5'}
# Check archived global BH families without redefining them for selected genes.
qcheck={}
for family,z in d.groupby('family'):
    valid=np.isfinite(z.pvalue)
    q=false_discovery_control(z.loc[valid,'pvalue'].to_numpy(),method='bh')
    assert np.allclose(q,z.loc[valid,'q_pair_family'],rtol=1e-10,atol=1e-14)
    qcheck[family]={'rows_including_missing':len(z),'valid_pvalues':int(valid.sum()),
                    'archived_global_BH_matches':True}
save(d,'all_heterozygote_contrasts.csv.gz')
keys=['tissue','stage','gene_id']
fields=['effect','lfcSE','pvalue','q_pair_family','ci95_low','ci95_high','zero_count_convention']
left=d[d.comparator.eq('HomoDel')].copy()
right=d[d.comparator.eq('HomoIns')].copy()
w=left[keys+fields+['strict','baseMean','gene_name','gene_description']].rename(columns={c:c+'_D' for c in fields})
w=w.merge(right[keys+fields].rename(columns={c:c+'_I' for c in fields}),on=keys,validate='one_to_one')
assert len(w)==134403
w['fully_testable']=np.isfinite(w.q_pair_family_D)&np.isfinite(w.q_pair_family_I)
w['detected_D']=w.q_pair_family_D<.05;w['detected_I']=w.q_pair_family_I<.05
w['pattern']='Neither_detected'
w.loc[w.detected_D&~w.detected_I,'pattern']='Del_only'
w.loc[~w.detected_D&w.detected_I,'pattern']='Ins_only'
both=w.detected_D&w.detected_I
assert ((w.loc[both,'effect_D']*w.loc[both,'effect_I'])!=0).all()
w.loc[both&(w.effect_D*w.effect_I>0),'pattern']='Both_same'
w.loc[both&(w.effect_D*w.effect_I<0),'pattern']='Both_opposite'
w.loc[~w.fully_testable,'pattern']='Incomplete'
w['direction_D']=np.sign(w.effect_D);w['direction_I']=np.sign(w.effect_I)
w['zero_count_convention_any']=w.zero_count_convention_D|w.zero_count_convention_I
w['primary_annotation']=w.gene_description.fillna('').str.split(' && ').str[0].map(unquote)
w['informative_primary_annotation']=~w.primary_annotation.str.contains('uncharacterized|hypothetical|unknown|^$',case=False,regex=True)
save(w,'A_all_genes_paired_comparisons.csv.gz')
summary=w.groupby(['tissue','stage','pattern'],sort=False).size().rename('genes').reset_index()
save(summary,'B_pattern_counts.csv')
totals=[]
for (tissue,stage),z in w.groupby(['tissue','stage']):
    any_detected=z.fully_testable&(z.detected_D|z.detected_I)
    totals.append({'tissue':tissue,'stage':stage,'all_genes':len(z),
        'fully_testable':int(z.fully_testable.sum()),'incomplete_tests':int((~z.fully_testable).sum()),
        'detected_Hetero_vs_Del':int(z.detected_D.sum()),'detected_Hetero_vs_Ins':int(z.detected_I.sum()),
        'union_detected':int(any_detected.sum()),'both_detected':int((z.detected_D&z.detected_I).sum()),
        'neither_detected':int((z.pattern=='Neither_detected').sum()),
        'zero_count_convention_records':int(z.zero_count_convention_any.sum())})
totals=pd.DataFrame(totals);save(totals,'B_totals_and_denominators.csv')

# Representatives illustrate evidence patterns, not presumed protection functions.
# Prioritize opposite-sign candidates first to avoid reusing their few gene IDs.
patterns=['Both_opposite','Both_same','Del_only','Ins_only']
selected=[]
for tissue in ['Spleen','Intestine']:
    used=set()
    for pattern in patterns:
        candidates=w[(w.tissue==tissue)&w.stage.ne('Baseline')&(w.pattern==pattern)&w.strict&(w.baseMean>=20)&~w.zero_count_convention_any&~w.gene_id.isin(used)].copy()
        assert len(candidates)>0,(tissue,pattern)
        stage_counts=w[(w.tissue==tissue)&w.stage.ne('Baseline')&(w.pattern==pattern)].groupby('stage').size()
        stages=sorted(stage_counts.index,key=lambda s:(-int(stage_counts[s]),s))
        stage=next(s for s in stages if (candidates.stage==s).any())
        candidates=candidates[candidates.stage==stage].copy()
        center=candidates[['effect_D','effect_I']].median().to_numpy()
        scale=np.median(np.abs(candidates[['effect_D','effect_I']].to_numpy()-center),axis=0)
        scale=np.maximum(scale,.25)
        candidates['distance_to_pattern_center']=np.sum(((candidates[['effect_D','effect_I']].to_numpy()-center)/scale)**2,axis=1)
        # An informative source annotation is a display preference, not validation.
        fallback=not candidates.informative_primary_annotation.any()
        pool=candidates if fallback else candidates[candidates.informative_primary_annotation]
        pool=pool.sort_values(['distance_to_pattern_center','gene_id'])
        s=pool.iloc[0].to_dict();s['annotation_display_fallback']=fallback
        s['selection_rule']='Largest response-stage pattern count with eligible candidates; nearest robust pattern center; informative primary annotation preferred.'
        selected.append(s);used.add(s['gene_id'])
selected=pd.DataFrame(selected)
save(selected,'C_selected_response_examples.csv')
nc=read('normalized_counts_long.csv.gz')
pts=nc.merge(selected[['tissue','gene_id','stage','pattern']].rename(columns={'stage':'selected_stage','pattern':'selected_pattern'}),on=['tissue','gene_id'],validate='many_to_one')
pts['log2_normalized_count_plus1']=np.log2(pts.normalized_count+1)
assert pts.groupby(['tissue','gene_id','genotype','time']).size().eq(3).all()
save(pts,'C_individual_expression_points.csv')
save(pts.groupby(['tissue','gene_id','genotype','time']).log2_normalized_count_plus1.agg(['mean','std','size']).reset_index(),'C_expression_summary.csv')
ef=d.merge(selected[['tissue','stage','gene_id','pattern']],on=['tissue','stage','gene_id'],validate='many_to_one')
save(ef,'D_selected_effects_and_intervals.csv')
selected_all=d.merge(selected[['tissue','gene_id','pattern']],on=['tissue','gene_id'],validate='many_to_one')
save(selected_all,'supplement_selected_examples_all_stages.csv')
om=read('omnibus_tests.csv.gz')
save(om.merge(selected[['tissue','gene_id','pattern']],on=['tissue','gene_id']),'supplement_selected_examples_omnibus.csv')
hashes=[]
for f in sorted(inp.iterdir()):hashes.append({'file':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
save(pd.DataFrame(hashes),'input_checksums.csv')
qa={'model_refitted':False,'extra_fold_change_cutoff':None,'IUT_used_as_primary_filter':False,
    'omnibus_used_as_additional_gate':False,'BH_families':qcheck,
    'selection_count':len(selected),'totals':totals.to_dict('records'),
    'category_note':'Both_same/opposite describe two individually detected comparisons and estimated signs; neither is a newly calibrated joint test.',
    'single_note':'A non-detected comparison is not evidence of equality.',
    'selected_examples':json.loads(selected[['tissue','stage','gene_id','pattern','primary_annotation','q_pair_family_D','q_pair_family_I']].to_json(orient='records',double_precision=15)),
    'versions':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}}
(out/'numerical_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(totals.to_string(index=False))
print(selected[['tissue','stage','pattern','gene_id','primary_annotation','q_pair_family_D','q_pair_family_I']].to_string(index=False))
