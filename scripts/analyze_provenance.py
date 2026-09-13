from pathlib import Path
import anndata as ad, numpy as np,pandas as pd,scipy.stats as st,json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'reports';S=R/'provenance_sources'
a=ad.read_h5ad(ROOT/'data/norman/perturb_processed.h5ad');o=a.obs.copy();o=o.astype({'condition':str,'cell_type':str,'dose_val':str,'condition_name':str});o['barcode_group']=o.index.str.rsplit('-',n=1).str[-1].astype(int)
def canon(s):return '+'.join(sorted(set(g for g in s.split('+') if g!='ctrl'))) or 'ctrl'
o['canonical_perturbation']=o.condition.map(canon);o['n_targets']=o.canonical_perturbation.map(lambda s:0 if s=='ctrl' else len(s.split('+')))
s=pd.read_csv(S/'GSE133344_raw_cell_identities.csv.gz').set_index('cell_barcode');filtered=pd.read_csv(S/'GEO_cell_identities_alternate.csv.gz').set_index('cell_barcode');assert s.index.is_unique
j=o.join(s,validate='one_to_one');assert j.guide_identity.notna().all()
def guide_condition(x):
 tokens=x.split('__')[1].split('_')[:2]
 tokens=['ctrl' if t.startswith('NegCtrl') else t for t in tokens]
 return 'ctrl' if tokens==['ctrl','ctrl'] else '+'.join(tokens)
j['source_condition_unmapped']=j.guide_identity.map(guide_condition)
j['source_condition']=j.source_condition_unmapped.map(lambda x:'+'.join('RHOXF2BB' if g=='RHOXF2' else g for g in x.split('+')));assert (j.condition==j.source_condition).all();assert (j.barcode_group==j.gemgroup).all();assert j.good_coverage.all() and (j.number_of_cells==1).all()
j.to_csv(R/'verified_cell_metadata.csv')
guide=j.groupby(['condition','canonical_perturbation','guide_identity'],observed=True).size().rename('n_cells').reset_index();guide.to_csv(R/'guide_condition_mapping.csv',index=False)
fields={c:{'nunique':int(j[c].nunique()),'missing':int(j[c].isna().sum()),'examples':j[c].dropna().unique()[:6].tolist()} for c in s.columns}
# Descriptive and non-parametric checks only. No model fitting.
ct=pd.crosstab(o.condition,o.barcode_group);chi,pval,df,expected=st.chi2_contingency(ct)
controlct=pd.crosstab(o.control,o.barcode_group);cc,cp,cd,ce=st.chi2_contingency(controlct)
counts=np.asarray(a.layers['counts'].sum(axis=1,dtype=np.float64));det=a.X.getnnz(axis=1);o['retained_gene_counts']=counts;o['detected_features']=det
rows=[];gene_rows=[];gene_set=sorted({g for v in o.canonical_perturbation if v!='ctrl' for g in v.split('+')})
for group,b in o.groupby('barcode_group'):
 control=b.control==1;pert=b.loc[~control];idx=o.index.get_indexer(b.index)
 er={'barcode_group':int(group),'expression_sample':f'GSM{3906019+group}','guide_sample':f'GSM{3906027+group}','cells':len(b),'controls':int(control.sum()),'control_fraction':float(control.mean()),'noncontrol_cells':int((~control).sum()),'original_conditions_including_control':int(b.condition.nunique()),'target_sets_including_control':int(b.canonical_perturbation.nunique()),'single_target_cells':int((b.n_targets==1).sum()),'double_target_cells':int((b.n_targets==2).sum()),'unique_target_genes':len({g for s0 in pert.canonical_perturbation for g in s0.split('+')}),'median_counts_retained_genes':float(np.median(counts[idx])),'median_detected_genes':float(np.median(det[idx])),'control_median_counts_retained_genes':float(b.loc[control,'retained_gene_counts'].median()),'control_median_detected_genes':float(b.loc[control,'detected_features'].median())}
 rows.append(er)
 for gene in gene_set:
  target=b.canonical_perturbation.map(lambda x:gene in x.split('+'))
  gene_rows.append({'barcode_group':int(group),'target_gene':gene,'n_cells':int(target.sum()),'single_target_cells':int((target&(b.n_targets==1)).sum()),'double_target_cells':int((target&(b.n_targets==2)).sum()),'fraction_of_all_group_cells':float(target.mean()),'fraction_of_noncontrol_group_cells':float(target.sum()/len(pert))})
summary=pd.DataFrame(rows);summary.to_csv(R/'barcode_group_summary.csv',index=False);tg=pd.DataFrame(gene_rows);tg.to_csv(R/'target_gene_by_barcode_group.csv',index=False)
tg.pivot(index='target_gene',columns='barcode_group',values='n_cells').to_csv(R/'target_gene_by_barcode_group_wide.csv')
# Control-only summaries avoid perturbation-composition confounding; no causal batch-effect claim.
ctrl=o.loc[o.control==1];kw={}
for col in ['retained_gene_counts','detected_features']:
 k=st.kruskal(*[b[col].to_numpy() for _,b in ctrl.groupby('barcode_group')]);kw[col]={'H':float(k.statistic),'p_value':float(k.pvalue),'epsilon_squared':float(max(0,(k.statistic-7)/(len(ctrl)-8)))}
means=np.stack([np.asarray(a.X[(o.control.to_numpy()==1)&(o.barcode_group.to_numpy()==g)].mean(axis=0)).ravel() for g in range(1,9)])
corr=np.corrcoef(means);pd.DataFrame(corr,index=range(1,9),columns=range(1,9)).to_csv(R/'control_profile_correlations.csv')
# Source disagreement is retained, not silently harmonized.
fj=o.join(filtered,rsuffix='_filtered');present=fj.guide_identity.notna();diff=present&(fj.guide_identity!=j.guide_identity)
fj.loc[diff,['condition','guide_identity','good_coverage','number_of_cells']].assign(raw_guide_identity=j.loc[diff,'guide_identity']).to_csv(R/'filtered_source_conflicts.csv')
# Check exact sparse duplicate expression rows: never use arbitrary synthetic paired cells as truth.
seen={};dups=0
for i in range(a.n_obs):
 lo,hi=a.X.indptr[i:i+2];sig=hashlib.sha256(a.X.indices[lo:hi].tobytes()+a.X.data[lo:hi].tobytes()).digest()
 if sig in seen:dups+=1
 else:seen[sig]=i
# Feasible, frozen proposed split assignments. No models are trained.
seed=20260911
def rank(values,tag):return sorted(values,key=lambda x:hashlib.sha256(f'{seed}|{tag}|{x}'.encode()).hexdigest())
manifest=o[['condition','canonical_perturbation','barcode_group','n_targets']].copy()
manifest['A_random_cell']=''
for condition,b in o.groupby('condition'):
 ids=rank(list(b.index),'A');n=len(ids);nt=max(1,int(np.floor(.1*n)));nv=max(1,int(np.floor(.1*n)));assert n>nt+nv
 manifest.loc[ids[:nt],'A_random_cell']='test';manifest.loc[ids[nt:nt+nv],'A_random_cell']='validation';manifest.loc[ids[nt+nv:],'A_random_cell']='train'
pairs=rank(sorted(o.loc[o.n_targets==2,'canonical_perturbation'].unique()),'B_pairs');assert len(pairs)==131
pairmap={p0:('test' if k<26 else 'validation' if k<52 else 'train') for k,p0 in enumerate(pairs)}
manifest['B_heldout_combination']=o.canonical_perturbation.map(pairmap).fillna('train')
for group,b in o.loc[o.control==1].groupby('barcode_group'):
 ids=rank(list(b.index),'B_control');n=len(ids);nt=int(np.floor(.1*n));nv=nt
 manifest.loc[ids[:nt],'B_heldout_combination']='test';manifest.loc[ids[nt:nt+nv],'B_heldout_combination']='validation'
genes=rank(gene_set,'C_genes');testgenes=set(genes[:16]);valgenes=set(genes[16:32]);traingenes=set(genes[32:])
def cs(x):
 g=set(x.split('+'))-{'ctrl'}
 if g&testgenes and g&valgenes:return 'excluded_val_test_bridge'
 if g&testgenes:return 'test'
 if g&valgenes:return 'validation'
 return 'train'
manifest['C_heldout_target_gene']=o.canonical_perturbation.map(cs)
ctrlidx=o.index[o.control==1];manifest.loc[ctrlidx,'C_heldout_target_gene']=manifest.loc[ctrlidx,'B_heldout_combination']
manifest['D_heldout_barcode_group']=o.barcode_group.map(lambda g:'train' if g<=6 else 'validation' if g==7 else 'test')
manifest.to_csv(R/'proposed_split_assignments.csv',index_label='cell_barcode')
splits={}
for col in ['A_random_cell','B_heldout_combination','C_heldout_target_gene','D_heldout_barcode_group']:
 splits[col]={}
 for split,b in manifest.groupby(col):
  nonctrl=b.loc[b.n_targets>0];gg={g for x in nonctrl.canonical_perturbation for g in x.split('+')}
  splits[col][split]={'cells':len(b),'controls':int((b.n_targets==0).sum()),'noncontrol_cells':len(nonctrl),'noncontrol_target_sets':int(nonctrl.canonical_perturbation.nunique()),'single_target_sets':int(b.loc[b.n_targets==1,'canonical_perturbation'].nunique()),'double_target_sets':int(b.loc[b.n_targets==2,'canonical_perturbation'].nunique()),'unique_targets_in_conditions':len(gg)}
# Assertions for leakage prevention.
for key in ['B_heldout_combination','C_heldout_target_gene']:
 sets=[set(manifest.loc[(manifest[key]==s0)&(manifest.n_targets>0),'canonical_perturbation']) for s0 in ['train','validation','test']]
 assert all(not (sets[i]&sets[j0]) for i in range(3) for j0 in range(i))
for split,forbidden in [('train',testgenes|valgenes),('validation',testgenes),('test',valgenes)]:
 used={g for c0 in manifest.loc[(manifest.C_heldout_target_gene==split)&(manifest.n_targets>0),'canonical_perturbation'] for g in c0.split('+')};assert not used&forbidden
assert manifest.index.is_unique and manifest.notna().all().all()
geo_excluded='NegCtrl1_NegCtrl0__NegCtrl1_NegCtrl0'
source_good=s.loc[(s.good_coverage==True)&(s.number_of_cells==1)]
r={'dataset_sha256':hashlib.file_digest((ROOT/'data/norman/perturb_processed.h5ad').open('rb'),'sha256').hexdigest(),'source_raw_rows':len(s),'source_raw_barcode_matches':int(o.index.isin(s.index).sum()),'source_condition_matches_before_documented_RHOXF2_mapping':int((j.condition==j.source_condition_unmapped).sum()),'source_condition_exact_matches':int((j.condition==j.source_condition).sum()),'gemgroup_exact_matches':int((j.barcode_group==j.gemgroup).sum()),'source_fields':fields,'source_good_coverage_counts':j.good_coverage.value_counts().to_dict(),'source_number_of_cells_counts':j.number_of_cells.value_counts().to_dict(),'source_cellranger_called_counts':j.cellranger_called.value_counts().to_dict(),'source_filtered_matches':int(present.sum()),'source_filtered_missing':int((~present).sum()),'source_filtered_guide_conflicts':int(diff.sum()),'source_good_singlet_rows':len(source_good),'source_excluded_control_good_singlets':int((source_good.guide_identity==geo_excluded).sum()),'raw_good_singlet_minus_excluded_equals_local_barcodes':set(source_good.loc[source_good.guide_identity!=geo_excluded].index)==set(o.index),'guide_identities':int(j.guide_identity.nunique()),'guide_variants_per_condition':guide.groupby('condition').size().value_counts().to_dict(),'control_guide_counts':j.loc[j.control==1].guide_identity.value_counts().to_dict(),'group_condition_chi_square':{'statistic':float(chi),'df':int(df),'p_value':float(pval),'cramers_v':float(np.sqrt(chi/(len(o)*7))),'min_expected':float(expected.min())},'group_control_chi_square':{'statistic':float(cc),'df':int(cd),'p_value':float(cp),'cramers_v':float(np.sqrt(cc/len(o)))},'control_only_kruskal':kw,'control_mean_profile_pairwise_correlation_min':float(corr[np.triu_indices(8,1)].min()),'control_mean_profile_pairwise_correlation_max':float(corr[np.triu_indices(8,1)].max()),'exact_duplicate_X_rows':dups,'target_genes':gene_set,'split_seed':seed,'splits':splits,'B_pair_assignments':pairmap,'C_train_genes':sorted(traingenes),'C_validation_genes':sorted(valgenes),'C_test_genes':sorted(testgenes),'group_summary':rows}
(R/'provenance_evidence.json').write_text(json.dumps(r,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf8')
print(json.dumps({k:v for k,v in r.items() if k not in ['source_fields','B_pair_assignments','target_genes','group_summary']},indent=2,default=str))
