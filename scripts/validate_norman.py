from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,json
root=Path(__file__).resolve().parents[1];p=root/'reports';a=ad.read_h5ad(root/'data/norman/perturb_processed.h5ad');x=a.X;c=a.layers['counts'];n=a.n_obs
factors=np.empty(n);worst=0.;bad=0;pattern_bad=0;examples=[]
for start in range(0,n,512):
 end=min(n,start+512);b=x[start:end];rows=np.repeat(np.arange(end-start),np.diff(b.indptr));cv=c[start:end][rows,b.indices].astype('float64');vals=b.data.astype('float64')
 pattern_bad+=int(np.count_nonzero(cv<=0));ratios=np.expm1(vals)/cv
 first=b.indptr[:-1];scales=ratios[first];factors[start:end]=scales
 err=np.abs(vals-np.log1p(cv*scales[rows]));worst=max(worst,float(err.max()));bad+=int(np.count_nonzero(err>2e-6))
 if start==0:
  for j in range(12):examples.append({'barcode':a.obs_names[rows[j]],'gene_id':a.var_names[b.indices[j]],'gene_name':str(a.var.iloc[b.indices[j]].gene_name),'count':float(cv[j]),'X':float(vals[j]),'expm1_X_over_count':float(ratios[j])})
pattern_bad+=int(np.count_nonzero(c))-x.nnz
suffix=a.obs_names.str.rsplit('-',n=1).str[-1];ct=pd.crosstab(a.obs.condition,suffix);ct.to_csv(p/'condition_by_barcode_suffix.csv')
bs=pd.DataFrame({'cells':pd.Series(suffix).value_counts().sort_index(),'control_cells':pd.Series(suffix[a.obs.control.to_numpy()==1]).value_counts().sort_index(),'condition_labels':(ct>0).sum()});bs.to_csv(p/'barcode_suffix_counts.csv',index_label='barcode_suffix')
uns={}
for k,d in a.uns.items():
 arrays=list(d.values());keys=set(d);cond=set(a.obs.condition_name.astype(str));first_key=next(iter(d));vs=np.concatenate(arrays)
 numeric=np.issubdtype(vs.dtype,np.number)
 uns[k]={'n_keys':len(d),'array_length_min':min(map(len,arrays)),'array_length_max':max(map(len,arrays)),'keys_not_in_condition_name':sorted(keys-cond),'conditions_without_key':sorted(cond-keys),'example_key':first_key,'example_values':d[first_key][:10].tolist(),'value_type':str(vs.dtype),'values_in_var_names':None if numeric else set(vs).issubset(set(a.var_names)),'min_value':int(vs.min()) if numeric else None,'max_value':int(vs.max()) if numeric else None}
sums=c.sum(axis=1,dtype=np.float64);normalized_totals=sums*factors
r={'checked_nonzero_entries':int(x.nnz),'zero_pattern_mismatches':pattern_bad,'max_abs_log1p_scaled_count_error':worst,'entries_log1p_error_above_2e-6':bad,'scale_quantiles':np.quantile(factors,[0,.25,.5,.75,1]).tolist(),'count_row_sum_times_scale_quantiles':np.quantile(normalized_totals,[0,.25,.5,.75,1]).tolist(),'normalization_10000_hypothesis':{'implied_original_total_quantiles':np.quantile(10000/factors,[0,.25,.5,.75,1]).tolist(),'max_distance_to_nearest_integer':float(np.abs(10000/factors-np.rint(10000/factors)).max()),'fraction_within_0_01_of_integer':float((np.abs(10000/factors-np.rint(10000/factors))<.01).mean())},'api_layers_None_alias_equals_X':bool((x!=a.layers[None]).nnz==0),'empty_expression_rows':int((np.diff(x.indptr)==0).sum()),'empty_expression_columns':int((x.getnnz(axis=0)==0).sum()),'examples':examples,'barcode_suffix':bs.to_dict(orient='index'),'barcode_without_suffix_duplicates':int(a.obs_names.str.rsplit('-',n=1).str[0].duplicated().sum()),'uns':uns}
(p/'validation_evidence.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
