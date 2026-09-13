from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp, pandas as pd, json, hashlib, importlib.metadata, h5py
root=Path(__file__).resolve().parents[1]; p=root/'data/norman/perturb_processed.h5ad'; report=root/'reports';report.mkdir(exist_ok=True)
a=ad.read_h5ad(p)
print(a,flush=True)
def plain(x):
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,np.ndarray):return x.tolist()
    return str(x)
def profile(df):
    o={}
    for c in df:
        s=df[c];d={'dtype':str(s.dtype),'missing':int(s.isna().sum()),'nunique':int(s.nunique()),'examples':[plain(x) for x in s.dropna().unique()[:12]]}
        if s.nunique()<=500:d['value_counts']={str(k):int(v) for k,v in s.value_counts(dropna=False).items()}
        if pd.api.types.is_numeric_dtype(s.dtype):d['describe']=s.describe().to_dict()
        o[c]=d
    return o
def matrix(m):
    v=m.data if sp.issparse(m) else np.asarray(m).ravel()
    sums=np.asarray(m.sum(axis=1)).ravel()
    return {'shape':list(m.shape),'type':type(m).__name__,'dtype':str(m.dtype),'stored_entries':int(v.size),'nonzero':int(m.count_nonzero() if sp.issparse(m) else np.count_nonzero(m)),'min_stored':float(v.min()),'max_stored':float(v.max()),'negative_entries':int(np.count_nonzero(v<0)),'nonfinite_entries':int(np.count_nonzero(~np.isfinite(v))),'fractional_entries_tolerance_1e-6':int(np.count_nonzero(np.abs(v-np.rint(v))>1e-6)),'value_quantiles':np.quantile(v,[0,.25,.5,.75,.99,1]).tolist(),'row_sum_quantiles':np.quantile(sums,[0,.25,.5,.75,1]).tolist(),'first_20_stored_values':v[:20].tolist()}
r={'file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'versions':{k:importlib.metadata.version(k) for k in ['anndata','numpy','scipy','pandas','h5py']},'shape':list(a.shape),'obs_index':{'name':a.obs.index.name,'unique':bool(a.obs_names.is_unique),'examples':a.obs_names[:12].tolist()},'var_index':{'name':a.var.index.name,'unique':bool(a.var_names.is_unique),'examples':a.var_names[:12].tolist()},'obs':profile(a.obs),'var':profile(a.var),'X':matrix(a.X),'layers':{k:matrix(v) for k,v in a.layers.items() if k is not None},'raw':None if a.raw is None else {'matrix':matrix(a.raw.X),'var':profile(a.raw.var)},'uns':{k:{'type':type(v).__name__,'shape':getattr(v,'shape',None),'repr':repr(v)[:3000]} for k,v in a.uns.items()},'obsm':{k:list(v.shape) for k,v in a.obsm.items()},'varm':{k:list(v.shape) for k,v in a.varm.items()},'obsp':{k:list(v.shape) for k,v in a.obsp.items()}}
if sp.issparse(a.X):
    inv=a.X.astype(np.float64).copy();inv.data=np.expm1(inv.data)
    r['X']['expm1_row_sum_quantiles']=np.quantile(np.asarray(inv.sum(axis=1)).ravel(),[0,.25,.5,.75,1]).tolist()
    r['X']['expm1_fractional_entries_tolerance_1e-4']=int(np.count_nonzero(np.abs(inv.data-np.rint(inv.data))>1e-4))
a.obs.to_csv(report/'obs_metadata.csv');a.var.to_csv(report/'var_metadata.csv')
with h5py.File(p, 'r') as h:
    r['hdf5_root_keys']=list(h.keys())
    r['hdf5_stored_layer_keys']=list(h['layers'].keys())
r['api_note']='AnnData 0.13 exposes layers[None] as an alias for X; this is not an additional stored layer.'
(report/'audit_evidence.json').write_text(json.dumps(r,indent=2,default=plain),encoding='utf-8')
print(json.dumps(r,indent=2,default=plain),flush=True)
