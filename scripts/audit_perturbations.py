from pathlib import Path
import pandas as pd, json
p=Path(__file__).resolve().parents[1]/'reports';o=pd.read_csv(p/'obs_metadata.csv',index_col=0);v=pd.read_csv(p/'var_metadata.csv',index_col=0)
def targets(s):return sorted(set(t for t in s.split('+') if t!='ctrl'))
counts=o.groupby('condition',observed=True).size().rename('n_cells').reset_index()
counts['target_genes']=counts.condition.map(lambda s:'+'.join(targets(s)))
counts['n_targets']=counts.condition.map(lambda s:len(targets(s)))
counts['canonical_perturbation']=counts.target_genes.replace('', 'ctrl')
counts.to_csv(p/'perturbation_counts.csv',index=False)
canon=counts.groupby(['canonical_perturbation','n_targets'],as_index=False).agg(n_cells=('n_cells','sum'),n_original_labels=('condition','size'),original_labels=('condition',lambda x:';'.join(sorted(x))))
canon.to_csv(p/'canonical_perturbation_counts.csv',index=False)
genes=sorted(set(g for s in counts.condition for g in targets(s)))
rows=[]
for g in genes:
 hit=counts.condition.map(lambda s:g in targets(s));single=hit & (counts.n_targets==1);pair=hit & (counts.n_targets==2)
 rows.append({'target_gene':g,'gene_id':';'.join(v.index[v.gene_name==g]),'measured_in_X':bool((v.gene_name==g).any()),'cells_single':int(counts.loc[single,'n_cells'].sum()),'cells_combinatorial':int(counts.loc[pair,'n_cells'].sum()),'cells_any':int(counts.loc[hit,'n_cells'].sum()),'n_original_labels':int(hit.sum())})
pd.DataFrame(rows).to_csv(p/'target_genes.csv',index=False)
r={'condition_label_counts_by_n_targets':counts.groupby('n_targets').size().to_dict(),'cell_counts_by_n_targets':counts.groupby('n_targets').n_cells.sum().to_dict(),'canonical_counts_by_n_targets':canon.groupby('n_targets').size().to_dict(),'n_target_genes':len(genes),'target_genes':genes,'targets_missing_in_var':[g for g in genes if g not in set(v.gene_name)],'control_agrees':bool(((o.condition=='ctrl')==(o.control==1)).all()),'condition_name_agrees':bool((o.condition_name==(o.cell_type+'_'+o.condition+'_'+o.dose_val)).all()),'dose_agrees':bool((o.dose_val==o.condition.map(lambda s:'1' if s=='ctrl' else '1+1')).all()),'barcode_suffix_counts':o.index.str.rsplit('-',n=1).str[-1].value_counts().to_dict(),'barcode_length_counts':o.index.str.len().value_counts().to_dict(),'count_summary_all':counts.n_cells.describe().to_dict(),'count_summary_noncontrol':counts.loc[counts.n_targets>0,'n_cells'].describe().to_dict(),'canonical_count_summary_noncontrol':canon.loc[canon.n_targets>0,'n_cells'].describe().to_dict(),'merged_canonical_groups':canon.loc[canon.n_original_labels>1].to_dict('records'),'lowest_counts':counts.nsmallest(10,'n_cells').to_dict('records'),'highest_counts':counts.nlargest(10,'n_cells').to_dict('records')}
(p/'perturbation_evidence.json').write_text(json.dumps(r,indent=2))
print(json.dumps(r,indent=2))
