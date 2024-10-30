import pandas as pd

x=pd.read_csv('sqantifiltered_monoexonicfiltered_2reads2samples_classification.txt', sep='\t')

x=x[x['chrom'].isna()==False]

for i in x['chrom'].unique():
    x.loc[x['associated_gene'].str.contains('novelGene') & (x['chrom']==i),'associated_gene']=x[x['associated_gene'].str.contains('novelGene') & (x['chrom']==i)]['associated_gene'].str.replace('novelGene_', 'novelGene_'+i+'_')

x.to_csv('sqantifiltered_monoexonicfiltered_2reads2samples_classification_finalversion.txt', sep='\t', index=False, na_rep='NA')
