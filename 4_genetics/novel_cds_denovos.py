import pandas as pd
import pyranges as pr
from tqdm import tqdm
import numpy as np
#import pyBigWig
#
#def get_phylop(ncds,bw):
#    all_vals=np.array([])
#    for i in tqdm(range(ncds.shape[0])):
#        tmp=ncds.iloc[i]
#        vals=np.array(bw.values(tmp['Chromosome'], tmp['Start'], tmp['End']+1))
#        all_vals=np.append(all_vals, vals)
#    return(all_vals)
#
#bw = pyBigWig.open("/lustre/home/vc362/MRC148213/vc362/conservation/phylop/hg38.cactus241way.phyloP.bw")
novel_cds=pr.read_bed('../conservation/gencode_defined_novel_cds.bed')
novel_cds=novel_cds.drop_duplicate_positions()
gencode_introns=pd.read_csv('gencode_introns.txt')
#gencode_introns=gencode_introns[['seqnames', 'start', 'end', 'gene_id', 'transcript_id', 'strand']]
#gencode_introns.columns=['Chromosome', 'Start', 'End', 'gene_id', 'transcript_id', 'Strand']
gencode_introns=pr.PyRanges(gencode_introns)
#gencode=pr.read_gtf('/lustre/home/vc362/resources/gencode.v43.annotation.gtf')
#gencode_cds=pr.read_gff3('../conservation/gencode_cds.gff')
gencode=pr.read_gtf('/lustre/home/vc362/resources/gencode.v44.annotation.gtf')
gencode_cds=gencode[gencode.Feature=='CDS']
novel_cds=novel_cds.intersect(gencode_cds, invert=True)
novel_cds=novel_cds.drop_duplicate_positions()
gencode_exons=gencode[(gencode.Feature=='exon')]
#gencode_utr=gencode_exons.intersect(gencode_cds, invert=True)
#gencode_utr=gencode_utr.intersect(novel_cds, invert=True)
gencode_utr=gencode[gencode.Feature=='UTR']
gencode_introns=gencode_introns.intersect(novel_cds, invert=True)
gencode_introns=gencode_introns.drop_duplicate_positions()
gencode_utr=gencode_utr.drop_duplicate_positions()
gencode_cds=gencode_cds.drop_duplicate_positions()

dn=pd.read_csv('denovo_db_hg38_phylop.tsv')
dn=dn[dn.StudyName.isin(['Turner_2017'])]
dn=pr.PyRanges(dn)

dnm=dn
dnm=dn[dn.Gene.isin(sfari['gene-symbol'])]
#dnm=dn[dn.Gene.isin(loeuf[loeuf['pLI']>0.9]['gene'])]
dnm=dn[dn.Gene.isin(sfari[sfari['gene-score']<3.0]['gene-symbol'])]
dnm=dnm[['Chromosome', 'Start', 'End', 'PrimaryPhenotype', 'SampleID', 'CaddScore', 'phylop']]


count={'introns':[], 'cds':[], 'utr':[], 'novel_cds':[]}
mean_phylop={'case':{'introns':[], 'cds':[], 'utr':[], 'novel_cds':[]}, 'control':{'introns':[], 'cds':[], 'utr':[], 'novel_cds':[]}}
mean_cadd={'case':{'introns':[], 'cds':[], 'utr':[], 'novel_cds':[]}, 'control':{'introns':[], 'cds':[], 'utr':[], 'novel_cds':[]}}
for i in tqdm(range(100)):
    tmp=dnm.sample(len(dnm), replace=True).intersect(gencode_introns).as_df().drop_duplicates()
    count['introns']+=[(np.sum(tmp.PrimaryPhenotype!="control")/516)/(np.sum(tmp.PrimaryPhenotype=="control")/516)]
    mean_phylop['case']['introns']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['autism']]
    mean_phylop['control']['introns']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['control']]
    mean_cadd['case']['introns']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['autism']]
    mean_cadd['control']['introns']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['control']]
    tmp=dnm.sample(len(dnm), replace=True).intersect(gencode_utr).as_df().drop_duplicates()
    count['utr']+=[(np.sum(tmp.PrimaryPhenotype!="control")/516)/(np.sum(tmp.PrimaryPhenotype=="control")/516)]
    mean_phylop['case']['utr']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['autism']]
    mean_phylop['control']['utr']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['control']]
    mean_cadd['case']['utr']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['autism']]
    mean_cadd['control']['utr']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['control']]
    tmp=dnm.sample(len(dnm), replace=True).intersect(gencode_cds).as_df().drop_duplicates()
    count['cds']+=[(np.sum(tmp.PrimaryPhenotype!="control")/516)/(np.sum(tmp.PrimaryPhenotype=="control")/516)]
    mean_phylop['case']['cds']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['autism']]
    mean_phylop['control']['cds']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['control']]
    mean_cadd['case']['cds']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['autism']]
    mean_cadd['control']['cds']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['control']]
    tmp=dnm.sample(len(dnm), replace=True).intersect(novel_cds).as_df().drop_duplicates()
    if len(tmp)==0:
        count['novel_cds']+=[np.inf]
    elif tmp[tmp['PrimaryPhenotype']=='control'].shape[0]==0:
        count['novel_cds']+=[np.inf]
    else:
        count['novel_cds']+=[(np.sum(tmp.PrimaryPhenotype!="control")/516)/(np.sum(tmp.PrimaryPhenotype=="control")/516)]
        mean_phylop['case']['novel_cds']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['autism']]
        mean_phylop['control']['novel_cds']+=[tmp.groupby('PrimaryPhenotype')['phylop'].mean().loc['control']]
        mean_cadd['case']['novel_cds']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['autism']]
        mean_cadd['control']['novel_cds']+=[tmp.groupby('PrimaryPhenotype')['CaddScore'].mean().loc['control']]


print([dnm.intersect(gencode_introns).as_df().drop_duplicates()['PrimaryPhenotype'].value_counts(), (np.sum(dnm.intersect(gencode_introns).as_df().drop_duplicates()['PrimaryPhenotype']!="control")/516)/(np.sum(dnm.intersect(gencode_introns).as_df().drop_duplicates()['PrimaryPhenotype']=="control")/516), np.quantile(count['introns'], 0.05), np.quantile(count['introns'], 0.95)])
print([dnm.intersect(gencode_cds).as_df().drop_duplicates()['PrimaryPhenotype'].value_counts(), (np.sum(dnm.intersect(gencode_cds).as_df().drop_duplicates()['PrimaryPhenotype']!="control")/516)/(np.sum(dnm.intersect(gencode_cds).as_df().drop_duplicates()['PrimaryPhenotype']=="control")/516), np.quantile(count['cds'], 0.05), np.quantile(count['cds'], 0.95)])
print([dnm.intersect(gencode_utr).as_df().drop_duplicates()['PrimaryPhenotype'].value_counts(), (np.sum(dnm.intersect(gencode_utr).as_df().drop_duplicates()['PrimaryPhenotype']!="control")/516)/(np.sum(dnm.intersect(gencode_utr).as_df().drop_duplicates()['PrimaryPhenotype']=="control")/516), np.quantile(count['utr'], 0.05), np.quantile(count['utr'], 0.95)])
print([dnm.intersect(novel_cds).as_df().drop_duplicates()['PrimaryPhenotype'].value_counts(), (np.sum(dnm.intersect(novel_cds).as_df().drop_duplicates()['PrimaryPhenotype']!="control")/516)/(np.sum(dnm.intersect(novel_cds).as_df().drop_duplicates()['PrimaryPhenotype']=="control")/516),np.quantile(count['novel_cds'], 0.05), np.quantile(count['novel_cds'], 0.95)])

print([[dnm.intersect(gencode_cds).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['autism'], np.quantile(mean_phylop['case']['cds'], 0.05), np.quantile(mean_phylop['case']['cds'], 0.95)], [dnm.intersect(gencode_cds).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['control'],np.quantile(mean_phylop['control']['cds'], 0.05), np.quantile(mean_phylop['control']['cds'],0.95)]])
print([[dnm.intersect(gencode_utr).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['autism'], np.quantile(mean_phylop['case']['utr'], 0.05), np.quantile(mean_phylop['case']['utr'], 0.95)], [dnm.intersect(gencode_utr).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['control'],np.quantile(mean_phylop['control']['utr'], 0.05), np.quantile(mean_phylop['control']['utr'],0.95)]])
print([[dnm.intersect(gencode_introns).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['autism'], np.quantile(mean_phylop['case']['introns'], 0.05), np.quantile(mean_phylop['case']['introns'], 0.95)], [dnm.intersect(gencode_introns).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['control'],np.quantile(mean_phylop['control']['introns'], 0.05), np.quantile(mean_phylop['control']['introns'],0.95)]])
print([[dnm.intersect(novel_cds).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['autism'], np.quantile(mean_phylop['case']['novel_cds'], 0.05), np.quantile(mean_phylop['case']['novel_cds'], 0.95)], [dnm.intersect(novel_cds).as_df().drop_duplicates().groupby('PrimaryPhenotype')['phylop'].mean().loc['control'],np.quantile(mean_phylop['control']['novel_cds'], 0.05), np.quantile(mean_phylop['control']['novel_cds'],0.95)]])
