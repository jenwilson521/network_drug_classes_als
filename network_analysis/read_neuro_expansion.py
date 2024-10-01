# written to assess pathways discovered
# in the neuroexpansion
# writte 9-8-23 JLW

import pickle,os,csv
import pandas as pd
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

rdir = '../results/neuro_expansion/'
rf = os.path.join(rdir,'all_drug_summary_neuro_expansion.xlsx')
df = pd.read_excel(rf) # df.shape (4699, 10)

udrugs = len(set(df.DrugBankID)) # 1257
ucuis = len(set(df.cui)) # 20

df.groupby('phenotype')['phenotype'].count()
#ALZHEIMER DISEASE 5                          13
#Alzheimer Disease                            99
#Alzheimer Disease, Early Onset              170
#Alzheimer Disease, Late Onset               336
#Alzheimer disease                          1012
#Alzheimer's Disease, Focal Onset            266
#Amyotrophic Lateral Sclerosis, Sporadic       6
#Amyotrophic lateral sclerosis                58
#Amyotrophic lateral sclerosis type 1          4
#Cerebral Infarction                         549
#Familial Alzheimer Disease (FAD)            442
#Ischemic stroke                             329
#Multiple Sclerosis                          534
#Multiple Sclerosis, Acute Fulminating         7
#Multiple Sclerosis, Relapsing-Remitting       1
#Myasthenia Gravis                           376
#Parkinson disease                           421
#Parkinson disease 2                          39
#Parkinson disease, late-onset                12
#Parkinsonian Disorders                       25

# group diseases 
dis_groups = ['AZ','AZ','AZ','AZ','AZ','AZ','ALS','ALS','ALS','STK','AZ','STK','MS','MS','MS','MYGR','PD','PD','PD','PD']
ph2gp = dict(zip(sorted(set(df.phenotype)),dis_groups))
df['disGroup'] = df.phenotype.map(ph2gp)

# Jaccard similarity between these diseases
# based on drug networks
dg2dbs = defaultdict(set)
for (db,dg) in zip(df['DrugBankID'],df['disGroup']):
	dg2dbs[dg].add(db)

all_rows = []
for (d1,drugs1) in dg2dbs.items():
	row_data = {'index':d1}
	for (d2,drugs2) in dg2dbs.items():
		num_inter = len(drugs1.intersection(drugs2))
		num_union = len(drugs1.union(drugs2))
		row_data[d2] = num_inter/float(num_union)
	all_rows.append(row_data)

pdf = pd.DataFrame(all_rows)
pdf = pdf.set_index('index')

sns.set(font_scale = 2)
fig,ax = plt.subplots()
sns.clustermap(pdf,cmap="binary").fig.suptitle('Disease groups by co-occurrence in drug networks')
ax.set_title('Shared drugs by disease groups')
plt.savefig(os.path.join(rdir,'neuro_expan_drugs_shared.png'),format='png')

# shared genes for drugs with each disease?
# do GO enrichment of shared genes
dg2netGenes = defaultdict(set)
for (dg,ng_str) in zip(df['disGroup'],df['NetworkGenes']):
	for g in ng_str.split(','):
		dg2netGenes[dg].add(g)

all_rows = []
for d1 in list(dg2netGenes.keys()):
	row_data = {'DisCat1':d1}
	for d2 in list(dg2netGenes.keys()):
		d1genes = dg2netGenes[d1]
		d2genes = dg2netGenes[d2]
		shared = d1genes.intersection(d2genes)
		row_data[d2] = ','.join(sorted(shared))
		outf = open(os.path.join(rdir, "_".join([d1,d2,'sharedNetGenes']) + ".txt"),"w")
		n = outf.write('\n'.join(sorted(shared)))
		outf.close()
	all_rows.append(row_data)

gdf = pd.DataFrame(all_rows)
gdf.to_excel(os.path.join(rdir,"shared_net_genes_dis_cat_.xlsx"),index=False)

# repeat for plotting
all_rows = []
for (d1,d1genes) in list(dg2netGenes.items()):
	row_data = {'index':d1}
	for (d2,d2genes) in list(dg2netGenes.items()):
		num_inter = len(d1genes.intersection(d2genes))
		num_union = len(d1genes.union(d2genes))
		row_data[d2] = num_inter/float(num_union)
	all_rows.append(row_data)

gdf = pd.DataFrame(all_rows)
gdf = gdf.set_index('index')

fig,ax = plt.subplots()
sns.clustermap(gdf,cmap='binary').fig.suptitle('Disease groups by shared drug network genes')
plt.savefig(os.path.join(rdir,'neuro_expan_disGroup_shared_net_genes.png'),format='png')

# look up all CUI genes, use 
# repeat, but now use all pathway genes, not just those in the networks
d2cuiGenes = defaultdict(set)
c2g = pickle.load(open('../rscs/Pfx050120_merged_unique_cuis2genes.pkl','rb'))
disGroups_cuis = set(zip(df['disGroup'],df['cui'])) #df.cui
for (dis_group,dis_cui) in disGroups_cuis:
	cui_genes = c2g[dis_cui]
	d2cuiGenes[dis_group]=set(cui_genes)

# calculate similarity
all_rows = []
for (d1,d1genes) in list(d2cuiGenes.items()):
        row_data = {'index':d1}
        for (d2,d2genes) in list(d2cuiGenes.items()):
                num_inter = len(d1genes.intersection(d2genes))
                num_union = len(d1genes.union(d2genes))
                row_data[d2] = num_inter/float(num_union)
        all_rows.append(row_data)

cgdf = pd.DataFrame(all_rows)
cgdf = cgdf.set_index('index')

fig,ax = plt.subplots()
sns.clustermap(cgdf,cmap='binary').fig.suptitle('Disease groups by shared pathway genes') 
plt.savefig(os.path.join(rdir,'neuro_expan_disGroup_shared_allPathway_genes.png'),format='png')
