# written to see if protein classes share network phenotypes
# written 3-24 JLW

import pandas as pd
import pickle,os
from collections import defaultdict
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import seaborn as sns

allf = [f for f in os.listdir('.') if '_pathways_0224.xlsx' in f]

prot_to_cuis = defaultdict(set)
cuis_names = []
for f in allf:
	pclass = f.replace('_pathways_0224.xlsx','')
	df = pd.read_excel(f)
	prot_to_cuis[pclass] = set(df.cui)
	cuis_names+=list(zip(df.cui,df.phenotype))

cuis_names_dic = dict(set(cuis_names))

all_rows = []
for (p1,cuis1) in prot_to_cuis.items():
	row_data = {'index':p1}
	for (p2,cuis2) in prot_to_cuis.items():
		num_inter = len(cuis1.intersection(cuis2))
		num_union = len(cuis1.union(cuis2))
		row_data[p2] = num_inter/float(num_union)
	all_rows.append(row_data)

pdf = pd.DataFrame(all_rows)
pdf = pdf.set_index('index')

fig,ax = plt.subplots()
sns.clustermap(pdf,cmap="YlOrBr")
ax.set_title('Shared disease pathways by protein classes')
plt.savefig(os.path.join('prot_class_shared_cuis.png'),format='png')

# write the phenotypes to excel
all_rows = []
for (p,cset) in prot_to_cuis.items():
	phen_list = [cuis_names_dic[c] for c in sorted(cset)]
	c_str = ' | '.join(sorted(cset))
	p_str = ' | '.join(phen_list)
	all_rows.append({'protein': p, 'cuis':c_str, 'disease_pathways':p_str})

cdf = pd.DataFrame(all_rows)
cdf.to_excel('net_classes_to_disease_pathways.xlsx',index=False)

prot_to_cuis['chemokine'].difference(prot_to_cuis['complement_system'])
#{'C0751967', 'C1862941', 'C3160718', 'SANms001', 'C0276496', 'C0007789', 'C0494463', 'TALS006'}
#[cuis_names_dic[c] for c in {'C0751967', 'C1862941', 'C3160718', 'SANms001', 'C0276496', 'C0007789', 'C0494463', 'TALS006'}
#['Amyotrophic Lateral Sclerosis, Sporadic', 'Multiple Sclerosis, Relapsing-Remitting', 'Parkinson disease, late-onset', 'MS CSF', 'Familial Alzheimer Disease (FAD)', 'Cerebral Palsy', 'Alzheimer Disease, Late Onset', 'Target ALS log2fc>1 SpinalCord']

prot_to_cuis['complement_system'].difference(prot_to_cuis['chemokine'])
#{'SANsn10', 'TALS007', 'SANsn9'}
#[cuis_names_dic[c] for c in {'SANsn10', 'TALS007', 'SANsn9'}]
#['pericytes_SMC snRNAseq pseudobulk', 'Target ALS log2fc>1 SpinalCordCervical', 'oligodendrocyte snRNAseq pseudobulk']
