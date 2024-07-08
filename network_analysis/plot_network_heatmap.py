# written 2-7-24 JLW
# to visualize a heatmap of network proteins
# and classes

import pickle,os,csv,sys, math
import pandas as pd
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# get ATC codes
dbdf = pd.read_excel('../../PathFXv2/data/Drugbank050120.xlsx')
d2atc = defaultdict(set)
for (d,a_str) in zip(dbdf.drugbank_id,dbdf.atc_codes):
	#print((d,a_str))
	#print(type(a_str))
	if type(a_str) == type(1.0):
		continue
	if '|' in a_str:
		atc_codes = a_str.split('|')
	else:
		atc_codes = [a_str]
	# print(atc_codes)
	for a in atc_codes:
		keep_a = a[0]
		d2atc[d].add(keep_a)

print('gathered atc codes')


# get approved drug list
app_drugs = set() # 3988 unique approved drugs
for (db,dg) in zip(dbdf.drugbank_id,dbdf.groups):
	for g in dg.split('|'):
		if g=='approved':
			app_drugs.add(db)

# load pathway results
f = 'all_drug_summary_als_san_0623.xlsx'
rdf = pd.read_excel(f)
r_dbids = set(rdf.DrugBankID) # 2984 drugs

# get all network proteins, define group classes
all_net_prot = set([g for gstr in rdf['NetworkGenes'].to_list() for g in gstr.split(',')])

cs_genes = ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C3', 'C3AR1', 'C5', 'C5AR1']
npy_genes = ['NPY','NPY5R', 'NPY2R', 'NPY1R']
cx_genes = [g for g in all_net_prot if 'CXC' in g or 'CCR' in g]

# create color map and mapping for network proteins of interest
net_classes = {'None':0,'Any':1,'CXCR5':2,'CXCR3':3,'Chemokine':4,'CS':5,'NPY,R':6,'CNR2':7}
net_class_colors = ['white','black','pink','lime','red','yellow','blue','green',] 
cmap = LinearSegmentedColormap.from_list('NetClass',net_class_colors,len(net_class_colors))

# loop through predictions and save phenotype-associated proteins
all_rows = []
unique_dbids = set(rdf.DrugBankID)
# for (db,ngs) in zip(rdf.DrugBankID,rdf.NetworkGenes):
for db in unique_dbids:
	if db not in app_drugs:
		continue
	df_short = rdf[rdf['DrugBankID']==db]
	row_data = {'DrugBankID':db}
	for ngs in df_short.NetworkGenes:
		for g in ngs.split(','):
			if g in cs_genes:
				row_data[g] = net_classes['CS']
			elif g in cx_genes:
				row_data[g] = net_classes['Chemokine']
			elif g in npy_genes:
				row_data[g] = net_classes['NPY,R']
			elif g in net_classes:
				row_data[g] = net_classes[g]
			else:
				row_data[g] = 1

	all_rows.append(row_data)

plot_data = pd.DataFrame(all_rows).fillna(0)
# print(plot_data.head())
plot_data = plot_data.set_index('DrugBankID')
plot_data = plot_data.loc[:,(plot_data.sum()>2)] # require at least two drug networks to be plotted
print('created plot data frame')

# map ATC codes to colors, drugs to colors, then make into dataframe
pfx_atc = dict([x for x in d2atc.items() if x[0] in r_dbids]) # 897 drugs with ATC cods (i.e. approved)
all_atc = sorted(set([x for aset in pfx_atc.values() for x in aset])) # associated with 14 unique ATC codes
atc_colors =['red','coral','peru','darkorange','gold','yellowgreen','green','lightseagreen','dodgerblue','slateblue','mediumorchid','violet','hotpink','lightpink'] 
atc_to_color = dict(zip(all_atc,atc_colors))
drug_to_color = dict([(d,atc_to_color[ac]) if len(ac_set)==1 else (d,'black') for (d,ac_set) in pfx_atc.items() for ac in ac_set]) # mixed class drugs in black,  drug-ATC codes, 3150 unique drugs
(dnames,dcolors) = zip(*[(k,v) for (k,v) in drug_to_color.items()])
rc_df = pd.DataFrame.from_dict({'Drugs': dnames,'ATC Codes': dcolors})
rc_df = rc_df.set_index('Drugs')

# bar chart of all ATC codes
atc_count = defaultdict(int)
for (d,aset) in pfx_atc.items():
	for acode in aset:
		atc_count[acode]+=1

(alabels,acounts) = zip(*sorted(atc_count.items()))
fig,ax = plt.subplots()
y_pos = np.arange(len(acounts))
bar_colors = [c for (a,c) in sorted(atc_to_color.items())]
hbars = ax.barh(y_pos,acounts,align='center',color=bar_colors)
ax.set_yticks(y_pos)
ax.set_yticklabels(alabels)
# old school way to force bar labels
for (yp,ac) in zip(y_pos,acounts):
	n=ax.text(ac+5, yp, str(ac))
ax.set_title('Number of active ingredients \nper level-1 ATC code')
ax.set_ylim([-1,14])
ax.get_xaxis().set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()  # labels read top-to-bottom
plt.savefig('Num_active_ingred_level1_atc_codes.png',format='png',dpi=300)

# plot using code that I know works
fig,ax = plt.subplots()
g = sns.clustermap(plot_data,cmap=cmap,yticklabels=False,xticklabels=False,row_colors=rc_df,cbar_kws={'label':'DrugClass','ticks':[0,1,2,3,4,5,6,7],},)
ax.set_xlabel('Network Proteins')
ax.set_ylabel('Drugs')
ax.set_title('Drug network proteins')

plt.savefig('approved_drugs_with_net_clusters.png',format='png')
print('saving figure')

# make a 2nd heatmap just of CS system drugs
# repeat code for selecting clusters to get DrugBankIds

# look at drugs associated with a particular CUI for defining network classes
drugs_by_cui = defaultdict(set)
for (dbid,cterm,gene_list) in zip(rdf.DrugBankID,rdf.cui,rdf.NetworkGenes):
        drugs_by_cui[cterm].add(dbid)

def get_cluster_data_frame(glist):
	gdf = rdf[rdf['NetworkGenes'].str.contains('|'.join(glist))]
	keep_drugs = set()
	for (d,gene_str) in zip(gdf.DrugBankID,gdf.NetworkGenes):
		genes = gene_str.split(',')
		if len(set(genes).intersection(glist)) >=1:
			keep_drugs.add(d)
	gdf = gdf[gdf.DrugBankID.isin(keep_drugs)]
	# just keep approved drugs, find comparator drugs by connection to CUIs
	gclus = set(gdf.DrugBankID).intersection(app_drugs)
	nong_clus = app_drugs.difference(gclus)
	gcui_terms = set(gdf.cui)
	cui_assoc_drugs = set([d for c in gcui_terms for d in drugs_by_cui[c]])
	nong_cui_assoc = nong_clus.intersection(cui_assoc_drugs)
	return (gclus,nong_cui_assoc)

def plot_class_htmp(glist,sname,short_name):
	# slice old way, and then force exact matches
	# (glist,sname,short_name) = (cs_genes,'Complement System','CS')
	print('cluster-specific heatmap ',short_name)
	(gclus,nong_cui_assoc) = get_cluster_data_frame(glist)

	# assign row colors
	drug_to_color = dict([(d,'lightskyblue') for d in gclus] + [(d,'darkgrey') for d in nong_cui_assoc]) 
	(dnames,dcolors) = zip(*[(k,v) for (k,v) in drug_to_color.items()])
	rc_df = pd.DataFrame.from_dict({'Drugs': dnames,'NetClass': dcolors})
	rc_df = rc_df.set_index('Drugs')

	# subset from other plot to keep the same values
	g_plot_data = plot_data[plot_data.index.isin(drug_to_color)].fillna(0)
	g_plot_data = g_plot_data.loc[:,(g_plot_data.sum()>1)]

	# highlight class of interest, this step can be used once in a loop
	base_colors = ['white','black','black','black','black','black','black','black',]
	update_index = list(net_classes).index(short_name)
	update_color = net_class_colors[update_index]
	base_colors[update_index] = update_color
	cmap = LinearSegmentedColormap.from_list('NetClass',base_colors,len(base_colors))

	# plot most common net prots
	fig,ax = plt.subplots()
	g = sns.clustermap(g_plot_data,cmap=cmap,yticklabels=False,xticklabels=False,row_colors=rc_df,cbar_kws={'label':'DrugClass','ticks':[0,1,2,3,4,5,6,7],},)
	ax.set_xlabel('Network Proteins')
	ax.set_ylabel('Drugs')
	ax.set_title('Drug network proteins\n' + sname + ' Class')
	plt.savefig('NetClassHeatmap_'+short_name+'.png',format='png')

	# repeat for just class prots
	base_colors = ['white',update_color]
	cmap = LinearSegmentedColormap.from_list('NetClass',base_colors,len(base_colors))
	class_col = [c for c in g_plot_data.columns if c in glist]
	gclass_plot = g_plot_data[class_col]
	fig,ax = plt.subplots()
	if len(class_col) > 1:
		g = sns.clustermap(gclass_plot,cmap=cmap,yticklabels=False,row_colors=rc_df,)
	else:
		g = sns.clustermap(gclass_plot,cmap=cmap,yticklabels=False,row_colors=rc_df,col_cluster=False,)
		
	ax.set_xlabel('Network Proteins')
	ax.set_ylabel('Drugs')
	ax.set_title(sname+' Proteins')
	plt.savefig('JustClassProt_'+short_name+'.png',format='png')

# eventually loop through all of them
plot_class_htmp(cs_genes,'Complement System','CS')
plot_class_htmp(cx_genes,'Chemokine System','Chemokine')
plot_class_htmp(['NPY','NPY1R','NPY2R','NPY5R'],'Neuropeptide and Receptors','NPY,R')
plot_class_htmp(['CXCR5'],'CXCR5','CXCR5')
plot_class_htmp(['CXCR3'],'CXCR3','CXCR3')
plot_class_htmp(['CNR2'],'CNR2','CNR2')




# dict keys: 'CXCR5':2,'CXCR3':3,'Chemokine':4,'CS':5,'NPY,R':6,'CNR2':7}
# developmet, doesn't all work
#cbar_ax = fig.axes[-1]
#cbar_ax.set_xticklabels(['none','any', 'CXCR5', 'CXCR3','Chemokine', 'CS', 'NPY', 'CNR2'])
#cbar = fig.colorbar(g)
#cbar.set_ticklabels(['any', 'CXCR5', 'CXCR3','Chemokine', 'CS', 'NPY', 'CNR2'])
# g.cax.set_visible(False)
# test out code that I'm not sure works
# Manually specify colorbar labelling after it's been generated
# colorbar = g.collections[0].colorbar
# colorbar.set_ticks([0,1,2,3,4,5,6,7])
# colorbar.set_ticklabels(['none','any', 'CXCR5', 'CXCR3','Chemokine', 'CS', 'NPY', 'CNR2'])


