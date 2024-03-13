# written to search DrugBank results
# updated 6-20-23 JLW for new ALS pathways

import pickle,os,csv
from collections import defaultdict
import pandas as pd
import argparse

# first gather relevant phenotypes
## Move this to the call_search_phen script
#p2c = pickle.load(open('../rscs/Pfx050120_all_phens_to_cuis.pkl','rb'))
d2n = pickle.load(open('../rscs/pfxDB050620_dbid2name.pkl','rb'))
n2db = dict([(v.lower(),k) for (k,v) in d2n.items()])
# gather arguments
parser = argparse.ArgumentParser(description='Parameters for the analysis')
parser.add_argument('-aname', action='store', dest='aname', help='Analysis Name, no spaces')
parser.add_argument('-rdir', action='store', dest='rdir', help='The path of the directory where the cui list is pickled')
parser.add_argument('-cuipkl', action='store', dest='cpkl', help='The pickled file name containing the list of CUI terms')
args = parser.parse_args()

rdir = args.rdir
keep_cuis = pickle.load(open(os.path.join(rdir,args.cpkl),'rb'))
aname = args.aname
out_excel_name = os.path.join(rdir,"all_drug_summary_"+aname+".xlsx") 
writer = pd.ExcelWriter(out_excel_name, engine='xlsxwriter')

# gather all association table files
rdir = '../results/alldrugbank_v4/' # updated for v4 results
allf = [(ssd,sd,flist) for (ssd,sd,flist) in os.walk(rdir)]
asf_dic = dict([(os.path.split(ssd)[-1],os.path.join(ssd,f)) for (ssd,sd,flist) in allf for f in flist if '_merged_neighborhood__assoc_table_.txt' in f])

print('Reading all Drug Bank')
res_df = pd.DataFrame()
for (dbid,dfile) in asf_dic.items():
	dname = d2n[dbid]
	df = pd.read_table(dfile)
	if not df.empty:
		df['DrugName'] = dname.lower()
		df["DrugBankID"] = dbid
		df_short = df[df['cui'].isin(keep_cuis)]
		if not df_short.empty:
#			print(dname)
			res_df = pd.concat([res_df,df_short])
	 
# add drug targets as FYI
print('Adding Drug Target and Description Info')
dint = pickle.load(open('../rscs/pfxDB050620_dint.pkl','rb'))
dint_str = dict([(d,','.join(glist)) for (d,glist) in dint.items()])
res_df['DrugTargetProteins'] = res_df['DrugBankID'].map(dint_str)
res_df.drop(['rank','assoc in neigh', 'assoc in intom','probability','Unnamed: 8'], axis=1)
ord_cols = ['DrugName',"DrugBankID",'DrugTargetProteins','phenotype','Benjamini-Hochberg', 'cui','genes']
res_df = res_df = res_df[ord_cols]
res_df = res_df.rename(columns={"genes": "NetworkGenes",})


# add some more info from DrugBank as fyi
dbank = pd.read_excel('../data/Drugbank050120.xlsx')
dbank = dbank.rename(columns={'drugbank_id':'DrugBankID'})

res_df = res_df.merge(dbank[['DrugBankID','type','description']],on='DrugBankID')
print(res_df.head)

# look up labeled safety concerns, convert to DB id
print('Looking up safety data')
# These are DMEs from the drug's labels
db_to_dme_str = pickle.load(open('../data/drugbankid_to_labeled_dme_str.pkl','rb'))
res_df['LabeledSafetyWarnings'] = res_df['DrugBankID'].map(db_to_dme_str)
res_df.to_excel(writer,index=False,sheet_name="AllPredicted")
relevant_DBids = list(set(res_df["DrugBankID"].to_list()))


## added to summarize number of drugs associated with multiple phenotypes
## count drug and phenotype associations
#drug_phs = defaultdict(set)
#for (dn,db,ph) in zip(res_df.DrugName, res_df.DrugBankID, res_df.phenotype):
#        drug_phs[(dn,db)].add(ph)
#
## assemble into data table
#all_rows = []
#for ((dn,db),ph_set) in drug_phs.items():
#        ph_str = ','.join(sorted(ph_set))
#        ph_count = len(ph_set)
#        row_data = {'DrugName':dn,'DrugBankID':db,'PhenotypeCount':ph_count,'Phenotypes':ph_str}
#        all_rows.append(row_data)
#
## create table and save to excel
#ph_count_df = pd.DataFrame(all_rows)
#ph_count_df = ph_count_df.sort_values(by='PhenotypeCount', ascending=False)
#ph_count_df.to_excel(writer,index=False,sheet_name="DrugsByPhenNum")
#
## these are predicted from PathFX
## pred_dme_df = pickle.load(open('../data/predicted_pathfx_dmes.pkl','rb'))
#pred_dme_df = pd.read_excel('../data/predicted_pathfx_dmes.xlsx')
#pred_dme_df = pred_dme_df[pred_dme_df["DrugBankID"].isin(relevant_DBids)]
#pred_dme_df.to_excel(writer,index=False,sheet_name="PathFXpredDME")


# look at common network mechanisms
net_genes = res_df['NetworkGenes'].to_list()
gene_counts = defaultdict(int)
for gentry in net_genes:
	glist = gentry.split(',')
	for g in glist:
		gene_counts[g]+=1

# keep the top 20 genes associated wih a drug
# pickle.dump(gene_counts,open(os.path.join(rdir,'top_net_genes_all.pkl'),'wb'))
top_mechs = sorted([x for x in gene_counts.items()],key=lambda x:x[1],reverse=True)[0:20]
top_genes = [x[0] for x in top_mechs]
keep_gene_strs = [gentry for gentry in net_genes for tg in top_genes if tg in gentry]
count_top_mechs = defaultdict(int)
for gentry in list(set(net_genes)):
	for tg in top_genes:
		if tg in gentry:
			count_top_mechs[gentry]+=1
top_mech_df = res_df[res_df['NetworkGenes'].isin(keep_gene_strs)]
top_mech_df['TopMechCount'] = top_mech_df['NetworkGenes'].map(count_top_mechs)
top_mech_df.to_excel(writer,index=False,sheet_name="TopPredicted")

top_mech_sum = pd.DataFrame(top_mechs,columns=['GeneName','Occured in X drug networks'])
top_mech_sum.to_excel(writer,index=False,sheet_name="SummOfMostCommonGenes")

writer.save()
print('Results saved in '+out_excel_name)


