# modified from the original version to search v4
# and have a more streamlined code base
# written 6-20-23 JLW

import pickle,os,csv
import pandas as pd
from collections import defaultdict

##### DRUGBANK APPROVAL, SYNONYMS #######
# get approved status from DrugBank file
dbf = '../data/Drugbank050120.xlsx'
dbdf = pd.read_excel(dbf)
approved_dbids = set()
dbid_to_groups = dict()
for (dbid,dgroups) in zip(dbdf.drugbank_id,dbdf.groups):
        dbid_to_groups[dbid] = dgroups
        if 'approved' in dgroups:
                approved_dbids.add(dbid)

# data object for later
dbid2name = dict(zip(dbdf.drugbank_id,dbdf.name))

# save synonyms
dbvocabf = '../data/drugbank_vocabulary.csv'
dbvocab = pd.read_csv(dbvocabf)
db2syns = defaultdict(set)
for (db,syn_ent) in zip(dbvocab['DrugBank ID'].to_list(),dbvocab['Synonyms'].to_list()):
#       print((db,syn_ent))
        if type(syn_ent)==type(1.):
                continue
        if "|" in syn_ent:
                syn_list = syn_ent.split(" | ")
                for s in syn_list:
                        db2syns[db].add(s)
        else:   
                db2syns[db].add(syn_ent)


#### ALL DRUGBANK RESULTS ####
rdir = 'UPDATE'
rf = os.path.join(rdir,'all_drug_summary_als_san_0623.xlsx')
rdf = pd.read_excel(rf)

# all network proteins
all_net_prot = set([g for gstr in rdf['NetworkGenes'].to_list() for g in gstr.split(',')])

#### collect phenoptypes supported by each network protein ####
genes_to_cuis = defaultdict(set)
drugs_by_cui = defaultdict(set)
cuis_by_drugs = defaultdict(set)
for (dbid,cterm,gene_list) in zip(rdf.DrugBankID,rdf.cui,rdf.NetworkGenes):
        drugs_by_cui[cterm].add(dbid)
        cuis_by_drugs[dbid].add(cterm)
        for net_gene in gene_list.split(','):
                genes_to_cuis[net_gene].add(cterm)

# pre-select for approved drugs
ardf = rdf[rdf.DrugBankID.isin(approved_dbids)]

def gene_in_str(gstr,g):
	gene_list = gstr.split(',')
	stat = False
	if g in gene_list:
		stat = True
	return stat

all_rows_long = []
for net_gene in all_net_prot:
	net_gene_slice = rdf[rdf['NetworkGenes'].str.contains(net_gene)]
	net_gene_slice['GeneActuallyThere'] = net_gene_slice.apply(lambda row: gene_in_str(row['NetworkGenes'],net_gene),axis=1)
	net_gene_slice = net_gene_slice[net_gene_slice['GeneActuallyThere']==True]

	num_drug_cui_pairs = net_gene_slice.shape[0]
	# num_san_paths = net_gene_slice[net_gene_slice['cui'].str.contains('SAN')].shape[0]
	net_gene_drugs = set(net_gene_slice.DrugBankID)
	non_net_gene_drugs = set(rdf.DrugBankID).difference(net_gene_drugs)
	# in cluster, format for output
	app_net_gene_drugs = net_gene_drugs.intersection(approved_dbids)
	num_app_in_cluster = len(app_net_gene_drugs)
	clus_app_dnames = sorted([dbid2name[x] for x in list(app_net_gene_drugs)])
	# non-gene cluster, format for output
	app_non_net_gene_drugs = non_net_gene_drugs.intersection(approved_dbids)
	num_app_not_in_cluster = len(app_non_net_gene_drugs)
	nonclus_app_dnames = sorted([dbid2name[x] for x in list(app_non_net_gene_drugs)])
	row_data = {'GeneName':net_gene,'NumDrugPathAssocs':num_drug_cui_pairs,'NumAppInCluster':num_app_in_cluster,'DrugsInCluster':','.join(clus_app_dnames),'NumAppOutOfCluster':num_app_not_in_cluster,'DrugsOutOfCluster':','.join(nonclus_app_dnames)}
	all_rows_long.append(row_data)


def ratio_fxn(nin,nout):
	return float(nin+1)/float(nout+1)

all_drug_clusters = pd.DataFrame(all_rows_long)
# all_drug_clusters['Ratio'] = all_drug_clusters.apply(lambda row: ratio_fxn(all_drug_clusters['NumAppInCluster'],all_drug_clusters['NumAppOutOfCluster']),axis=1)
all_drug_clusters.to_excel(os.path.join(rdir,'all_net_prot_clusters_0224.xlsx'),index=False)


### method for extracting cluster info for gene pathways of interest ###
# select pathway proteins
complement_proteins = [g for g in all_net_prot if g[0:2]=='C1' or g[0:2]=='C3' or g[0:2]=='C5'] # ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C3', 'C3AR1', 'C5', 'C5AR1']
chemokine_proteins = [g for g in all_net_prot if 'CXC' in g or 'CCR' in g]

def check_bad_str(dstr):
# clean out strings that cannot be mapped
        stat = False
        if '(' in dstr:
                stat=True
        if "\'" in dstr:
                stat = True
        if "-" in dstr:
                stat = True
        return stat

def get_drug_syns_list(dblist):
	all_names = []
	for d in dblist:
		all_names.append(dbid2name[d])
		if d in db2syns:
			dsyns = db2syns[d]
			all_names+=dsyns
	clean_names = [n for n in all_names if not check_bad_str(n)]
	full_list = sorted(clean_names)
	return ','.join(full_list)	

def char_single_gene(gname):
	gdf = rdf[rdf['NetworkGenes'].str.contains(gname)]
	glist = [gname]
	keep_drugs = set()
	error_drugs = set()
	for (d,gene_str) in zip(gdf.DrugBankID,gdf.NetworkGenes):
		genes = gene_str.split(',')
		if len(set(genes).intersection(glist)) >=1:
			keep_drugs.add(d)
		else:
			error_drugs.add(d)
	errors_not_repeats = error_drugs.difference(keep_drugs)
	gdf = gdf[gdf.DrugBankID.isin(keep_drugs)]
	gdf.to_excel(os.path.join(rdir,gname+'_pathways_0224.xlsx'),index=False)
	gclus = set(gdf.DrugBankID).intersection(approved_dbids) 
	nong_clus = approved_dbids.difference(gclus)
	gcui_terms = set(gdf.cui) 
	cui_assoc_drugs = set([d for c in gcui_terms for d in drugs_by_cui[c]])
	nong_cui_assoc = nong_clus.intersection(cui_assoc_drugs) 
	print('No of drug-cui pairs: ',gdf.shape[0])
	print('No of cuis: ',len(gcui_terms))
	print('No of in-class, approved: ',len(gclus))
	print('No of cui-assoc, approved, out-of-class: ',len(nong_cui_assoc))
	g_drug_str = get_drug_syns_list(gclus)
	nong_drug_str = get_drug_syns_list(nong_cui_assoc)
	outf = open(os.path.join(rdir,gname+'_drug_strings_0224.txt'),'w')
	n = outf.write(gname+"_syns=\'" + g_drug_str + "\'\n" + "non_"+gname+"_syns=\'" + nong_drug_str + "\'")
	outf.close()
	return (g_drug_str,nong_drug_str)


def char_gene_list(glist,sname,short_name):
	gdf = rdf[rdf['NetworkGenes'].str.contains('|'.join(glist))]
	keep_drugs = set()
	error_drugs = set()
	for (d,gene_str) in zip(gdf.DrugBankID,gdf.NetworkGenes):
		genes = gene_str.split(',')
		if len(set(genes).intersection(glist)) >=1:
			keep_drugs.add(d)
		else:
			error_drugs.add(d)
	errors_not_repeats = error_drugs.difference(keep_drugs)

	gdf = gdf[gdf.DrugBankID.isin(keep_drugs)]

	gdf.to_excel(os.path.join(rdir,sname+'_pathways_0224.xlsx'),index=False)
	gclus = set(gdf.DrugBankID).intersection(approved_dbids) 
	nong_clus = approved_dbids.difference(gclus)
	gcui_terms = set(gdf.cui) 
	cui_assoc_drugs = set([d for c in gcui_terms for d in drugs_by_cui[c]])
	nong_cui_assoc = nong_clus.intersection(cui_assoc_drugs) 
	print('No of drug-cui pairs: ',gdf.shape[0])
	print('No of cuis: ',len(gcui_terms))
	print('No of in-class, approved: ',len(gclus))
	print('No of cui-assoc, approved, out-of-class: ',len(nong_cui_assoc))
	g_drug_str = get_drug_syns_list(gclus)
	nong_drug_str = get_drug_syns_list(nong_cui_assoc)
	outf = open(os.path.join(rdir,short_name+'_drug_strings_0224.txt'),'w')
	n = outf.write(short_name+"_syns=\'" + g_drug_str + "\'\n" + "non_"+short_name+"_syns=\'" + nong_drug_str + "\'")
	outf.close()
	return (g_drug_str,nong_drug_str)

print('NPY')
(npys,nonnpys) = char_single_gene('NPY')
print('CNR2')
(cs,ncs) = char_single_gene('CNR2')
print('CXCR3')
(cx3,ncx3) = char_single_gene('CXCR3')
print('CXCR5')
(cx5,ncx5) = char_single_gene('CXCR5')
print('chemokine')
(cxs,ncxs) = char_gene_list(chemokine_proteins,'chemokine','cxc')
print('complement system')
(cps,ncps) = char_gene_list(complement_proteins,'complement_system','cs')
