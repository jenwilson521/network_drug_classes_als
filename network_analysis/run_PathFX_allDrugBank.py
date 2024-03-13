# re-written 3-19-21 JLW
import pickle,os,csv

analysis_name = 'allDrugBank'
dint = pickle.load(open('../rscs/drugs_to_drugbank_targets.pkl','rb')) # use symbolic links to get drug targets in each version
all_dbids = [k for k in dint.keys() if k[0:2]=='DB']
for drug_name in all_dbids:
	### call the algorithm without phenotype clustering
	cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
	os.system(cmd)

