# looks at all drugs associated with a search
# phenotype, creates some summary tables
# re-written 9-8-23 JLW
# applying search to neurodegen. disease expansion

import pickle,os,csv
from collections import defaultdict

# first gather relevant phenotypes
# Move this to the call_search_phen script
p2c = pickle.load(open('../rscs/Pfx050120_all_phens_to_cuis.pkl','rb'))
c2g = pickle.load(open('../rscs/Pfx050120_merged_unique_cuis2genes.pkl','rb'))

# seed search with neuro key words
key_terms = {'alzheimer','parkinson','myasthenia','stroke','cerebral infarction', 'multiple sclerosis','gravis','amyotrophic'}
key_words_to_phens = defaultdict(set)
for k in p2c.keys():
	for kt in key_terms:
		if kt in k.lower():
			key_words_to_phens[kt].add(k)

# manually review and remove the following phenotypes
rem_phens = ['Autosomal Dominant Juvenile Parkinson Disease','Familial Juvenile Parkinsonism','Parkinsonism, Juvenile','Wolff-Parkinson-White pattern','Heat Stroke','Familial Alzheimer-like prion disease','Lewy Body Variant of Alzheimer Disease','Junctional epidermolysis bullosa gravis of Herlitz','Icterus Gravis Neonatorum']

keep_cuis = [p2c[k] for (kt,phen_set) in key_words_to_phens.items() for k in phen_set if k not in rem_phens]
print('Keep Phenotypes ',len(keep_cuis)) #220
keep_cuis_25 = [c for c in keep_cuis if len(c2g[c])>=25]
print('Phens with > 25 genes ',len(keep_cuis_25)) #45

aname = "neuro_expansion"
rdir = os.path.join("../results/",aname)
if not os.path.exists(rdir):
	os.makedirs(rdir)
cui_file =os.path.join(rdir,aname+"_keep_cuis.pkl") 
pickle.dump(keep_cuis,open(cui_file,'wb'))

cmd = "python search_DrugBank_for_phen_abbrev.py -aname %s -rdir %s -cuipkl %s"%(aname,rdir,aname+"_keep_cuis.pkl")
print(cmd)
os.system(cmd)


