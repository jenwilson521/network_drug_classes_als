# this code shows how we curated phenotypes
# from the PathFX database
# written 3-14-24 JLW

#### PATHFX pathways ####       
# ALS and other cognitive pathways from existing PathFX
p2cf = os.path.join(ddir,'Pfx050120_all_phens_to_cuis.pkl')
p2c = pickle.load(open(p2cf,'rb'))
c2gf = os.path.join(ddir,'Pfx050120_merged_unique_cuis2genes.pkl')
c2g = pickle.load(open(c2gf,'rb'))
c2pf = os.path.join(ddir,'Pfx050120_cui_to_phens.pkl')
c2p = pickle.load(open(c2pf,'rb'))

# modify this per disease
key_terms = [ 'Lou Gehrig','amyotrophic','sclerosis']
key_terms += ['alzheimer', 'ataxia', 'huntington', 'parkinson', 'palsy', 'motor neuron']
keep_phens = [k for k in p2c.keys() for kt in key_terms if kt.lower() in k.lower()] #733 phenotypes

# remove erroneous matches
err_terms = ['Atherosclerosis','Arteriosclerosis','glomerulosclerosis','Osteosclerosis','Scleroderma','Wolff-Parkinson-White','Tuberous sclerosis']
def check_err_terms(s):
        has_err_term = False
        for r in err_terms:
                if r.lower() in s.lower():
                        has_err_term = True
        return has_err_term

short_list = [p for p in keep_phens if not (check_err_terms(p))] # 669 phenotypes
keep_cuis = set([p2c[k] for k in short_list]) # 438 phenotypes
c25g = [c for c in keep_cuis if len(c2g[c]) >=25] # 60 phenotypes
pickle.dump(c25g,open('../data/pfx_0623_c25g.pkl','wb'))

# remove 443 phenotypes, random sample an additional ~1000 phenotypes to add
rem_cuis = set(c2g.keys()).difference(set(c25g)) # remaining cuis
rem_cuis_25 = [c for c in rem_cuis if len(c2g[c]) >=25]
rem_cuis_ran_1000 = np.random.choice(rem_cuis_25, size = 1000, replace=False)
