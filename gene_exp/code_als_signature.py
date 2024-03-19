
####################

# libraries

import os
import pandas as pd
import pickle

# directory

data_dir = ''
os.chdir(data_dir)

print('####################')
print('')

####################

# analyze the drug signature data

print('Drug Signature Data:')
print('')

# PharmOmics_drug_signature_database

df_drug_sig = pd.read_csv('PharmOmics_drug_signature_database.txt', sep='\t')
print('All species.shape:', df_drug_sig.shape)
df_hs = df_drug_sig.loc[(df_drug_sig['species'] == 'Homo sapiens')]
print('Homo sapiens.shape:', df_hs.shape)
print('')
df_dr = df_hs['drug'].tolist()
df_dr_ls = [[x for x in word.lower().split()] for word in df_dr]
df_drug_ls = [''.join(x) for x in df_dr_ls]
print('# of drugs in the drug sig ls:', len(df_drug_ls))
print('')

db2n = pickle.load(open(('drugbankid_to_name.pkl'),'rb'))
drugbank_names = []
for name, value in db2n.items():
  drugbank_names.append(value.lower())
print('# of drugs in the drugbank:', len(drugbank_names))
print('')

sigdrug_ls = [x for x in df_drug_ls if x in drugbank_names]
print('# of mapped sig drugs in the drugbank:', len(sigdrug_ls))
print('')

print('####################')
print('')

####################

# analyze the ALS data: cs_drug_strings_0224.txt

print('CS ALS Data:')
print('')

# cs_drug_strings

with open('cs_drug_strings_0224.txt', 'r') as file:
    content = file.read()

cs_line = content.split('\n')[0]
cs_strings = cs_line.split('=')[1].strip(" '").split(',')
l1 = [s.lower() for s in cs_strings]
print('# of cs_syns drugs:', len(l1))
print('cs_syns drugs:', l1)
print('')

non_cs_line = content.split('\n')[1]
non_cs_strings = non_cs_line.split('=')[1].strip(" '").split(',')
l2 = [s.lower() for s in non_cs_strings]
print('# of non_cs_syns drugs:', len(l2))
print('non_cs_syns drugs:', l2)
print('')

print('####################')
print('')

####################

print('CS ALS & Drug Signature Data:')
print('')

# intersection

dr_ls1 = [x for x in l1 if x in sigdrug_ls]
print('# of cs_syns drugs in mapped sig drugs:', len(dr_ls1))
dr_ls2 = [x for x in l2 if x in sigdrug_ls]
print('# of non_cs_syns drugs in mapped sig drugs:', len(dr_ls2))
print('')

# drugs & up_&_down genes dataframes

map_dr_ls1 = [i.capitalize() for i in dr_ls1]
df1_data = df_hs.loc[df_hs['drug'].isin(map_dr_ls1)]
drug_gene1 = df1_data[['drug', 'genes_up', 'genes_down']]
print('drug_gene_cs.shape:', drug_gene1.shape)
map_dr_ls2 = [i.capitalize() for i in dr_ls2]
df2_data = df_hs.loc[df_hs['drug'].isin(map_dr_ls2)]
drug_gene2 = df2_data[['drug', 'genes_up', 'genes_down']]
print('drug_gene_non_cs.shape:', drug_gene2.shape)
print('')

# up_&down gene lists

genes_up_1 = drug_gene1['genes_up'].tolist()
genes_up1 = [x for x in genes_up_1 if pd.notnull(x)] # nan is not a string
genes_up1_all = []
for i in range(len(genes_up1)):
  agl = (genes_up1[i].split(','))
  genes_up1_all.extend(agl)
print('# of all genes_up in the cs_syns data:', len(genes_up1_all))
genes_up1_ls = list(set(genes_up1_all))
print('# of distinct genes_up in the cs_syns data:', len(genes_up1_ls))
print('')

genes_down_1 = drug_gene1['genes_down'].tolist()
genes_down1 = [x for x in genes_down_1 if pd.notnull(x)] # nan is not a string
genes_down1_all = []
for i in range(len(genes_down1)):
  agl = (genes_down1[i].split(','))
  genes_down1_all.extend(agl)
print('# of all genes_down in the cs_syns data:', len(genes_down1_all))
genes_down1_ls = list(set(genes_down1_all))
print('# of distinct genes_down in the cs_syns data:', len(genes_down1_ls))
print('')

genes_up_2 = drug_gene2['genes_up'].tolist()
genes_up2 = [x for x in genes_up_2 if pd.notnull(x)] # nan is not a string
genes_up2_all = []
for i in range(len(genes_up2)):
  ag2 = (genes_up2[i].split(','))
  genes_up2_all.extend(ag2)
print('# of all genes_up in the non_cs_syns data:', len(genes_up2_all))
genes_up2_ls = list(set(genes_up2_all))
print('# of distinct genes_up in the non_cs_syns data:', len(genes_up2_ls))
print('')

genes_down_2 = drug_gene2['genes_down'].tolist()
genes_down2 = [x for x in genes_down_2 if pd.notnull(x)] # nan is not a string
genes_down2_all = []
for i in range(len(genes_down2)):
  ag2 = (genes_down2[i].split(','))
  genes_down2_all.extend(ag2)
print('# of all genes_down in the non_cs_syns data:', len(genes_down2_all))
genes_down2_ls = list(set(genes_down2_all))
print('# of distinct genes_down in the non_cs_syns data:', len(genes_down2_ls))
print('')

# check if complement proteins are in different gene lists

complement_proteins = ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C3', 'C3AR1', 'C5', 'C5AR1']
print('# of complement proteins:', len(complement_proteins))
print('complement proteins:', (complement_proteins))
print('')

print('####################')
print('')

####################

# results: (# of) complement proteins in different lists

print('Results: (# of) complement proteins in different lists:')
print('')

gene_cp_up1 = [x for x in complement_proteins if x in genes_up1_ls]
print('# of complement proteins in genes_up cs_syns:', len(gene_cp_up1))
print('complement proteins in genes_up cs_syns:', (gene_cp_up1))
print('')

gene_cp_down1 = [x for x in complement_proteins if x in genes_down1_ls]
print('# of complement proteins in genes_down cs_syns:', len(gene_cp_down1))
print('complement proteins in genes_down cs_syns:', (gene_cp_down1))
print('')

gene_cp_up2 = [x for x in complement_proteins if x in genes_up2_ls]
print('# of complement proteins in genes_up non_cs_syns:', len(gene_cp_up2))
print('complement proteins in genes_up non_cs_syns:', (gene_cp_up2))
print('')

gene_cp_down2 = [x for x in complement_proteins if x in genes_down2_ls]
print('# of complement proteins in genes_down non_cs_syns:', len(gene_cp_down2))
print('complement proteins in genes_down non_cs_syns:', (gene_cp_down2))
print('')

print('####################')
print('')

####################

# results: (# of) drugs associated with complement proteins in different lists

print('Results: (# of) drugs associated with complement proteins in different lists:')
print('')

def gene2drug (df, target_gene):

  import pandas as pd

  # create mappings of genes to drugs
  gene_to_drug_map_up = {}
  gene_to_drug_map_down = {}
  # iterate through each row in the DataFrame
  for index, row in df.iterrows():
    drug = row['drug']
    # process genes_up&down columns
    if not pd.isna(row['genes_up']):
        genes_up = row['genes_up'].split(',')
        for gene in genes_up:
            gene = gene.strip()
            if gene not in gene_to_drug_map_up:
                gene_to_drug_map_up[gene] = []
            gene_to_drug_map_up[gene].append(drug)
    if not pd.isna(row['genes_down']):
        genes_down = row['genes_down'].split(',')
        for gene in genes_down:
            gene = gene.strip()
            if gene not in gene_to_drug_map_down:
                gene_to_drug_map_down[gene] = []
            gene_to_drug_map_down[gene].append(drug)
  # query the gene_to_drug_map_up&down to find drugs associated with a specific gene
  associated_drugs_up = gene_to_drug_map_up.get(target_gene, [])
  associated_drugs_down = gene_to_drug_map_down.get(target_gene, [])

  return associated_drugs_up, associated_drugs_down

print('##########')
print('')

df1 = drug_gene1
df2 = drug_gene2
complement_proteins = ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C3', 'C3AR1', 'C5', 'C5AR1']

for target_gene in complement_proteins:
  associated_drugs_up1, associated_drugs_down1 = gene2drug (df1, target_gene)
  associated_drugs_up2, associated_drugs_down2 = gene2drug (df2, target_gene)
  print('# of drugs associated with %s in gene_up cs_syns:' %target_gene, len(associated_drugs_up1))
  print('drugs associated with %s in gene_up cs_syns:' %target_gene, (associated_drugs_up1))
  print('')
  print('# of drugs associated with %s in gene_down cs_syns:' %target_gene, len(associated_drugs_down1))
  print('drugs associated with %s in gene_down cs_syns:' %target_gene, (associated_drugs_down1))
  print('')
  print('# of drugs associated with %s in gene_up non_cs_syns:' %target_gene, len(associated_drugs_up2))
  print('drugs associated with %s in gene_up non_cs_syns:' %target_gene, (associated_drugs_up2))
  print('')
  print('# of drugs associated with %s in gene_down non_cs_syns:' %target_gene, len(associated_drugs_down2))
  print('drugs associated with %s in gene_down non_cs_syns:' %target_gene, (associated_drugs_down2))
  print('')
  print('##########')
  print('')

print('####################')
print('')

####################
