####################

'''
- ALS project
- Incorporating LINCS (meta)data
- Go to the corresponding directory
    #terminal: 'cd /als/lincs/'
- Activate the virtual environment in the directory
    #terminal: 'source alsproj/bin/activate'
- Run the code
    #terminal: 'python code_level51.py'
'''

####################

# libraries

import os
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
import cmapPy.pandasGEXpress.subset_gctoo as sg
import matplotlib.pyplot as plt
import seaborn as sns

# directory

data_dir = '/als/lincs/'
os.chdir(data_dir)

print('')
print('####################')
print('')

####################

# read cellinfo.txt: Metadata for cell lines

print('>>>>>>>>>>')
print('read cellinfo.txt: Metadata for cell lines')
print('>>>>>>>>>>')
print('')
cellinfo = pd.read_csv('cellinfo_beta.txt', sep='\t')
cell_lineage_ls = ['autonomic_ganglia', 'central_nervous_system']
primary_disease_ls = ['neuroblastoma', 'normal stem cell sample', 'brain cancer']
filtered_cellinfo = cellinfo[cellinfo['cell_lineage'].isin(cell_lineage_ls) | cellinfo['primary_disease'].isin(primary_disease_ls)]
cell_iname_cellinfo = filtered_cellinfo['cell_iname'].tolist()
#print('cell_iname_cellinfo:', cell_iname_cellinfo)
print('Done!')

#print('')
#print('####################')
#print('')

# read siginfo.txt: Metadata for level 5 signatures

#print('>>>>>>>>>>')
#print('read siginfo.txt: Metadata for level 5 signatures')
#print('>>>>>>>>>>')
#print('')
#sig_inf = pd.read_csv('siginfo_beta.txt', sep='\t', low_memory=False)
#siginf = sig_inf[['cell_mfc_name', 'pert_mfc_id', 'pert_id', 'sig_id', 'cell_iname', 'cmap_name']]
#filtered_siginf = siginf[siginf['cell_iname'].isin(cell_iname_cellinfo)]
#print('filtered_siginf:', filtered_siginf)
#print('')
#print('cell_iname / filtered_siginf:', list(set(filtered_siginf['cell_iname'].tolist())))

print('')
print('####################')
print('')

# read broad_lincs_sig_info: Metadata for each signature in the Level 5 matrix (metadata for the columns in the Level 5 data matrix)

print('>>>>>>>>>>')
print('read broad_lincs_sig_info.txt: Metadata for each signature in the Level 5 matrix (metadata for the columns in the Level 5 data matrix')
print('>>>>>>>>>>')
print('')
sig_info = pd.read_csv('broad_lincs_sig_info.txt', sep='\t', dtype=str)
siginfo = sig_info[['sig_id', 'pert_id', 'pert_iname', 'cell_id']]
filtered_siginfo_cell = siginfo[siginfo['cell_id'].isin(cell_iname_cellinfo)]
print('filtered_siginfo_cell:', filtered_siginfo_cell)
#print('')
#print('cell_id / filtered_siginfo_cell:', list(set(filtered_siginfo_cell['cell_id'].tolist())))
pert_iname_ls = siginfo['pert_iname'].tolist()

print('')
print('####################')
print('')

# read broad_lincs_gene_info: Metadata for each measured feature / gene (metadata for rows of the data matrices)

print('>>>>>>>>>>')
print('read broad_lincs_gene_info.txt: Metadata for each measured feature / gene (metadata for rows of the data matrices)')
print('>>>>>>>>>>')
print('')
geneinfo = pd.read_csv('broad_lincs_gene_info.txt', sep='\t', dtype=str)

complement_proteins = ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C3', 'C3AR1', 'C5', 'C5AR1']
print('# of complement proteins:', len(complement_proteins))
print('complement proteins:', (complement_proteins))
print('')
filtered_geneinfo = geneinfo[geneinfo['pr_gene_symbol'].isin(complement_proteins)]
print('filtered_geneinfo:', filtered_geneinfo)

print('')
print('####################')
print('')

# read cs_drug_strings: cs_drug_strings_0224.txt

print('>>>>>>>>>>')
print('read cs_drug_strings_0224.txt: cs_drug_strings')
print('>>>>>>>>>>')
print('')

with open('cs_drug_strings_0224.txt', 'r') as file:
    content = file.read()
cs_line = content.split('\n')[0]
cs_strings = cs_line.split('=')[1].strip(" '").split(',')
cs_strings_ls = [s.lower() for s in cs_strings]
#print('# of cs_syns drugs:', len(cs_strings_ls))
#print('cs_syns drugs:', cs_strings_ls)
non_cs_line = content.split('\n')[1]
non_cs_strings = non_cs_line.split('=')[1].strip(" '").split(',')
non_cs_strings_ls = [s.lower() for s in non_cs_strings]
#print('# of non_cs_syns drugs:', len(non_cs_strings_ls))
#print('non_cs_syns drugs:', non_cs_strings_ls)
print('Done!')

#print('')
#print('####################')
#print('')

# read the 'level5_beta_all.gctx' file: All Level 5 data

#print('>>>>>>>>>>')
#print('read the "level5_beta_all.gctx" file: All Level 5 data')
#print('>>>>>>>>>>')
#print('')

#my_level5_data_column = parse("level5_beta_all.gctx", cid=sample_column_ids)
#print('my_level5_data_column:', my_level5_data_column.data_df.shape)

#my_level5_data_row = parse("level5_beta_all.gctx", rid = gene_row_ids)
#print('my_level5_data_row:', my_level5_data_row.data_df.shape)

print('')
print('####################')
print('')

# read the 'level5_modz.gctx.gz' file: Level 5 data (moderated z-scores)

print('>>>>>>>>>>')
print('read the "level5_modz.gctx" file: Level 5 data (moderated z-scores)')
print('>>>>>>>>>>')
print('')

# >>>>>>>>>>>>>>>>>>>>
# rows

#gene_row_ids = filtered_geneinfo['pr_gene_id']
#print('# of gene_row_ids:', len(gene_row_ids))
#print('gene_row_ids:', gene_row_ids)
#print('')

#my_level5_modz_row = parse('level5_modz.gctx', rid = gene_row_ids)

#data_pr = my_level5_modz_row.data_df
#print('data_pr.shape:', data_pr.shape)
#print('data_pr.columns:', data_pr.columns)
#print('')
#df = my_level5_modz_row.data_df
#df.to_csv('parsed_complements.txt', sep='\t', index=False)

# >>>>>>>>>>>>>>>>>>>>
# columns

inter_pert_cs = [x for x in pert_iname_ls if x in cs_strings_ls]
pert_cs_ls = list(set(inter_pert_cs))
inter_pert_noncs = [x for x in pert_iname_ls if x in non_cs_strings_ls]
pert_noncs_ls = list(set(inter_pert_noncs))

rows_to_keep = [1496, 3221, 5334, 5975, 6699, 8689, 9077, 10694]
for gene_row_id in rows_to_keep:
    gene_expression_data_cs = []
    gene_expression_data_noncs = []

    for pert in pert_cs_ls + pert_noncs_ls:
        sample_column_ids = siginfo['sig_id'][siginfo['pert_iname'] == pert]

        if not sample_column_ids.empty:
            my_level5_modz_column = parse('level5_modz.gctx', cid=sample_column_ids)

            sample_sig_id_info = siginfo[siginfo['pert_iname'] == pert]
            sample_sig_id_info.set_index('sig_id', inplace=True)
            my_level5_modz_column.col_metadata_df = sample_sig_id_info

            cell_line_ls = my_level5_modz_column.col_metadata_df['cell_id'].unique()
            cell_iname_cellinfo = ['SHSY5Y', 'NEU', 'NPC', 'LN229', 'SKNSH', 'U251MG', 'YH13', 'GI1']
            inter_cell_lines = [x for x in cell_line_ls if x in cell_iname_cellinfo]

            for cell in inter_cell_lines:
                sample_cell_ids = my_level5_modz_column.col_metadata_df.index[my_level5_modz_column.col_metadata_df['cell_id'] == cell]
                if not sample_cell_ids.empty:
                    sample_cell_ids_list = sample_cell_ids.tolist()
                    sample_cell_gctoo = sg.subset_gctoo(my_level5_modz_column, cid=sample_cell_ids_list)

                    sample_cell_gene_data = sample_cell_gctoo.data_df.iloc[gene_row_id]

                    if pert in pert_cs_ls:
                        #gene_expression_data_cs.append(sample_cell_gene_data)
                        gene_expression_data_cs.extend(sample_cell_gene_data)
                    else:
                        #gene_expression_data_noncs.append(sample_cell_gene_data)
                        gene_expression_data_noncs.extend(sample_cell_gene_data)

    plt.hist(gene_expression_data_cs, bins=30, alpha=0.5, color='blue', label='cs_syns', density=True)
    plt.hist(gene_expression_data_noncs, bins=30, alpha=0.5, color='red', label='non_cs_syns', density=True)
    plt.xlabel('Gene Expression')
    plt.ylabel('Frequency')
    plt.legend(loc='best')

    plt.savefig(f'gene_{gene_row_id}_histogram.png')
    plt.clf()

print('')
print('####################')
print('')

####################
