# network_drug_classes_als

## Network Analysis
Example code for the network analysis is provided in the network_analysis folder.

### Finding PathFX phenotypes to incorporate into ALS-PathFX
We demonstrate how we leverage the published PathFX database to prioritize phenotypes to keep in the ALS-PathFX version in search_pathfx_database.py.

### Neurodegeneration disease indication expansion
Our code for creating the network heatmap images relies on a file Drugbank050120.xlsx that we cannot release. We recommend downloading DrugBank version 5.1.6. and we used the following code to parse the XML data base and create a .xlsx object that contains DrugBankID, DrugName, ATC codes, and description: ttps://github.com/dhimmel/drugbank/blob/gh-pages/parse.ipynb

With this version of DrugBank you can re-run plot_network_heatmap.py to create the heatmaps.

We provide examples for conducting the neuro-indication expansion analysis. This requires that you have a working installation of PathFX and PathFX networks for all drugs in DrugBank. You can clone the repository from here: https://github.com/jenwilson521/PathFX
And then move the following scripts into the 'PathFX/scripts/' directory to recreate these results:
run_PathFX_allDrugBank.py
search_neuro_expansion.py
read_neuro_expansion.py
### Discovery of ALS-associatied drug networks, generation of protein classes
We cannot release ALS-PathFX, but have provided code that can replicate most of these results, except assessing connections to sequencing-data pathways. We used the script, search_DrugBank_for_phen_v4.py, like above for the neuro-indication-expansion analysis, except that the search was run on all DrugBank networks created with ALS-PathFX instead of PathFXv2. Assuming that you have the PathFX repository and DrugBank files available you can re-create this analysis by moving the script into /PathFX/scripts/ and by changing the file paths in line 28 to reflect ../results/alldrugbank/ or whichever path you used. You may also comment out lines 64-69 which contain severe adverse events from drug's labels. If you are interested in that data, it is available through our other repository: https://github.com/jenwilson521/Designated-Medical-Event-Pathways

The results of this script are provided as Supplementary File 1 and are also included in this repository. 

To generate the drug-network classes, we used search_clusters_of_interst_v5_fixStrings.py. Again, this file requires the same DrugBank files from above, and the result from the prior script, all_drug_summary_als_san_0623.xlsx. You will need to update the path to the file in line 40. This will create supplementary file 2, all_net_prot_clusters_0224.xlsx. This script will also create slices of the above dataframe per network class and will generate .txt files of all drug strings and their synonyms which were needed for later clinical analysis. All of the drug strings files are included in drug_strings_0224.

## Gene expression analysis
We provided our code for accessing PharmOmics (code_als_signature.py). This script requires that you have access to the data as an object "PharmOmics_drug_signature_database.txt", the list of drug strings from the network analysis ("cs_drug_strings_0224.txt"), and a pickled dictionary of DrugBank ID to DrugName generated from DrugBank version 5.1.6 ("drugbankid_to_name.pkl").

We also provided code for accessing the LINCS dataset (code_level51.py). This requires that you have successfully downloaded the LINCS data and meta-data files: cellinfo.txt, siginfo.txt, read broad_lincs_sig_info.txt, read broad_lincs_gene_info.txt, level5_beta_all.gctx, and level5_modz.gctx.gz. This also requires the list of drug strings from the network analysis ("cs_drug_strings_0224.txt").

## Clinical data analysis
Example code for the clinical data analysis is provided in the clinical_analysis folder.

Again, we are unable to make clinical data available, but we accessed the Optum Market Clarity data through DataBricks. In this environment, the data were in PySpark tables and we provide anonymized code to demonstrate our pipeline and our use of the python sklearn and R survival packages. For users with access to Market Clarity, they can update the empty file paths and use the PySpark, python, and R code provided. We provide an example for the ALS patients, measuring effects of the CXCR5 class, but this notebook can easily be changed using the drug_strings (described above) for other network classes and the diagnosis codes provided in Supplementary File 4 to recreate all other analyses.
