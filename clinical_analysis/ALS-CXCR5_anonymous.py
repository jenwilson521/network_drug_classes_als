# This is meant to demonstrate the commands used, assuming
# access to Market Clarity data in Pyspark tables.
# the code can easily be modified with the provided diagnosis
# codes and drug_strings for other network classes
# written 3-13-23 JLW

# DBTITLE 1,PySpark Import Statements
from pyspark.sql.functions import col, when, count, min, date_sub, concat_ws, substring, lit, coalesce, collect_list, countDistinct, expr, datediff, max, row_number
from pyspark.sql.window import Window
import matplotlib.pyplot as plt 
import pandas as pd
from pyspark.sql.functions import sum as _sum
DRUG_CLASS = 'CXCR5'

# DBTITLE 1,Previous Analysis, creating filtered tables, only run once
## PRE-Analysis: load diagnosis table and get ALS patient IDs, write able
#mdf = spark.read.format("delta").load('Medical Diagnosis')

## extract list of patient ids for selecting other features
#mdf_als = mdf.filter(col('DIAG')==33520) # patients with diagnosis, contains repeats per patient

## convert patient IDs to list
#patids_list = mdf_als.dropDuplicates(["PATID",]).select(["PATID",]).collect()
#als_patids = [p.PATID for p in patids_list]

## save the filtered diagnosis table 
#mdf_als.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_als_diagnosis.csv")

## PRE-Analysis: re-load diagnosis and save codes for als_patids
#mds_other_als = mdf.filter(col('PATID').isin(als_patids))
#mds_other_als.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_mds_other_als.csv")

## PRE-Analysis: load RX claims data, save RX values for patients above
#rdf = spark.read.format("delta").load('RX Claims')
#rdf_als = rdf.filter(col('PATID').isin(als_patids))
#rdf_als.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_als_rx_claims.csv")

## diagnostics other than ALS for covariate analysis
#diagdf = spark.read.format("delta").load('/Medical Diagnosis/')
#diagdf_als = diagdf.filter(col('PATID').isin(als_patids))
#diagdf_als.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_als_med_diag.csv")

## procedures just in case, not currently in propensity model
#procdf = spark.read.format("delta").load('/Medical Procedures/')
#proc_als = procdf.filter(col('PATID').isin(als_patids))
#proc_als.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_als_med_proc.csv")


### COPY/PASTE contents of relevant drug_strings file here
CXCR5_syns=''
non_CXCR5_syns=''

# DBTITLE 1,Load filtered tables
# read a saved, filtered table - use this block to initialize environment before getting started
mdf_als = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_als_diagnosis.csv")
rdf_als = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_als_rx_claims.csv")
diagdf_als = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_als_med_diag.csv")
# proc_als = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_als_med_proc.csv")
mds_other_als = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_mds_other_als.csv")

# get a list of patient IDs, convert to list
patids_list = mdf_als.dropDuplicates(["PATID",]).select(["PATID",]).collect()
als_patids = [p.PATID for p in patids_list]

# DBTITLE 1,Load outcomes, demographic tables, look-up for diagnosis codes
# load outcomes and demographic information tables
doddf_als = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_als_dod.csv")
lud = spark.read.format("delta").load('/Lookup Diagnosis/')
memdf_als = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_member_info_als.csv")

# DBTITLE 1,Get initial treatment/comparator cohorts
# MODIFY DRUG STRINGS AND SAVE NAME
tdrug_str = CXCR5_syns # target
cdrug_str = non_CXCR5_syns # comparator
DRUG_CLASS = 'CXCR5'
col_list = ['PATID','BRND_NM','DAYS_SUP','FILL_DT','GNRC_IND']

# use the target and comparator drug strings
tdrug_list = [d.upper() for d in tdrug_str.split(',')]
trx_claims = rdf_als.filter((rdf_als.BRND_NM.isin(tdrug_list))).select(col_list)
trx_data = trx_claims.collect()
target_patids = set([r.PATID for r in trx_data])

cdrug_list = [d.upper() for d in cdrug_str.split(',')]
crx_claims = rdf_als.filter((rdf_als.BRND_NM.isin(cdrug_list))).select(col_list)
crx_data = crx_claims.collect()
comparator_patids = set([r.PATID for r in crx_data])

# look for shared patient IDs, remove duplicates
both_patids = target_patids.intersection(comparator_patids)
final_target_patids = target_patids.difference(both_patids) # remove those exposed to both
final_comparator_patids = comparator_patids.difference(both_patids)

# one or more target drugs, and als diagnosis, and removed if also in comparator group
target_rx_claims = trx_claims.filter(trx_claims.PATID.isin(final_target_patids))
target_rx_claims.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_target_rx_claims.csv")
comparator_rx_claims = crx_claims.filter(crx_claims.PATID.isin(final_comparator_patids))
comparator_rx_claims.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_rx_claims.csv")

# COMMAND ----------

# DBTITLE 1,Cohort characterization - number of prescriptions/exposures
# count exposure days and plot
DRUG_CLASS = 'CXCR5'
target_rx_claims = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_target_rx_claims.csv")
comparator_rx_claims = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_rx_claims.csv")

target_rx_claims = target_rx_claims.na.fill(1, subset=["DAYS_SUP"]) # use a dummy count of 1 when no days are provided
target_rx_claims = target_rx_claims.withColumn("DAYS_SUP_INT", col("DAYS_SUP").cast("integer"))
drug_sum_df = target_rx_claims.groupBy("PATID", "BRND_NM").agg(_sum("DAYS_SUP_INT").alias("total_days_supply")) # total days per patient
drug_count_df = target_rx_claims.groupBy("PATID", "BRND_NM").agg(count("*").alias("COUNT")) # total count per patient
# display(drug_sum_df)

targ_df = drug_sum_df.toPandas()
targ_data_series = targ_df.groupby("BRND_NM")["total_days_supply"].apply(list)
targ_sorted_drug_names = targ_data_series.keys()

# Create a boxplot per drug
fig_a, ax_a = plt.subplots(1, 2, figsize=(30, 6))
for (i,days_list) in enumerate(targ_data_series):
  x = [i+1] * len(days_list)
  ax_a[0].scatter(x,days_list)
ax_a[0].set_xlabel("Drug")
ax_a[0].set_ylabel("Exposure Days")
ax_a[0].set_title("Total Exposure Days per Drug\n TARGET, ALS, "+DRUG_CLASS)
ax_a[0].set_xticks(range(1,len(targ_sorted_drug_names)+1))
ax_a[0].set_xticklabels(targ_sorted_drug_names,rotation=90)

comparator_rx_claims = comparator_rx_claims.withColumn("DAYS_SUP_INT", col("DAYS_SUP").cast("integer"))
drug_sum_df_comp = comparator_rx_claims.groupBy("PATID", "BRND_NM").agg(_sum("DAYS_SUP_INT").alias("total_days_supply"))
comp_df = drug_sum_df_comp.toPandas()
comp_data_series = comp_df.groupby("BRND_NM")["total_days_supply"].apply(list)
comp_sorted_drug_names = comp_data_series.keys()

for (i,days_list) in enumerate(comp_data_series):
  x = [i+1] * len(days_list)
  ax_a[1].scatter(x,days_list)
ax_a[1].set_xlabel("Drug")
ax_a[1].set_ylabel("Exposure Days")
ax_a[1].set_title("Total Exposure Days per Drug\n COMPARATOR, ALS, "+DRUG_CLASS)
ax_a[1].set_xticks(range(1,len(comp_sorted_drug_names)+1))
ax_a[1].set_xticklabels(comp_sorted_drug_names,rotation=90)

plt.subplots_adjust(bottom=0.25)
plt.show()

print('Number of brand names in target group', len(targ_sorted_drug_names))
print('Number of brand names in comparator group', len(comp_sorted_drug_names))

tunique_patient_count = target_rx_claims.select(countDistinct("PATID").alias("unique_patient_count")).first()
count_t = tunique_patient_count["unique_patient_count"]
cunique_patient_count = comparator_rx_claims.select(countDistinct("PATID").alias("unique_patient_count")).first()
count_c = cunique_patient_count["unique_patient_count"]

print('Number of unique patients in target group', count_t)
print('Number of patients in the comparator group', count_c)

# COMMAND ----------

# DBTITLE 1,Consider experimental validation for drugs with most patients
# look at total unique patient and total exposure days across the cohort
drug_patient_sum_df_target = target_rx_claims.groupBy("BRND_NM").agg(countDistinct('PATID').alias("Unique_Patients"))
drug_patient_sum_df_target.display()
drug_patient_sum_df_comp = comparator_rx_claims.groupBy("BRND_NM").agg(countDistinct('PATID').alias("Unique_Patients"))
drug_patient_sum_df_comp.display()

drug_expDay_sum_df_target = target_rx_claims.groupBy("BRND_NM").agg(_sum("DAYS_SUP_INT").alias("total_days_supply"))
drug_expDay_sum_df_target.display()
drug_expDay_sum_df_comp = comparator_rx_claims.groupBy("BRND_NM").agg(_sum("DAYS_SUP_INT").alias("total_days_supply"))
drug_expDay_sum_df_comp.display()

# COMMAND ----------

# DBTITLE 1,Find the earliest prescription format DOD to be first of the month, save treatment effect
target_rx_claims = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_target_rx_claims.csv")
comparator_rx_claims = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_rx_claims.csv")
# join with DOD, format, calculate "min date" as first exposure to either drug class
col_list = ['PATID','dod_YMDOD','DAYS_SUP','FILL_DT'] # keep all of the rx rows, to find first cohort entry
doddf_als = doddf_als.select([col(c).alias("dod_" + c) for c in doddf_als.columns])
target_with_outcome = doddf_als.join(target_rx_claims, doddf_als["dod_PATID"] == target_rx_claims["PATID"], "right").select(col_list) 
comparator_with_outcome = doddf_als.join(comparator_rx_claims, doddf_als["dod_PATID"] == comparator_rx_claims["PATID"], "right").select(col_list)

# Define a window specification partitioned by patient ID and ordered by date
from pyspark.sql.functions import col, min
from pyspark.sql.window import Window
from pyspark.sql.functions import concat_ws, col, substring, lit
window_spec = Window.partitionBy("PATID").orderBy("FILL_DT")

# Add a column with the minimum date per patient, filter where fill date matches earliest date, modify dod column to be the first of the month
target_add_earliest_rx = target_with_outcome.withColumn("min_date", min(col("FILL_DT")).over(window_spec))
target_earliest_rx_only = target_add_earliest_rx.filter(col("FILL_DT") == col("min_date"))
target_earliest_rx_only = target_earliest_rx_only.withColumn("dod_year_month", concat_ws("-", substring(col("dod_YMDOD"), 1, 4), substring(col("dod_YMDOD"), 5, 2), lit("01")))

comp_add_earliest_rx = comparator_with_outcome.withColumn("min_date", min(col("FILL_DT")).over(window_spec))
comp_earliest_rx_only = comp_add_earliest_rx.filter(col("FILL_DT") == col("min_date"))
comp_earliest_rx_only = comp_earliest_rx_only.withColumn("dod_year_month", concat_ws("-", substring(col("dod_YMDOD"), 1, 4), substring(col("dod_YMDOD"), 5, 2), lit("01")))

from pyspark.sql.functions import datediff
target_earliest_rx_only = target_earliest_rx_only.withColumn("date_diff", datediff(target_earliest_rx_only.dod_year_month, target_earliest_rx_only.FILL_DT))
comp_earliest_rx_only = comp_earliest_rx_only.withColumn("date_diff", datediff(comp_earliest_rx_only.dod_year_month, comp_earliest_rx_only.FILL_DT))
display(comp_earliest_rx_only)

# add treatment effect and combine tables
target_earliest_rx_only = target_earliest_rx_only.withColumn("treatment", lit(1))
comp_earliest_rx_only = comp_earliest_rx_only.withColumn("treatment", lit(0))
pts_w_tr_ocs = target_earliest_rx_only.union(comp_earliest_rx_only) # patients with treatments and outcomes

target_earliest_rx_only.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_target_earliest_rx_only.csv")
comp_earliest_rx_only.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_comp_earliest_rx_only.csv")
pts_w_tr_ocs.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_pts_w_tr_ocs.csv")

# COMMAND ----------

# DBTITLE 1,Find all non-ALS diagnoses, keep the most frequent codes
# join on diagnosis covariates w/in -365, and then count, save the top-most frequent
pts_w_tr_ocs = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_pts_w_tr_ocs.csv")
# add a prior date for covariate joining
pts_w_tr_ocs = pts_w_tr_ocs.withColumn("prior_date", date_sub(pts_w_tr_ocs.min_date,365))

# add aliases to make merging easier
diagdf_als_alias = diagdf_als.select([col(c).alias("diag_" + c) for c in diagdf_als.columns])
diagdf_als_alias = diagdf_als_alias.filter(col("diag_DIAG") != 33520)

# first get diagnoses witin window, skip ALS
pts_w_tr_ocs_with_all_diag = pts_w_tr_ocs.join(diagdf_als_alias, 
  (pts_w_tr_ocs.PATID == diagdf_als_alias.diag_PATID) & 
  (diagdf_als_alias.diag_FST_DT >= pts_w_tr_ocs.prior_date) &
  (diagdf_als_alias.diag_FST_DT < pts_w_tr_ocs.min_date), "left")

diag_col = ['PATID','min_date','treatment','dod_year_month','date_diff','prior_date','diag_DIAG','diag_FST_DT']
pts_w_tr_ocs_with_all_diag = pts_w_tr_ocs_with_all_diag.select(diag_col)

# assess number of patients associated with each code
diag_count_df = pts_w_tr_ocs_with_all_diag.groupBy('diag_DIAG').agg(count("*").alias("COUNT")) # total rows (patients) per code
filtered_df = diag_count_df.filter(col("COUNT") > 100) #904 rows
filtered_df = filtered_df.join(lud, filtered_df.diag_DIAG == lud.DIAG_CD,'left').select(['diag_DIAG','COUNT','DIAG_DESC', 'DIAG_FST3_DESC']).distinct()
code_count = filtered_df.select("diag_DIAG").distinct().count()
display(filtered_df) # this previously had duplicated entries for diagnosis codes because multiple mappings
filtered_df.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_allOtherDiag_counts.csv")

# filter patient dataframe by highest counted diagnosis codes and save
filtered_df_alias = filtered_df.select([col(c).alias("fd_" + c) for c in filtered_df.columns])
pts_w_tr_ocs_with_all_diag_top = pts_w_tr_ocs_with_all_diag.join(filtered_df_alias.select('fd_diag_DIAG'), pts_w_tr_ocs_with_all_diag.diag_DIAG == filtered_df_alias.fd_diag_DIAG, "inner")

# display(pts_w_tr_ocs_with_all_diag_top.distinct())

# COMMAND ----------

# DBTITLE 1,Find all non-class drug prescriptions, keep the most frequent
# join on rx coviariates in -365 days prior, and then count
target_rx_claims = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_target_rx_claims.csv")
comparator_rx_claims = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_rx_claims.csv")

rdf_als_alias = rdf_als.select([col(c).alias("rx_" + c) for c in rdf_als.columns])
all_rx_claims = target_rx_claims.union(comparator_rx_claims) # all study drugs

# Join the tables and filter out the drug claims that are not "target" or "comparator"
rdf_als_filtered = rdf_als_alias.join(all_rx_claims, rdf_als_alias.rx_BRND_NM == all_rx_claims.BRND_NM, "left_anti")

# now join the table and filter for -365 days
pts_w_tr_ocs_with_all_rx = pts_w_tr_ocs.join(rdf_als_filtered, 
  (pts_w_tr_ocs.PATID == rdf_als_filtered.rx_PATID) & 
  (rdf_als_filtered.rx_FILL_DT >= pts_w_tr_ocs.prior_date) &
  (rdf_als_filtered.rx_FILL_DT < pts_w_tr_ocs.min_date), "left")

rx_col = ['PATID','min_date','treatment','dod_year_month','date_diff','prior_date','rx_BRND_NM','rx_FILL_DT']
pts_w_tr_ocs_with_all_rx = pts_w_tr_ocs_with_all_rx.select(rx_col)


# assess number of patients associated with each code
rx_count_df = pts_w_tr_ocs_with_all_rx.groupBy('rx_BRND_NM').agg(count("*").alias("COUNT"), count("PATID").alias("PATCOUNT")) # total patient-drug claims
rx_filtered_df = rx_count_df.filter(col("COUNT") > 100) # 
code_count = rx_filtered_df.select("rx_BRND_NM").distinct().count()
display(rx_filtered_df)

# filter patient dataframe by highest counted diagnosis codes and save
rx_filtered_df_alias = rx_filtered_df.select([col(c).alias("fd_" + c) for c in rx_filtered_df.columns])
pts_w_tr_ocs_with_all_rx_top = pts_w_tr_ocs_with_all_rx.join(rx_filtered_df_alias.select('fd_rx_BRND_NM'), pts_w_tr_ocs_with_all_rx.rx_BRND_NM == rx_filtered_df_alias.fd_rx_BRND_NM, "inner")
display(pts_w_tr_ocs_with_all_rx_top)

# COMMAND ----------

# DBTITLE 1,Use RX and diagnoses as covariates, create pivot table, add demographic information
# now take the union of covariates and expand for logisitic regression
pat_diag_dis = pts_w_tr_ocs_with_all_diag_top.select("PATID", "diag_DIAG").distinct().withColumnRenamed("diag_DIAG", "covar")
pat_rx_dis = pts_w_tr_ocs_with_all_rx_top.select("PATID",'rx_BRND_NM').distinct().withColumnRenamed("rx_BRND_NM","covar")
all_covar = pat_diag_dis.union(pat_rx_dis)

# Pivot the column with covariate IDs
pivoted_df = all_covar.groupBy("PATID").pivot("covar").agg(lit(1)).fillna(0)

# create alias and select relevant tables
mem_alias = memdf_als.select([col(c).alias("mem_" + c) for c in memdf_als.columns])
mem_col = ['mem_GDR_CD','mem_RACE','mem_YRDOB']
diag_rx_cov_cols = pivoted_df.columns
cov_cols = diag_rx_cov_cols + mem_col
covdf = pivoted_df.join(mem_alias,(pivoted_df.PATID == mem_alias.mem_PATID),'left').select(cov_cols)

# put the treatments values back
pts_w_tr_ocs_alias = pts_w_tr_ocs.select([col(c).alias("tr_" + c) for c in pts_w_tr_ocs.columns])
logr_col = covdf.columns + ["tr_treatment"]
covdf = covdf.join(pts_w_tr_ocs_alias,covdf.PATID == pts_w_tr_ocs_alias.tr_PATID).select(logr_col)
# display(covdf) # will have duplicates because of patients having multiple entries for the same diagnosis code (separate visits), or same RX (because refills)

covdf_distinct = covdf.distinct()
display(covdf_distinct)
# save the covariates data frame
covdf_distinct.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_covdf.csv")

# COMMAND ----------

# DBTITLE 1,Characterize demographic information
covdf.groupBy('mem_RACE').count().show()
covdf.groupBy('mem_GDR_CD').count().show()

# COMMAND ----------

# DBTITLE 1,Logistic Regression for Propensity Scores
# implement logistic regression, look at covariate balance, save propensity scores
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import mean_squared_error, r2_score, accuracy_score
import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np

covdf = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_covdf.csv")
# display(covdf)
df = covdf.toPandas()
cov_cols = [c for c in df.columns if c != 'tr_treatment' and c != 'PATID']

sex_dic = {'M':1,'F':0}
race_dic = {'B':1,'A':2,'W':3,'H':4}

X = df[cov_cols].fillna(0)
X['mem_RACE'] = X['mem_RACE'].map(race_dic)
X['mem_RACE'] = X['mem_RACE'].fillna(0) # fill those not recorded
X['mem_GDR_CD'] = X['mem_GDR_CD'].map(sex_dic)
y = df['tr_treatment']

psmodel = LogisticRegression(random_state=5)
psmodel.fit(X, y)
predicted_probs = psmodel.predict_proba(X)
predicted_classes = psmodel.predict(X)
prob_treat = [x[0] for x in predicted_probs]
df['propensity_score'] = prob_treat

# assessment, eventually convert and save covariates
print('Mean squared error: %.2f'% mean_squared_error(y, predicted_classes))
print('Coefficient of determination: %.2f'% r2_score(y, predicted_classes))

# TO-DO Map Coefficients and get the top values
reg_coef = list(zip(cov_cols, psmodel.coef_[0]))
#reg_coef_named = [(concepts_to_names[cv],mc) for (cv,mc) in reg_coef]
#sort_rc = sorted(reg_coef_named,key = lambda x:x[1],reverse=True)
sort_rc = sorted(reg_coef,key = lambda x:x[1],reverse=True)
# print(sort_rc[0:20])

# merge with outcomes data, save for survival analysis
pts_w_tr_ocs = spark.read.option("header","true").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_pts_w_tr_ocs.csv")
pts_w_tr_ocs = pts_w_tr_ocs.withColumn("event_occurred", when(col("date_diff").isNull(), 0).otherwise(1))
pdf = pts_w_tr_ocs.toPandas()
pdf['date_diff'] = pdf['date_diff'].astype(float).fillna(9999)
pdf_short = pdf[['PATID','date_diff']]

fulldf = pd.merge(df, pdf_short, on='PATID', how="left")
fulldf = fulldf[['PATID','tr_treatment','date_diff','propensity_score']]
fulldf['censored'] = np.where(fulldf['date_diff']== 9999, 1, 0)
display(fulldf)
fulldf_spark = spark.createDataFrame(fulldf)
fulldf_spark.write.format("csv").option("header","true").mode("OverWrite").csv("/dir_name/sub_dir_name/202302_"+DRUG_CLASS+"_als_survival_df.csv")

# COMMAND ----------

# DBTITLE 1,Look at top/bottom regression coefficients
reg_coef = list(zip(cov_cols, psmodel.coef_[0]))
rc_df = pd.DataFrame(reg_coef, columns =['Code', 'Coefficient'])
rcdf_spark = spark.createDataFrame(rc_df)
rcdf_with_desc = rcdf_spark.join(lud, rcdf_spark.Code == lud.DIAG_CD,'left').select(['Code', 'Coefficient','DIAG_DESC']).distinct()
display(rcdf_with_desc)


# COMMAND ----------

# DBTITLE 1,Plot propensity scores for treatment/comparator groups
fig,ax = plt.subplots(1,2)
fulldf['propensity_score'].hist(by=fulldf['tr_treatment'],ax=ax)
ax[0].set_xlabel("prop score")
ax[0].set_ylabel("num patients")
ax[1].set_xlabel("prop score")
ax[0].set_title("ALS, " + DRUG_CLASS + "\ncomparator")
ax[1].set_title("ALS, " + DRUG_CLASS + "\ntarget")

# COMMAND ----------

# DBTITLE 1,IPW Cox-Proportional Hazards
# MAGIC %r # re-read the PySpark table into R as a dataframe
# MAGIC library(survival)
# MAGIC library(ggplot2)
# MAGIC library(dplyr)
# MAGIC library(SparkR)
# MAGIC
# MAGIC DRUG_CLASS <- "CXCR5"
# MAGIC sparkR.session()
# MAGIC fullfilepath <- paste("/dir_name/sub_dir_name/202302_",DRUG_CLASS,"_als_survival_df.csv",sep="")
# MAGIC sdata = read.df(fullfilepath,"csv", sep = ",", inferSchema = TRUE, header = TRUE)
# MAGIC sdf = as.data.frame(sdata)
# MAGIC
# MAGIC # https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
# MAGIC km_fit <- survfit(Surv(date_diff, !censored) ~ tr_treatment, data=sdf)
# MAGIC summary(km_fit, times = c(1,30,60,90*(1:10)))
# MAGIC plot(km_fit, xlab="Days", main = paste(DRUG_CLASS,'\nKaplan Meyer Plot'), conf.int = TRUE, lty=1:2)
# MAGIC legend(5000, .9, c(paste("non", DRUG_CLASS, "drugs"), paste(DRUG_CLASS, "drugs")), lty = 1:2)
# MAGIC
# MAGIC # implement weighting
# MAGIC weights <- ifelse(sdf$tr_treatment==1, 1/sdf$propensity_score , 1/(1-sdf$propensity_score))
# MAGIC HR5 <- coxph(Surv(date_diff, !censored)~as.factor(tr_treatment), weights = weights, data = sdf) 
# MAGIC summary(HR5)
