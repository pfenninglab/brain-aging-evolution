import numpy as np
import pandas as pd
import anndata as ad
import scipy
import scanpy as sc
import pickle
import gc
import random

from sklearn.linear_model import ElasticNetCV
#from sklearn.linear_model import ElasticNet
#from sklearn.model_selection import sc_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_absolute_error
import scipy.stats as stats
from scipy.stats import pearsonr
from scipy.sparse import coo_matrix, hstack
from matplotlib import pyplot as plt

import statsmodels.api as sm
import statsmodels.formula.api as smf

### directories
base_dir = "/singleCell/Ruzicka_snRNA_Seq"
data_dir = base_dir + "/processed_standardized_data/"
data_dir2 = base_dir + "/data/"
out_dir = "/Aging/Aging_Clock/Output/"

#### the metadata of test and train.
train_df = pd.read_csv(data_dir2+ "22_11_02_df_train.csv")
val_df = pd.read_csv(data_dir2+ "22_11_02_df_val.csv")
train_val_df =pd.concat( [train_df, val_df], axis=0)

#### load corresponding data
## Data Documentation -- /singleCell/Ruzicka_snRNA_Seq/data/RAW/data_README.md
#* (17658 genes x 183056 cells, 69 controls)
sc=ad.read_h5ad(data_dir2+"22_08_12_seurat_raw_FILTERED_neurons_only_controls_only_lognormalized.h5ad")
print(sc.shape)
print(len(np.unique(sc.obs['ID'])))
sc = sc[sc.obs['ID'].isin(set(train_val_df['ID']))]
sc.var.index = sc.var['features']

print(sc.shape)
print(len(np.unique(sc.obs['ID'])))

### Add nCount_RNA / ID
sc.obs['nCount_RNA_ID'] = sc.obs['ID'].map(sc.obs.groupby('ID')['nCount_RNA'].mean())

## load iterative predicted age 
pred_age = pd.read_csv(out_dir+'24_10_02_model_M24_predictions_gene_filtering_mean_age_corrected_pred_age.csv', delimiter=",")
# left join model outputs with original metadata table
pred_age.index = pred_age['Id']
pred_age = pred_age[['predicted_age_0', 'predicted_age_24']]
sc.obs = sc.obs.merge(pred_age, left_index = True, right_index = True, how = "left")

#######################################
#######################################

###### perform feature normalisation
scaler = StandardScaler(with_mean=False)
scaler.fit(sc.X)  
cts = scaler.transform(sc.X)
#### scale nCount_RNA_ID
nscaler = StandardScaler()
nscaler.fit(sc.obs['nCount_RNA_ID'].values.reshape(-1, 1))
sc.obs['nCount_RNA_ID']=nscaler.transform(sc.obs['nCount_RNA_ID'].values.reshape(-1, 1))
#### scale the age
age_scaler = StandardScaler()
sc.obs["Age"]=age_scaler.fit_transform(sc.obs["Age"].values.reshape(-1, 1))
sc.obs["predicted_age_0"]=age_scaler.fit_transform(sc.obs["predicted_age_0"].values.reshape(-1, 1))
sc.obs["predicted_age_24"]=age_scaler.fit_transform(sc.obs["predicted_age_24"].values.reshape(-1, 1))

var = sc.var['features'] # features
df_coef = pd.DataFrame()
for i in range(var.shape[0]):
    print("feature"+str(i)+" "+var[i])
    # Create a DataFrame with the necessary variables
    data = pd.DataFrame({
    "response": cts[:, i].toarray().flatten(), # Convert sparse response matrix to a dense array
    "Age": sc.obs["Age"].values,  
    "nCount_RNA_ID": sc.obs["nCount_RNA_ID"].values, 
    "predicted_age_0": sc.obs["predicted_age_0"].values,
    "predicted_age_24": sc.obs["predicted_age_24"].values,
    "Celltype": sc.obs["Celltype"].values  # Grouping variable 
    })
    # Fit the mixed-effects model
    md = smf.ols("response ~ Age + nCount_RNA_ID", 
                 data)
    mdf = md.fit()
    
    coef = pd.DataFrame(mdf.params[1:3])
    coef.index = coef.index+"_coef"
    pval = pd.DataFrame(mdf.pvalues[1:3])
    pval.index = pval.index+"_pval"
    df = pd.concat([coef, pval], axis=0)
    df.columns= [var[i]]
    df_coef.index = df.index
    df_coef = df_coef.merge(df, left_index = True, right_index = True, how = "left")
    
datestring="25_02_09"
df_coef = df_coef.T
df_coef
df_coef.to_csv(out_dir+datestring+"_OLS_Ruzicka_age_fullset.csv")