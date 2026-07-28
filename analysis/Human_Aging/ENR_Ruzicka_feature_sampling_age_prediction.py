#!/usr/bin/env python3

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
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
import scipy.stats as stats
from scipy.stats import pearsonr
from scipy.sparse import coo_matrix, hstack
from matplotlib import pyplot as plt
from copy import deepcopy

### directories
base_dir = "/singleCell/Ruzicka_snRNA_Seq"
data_dir = base_dir + "/processed_standardized_data/"
data_dir2 = base_dir + "/data/"
out_dir = "/Aging/Aging_Clock/Output/"

#### the metadata of test and train.
train_df = pd.read_csv(data_dir2+ "22_11_02_df_train.csv")
val_df = pd.read_csv(data_dir2+ "22_11_02_df_val.csv")
train_val_df =pd.concat( [train_df, val_df], axis=0)
test_df = pd.read_csv(data_dir2+ "22_11_02_df_test.csv")

#### load corresponding data
## Data Documentation -- /singleCell/Ruzicka_snRNA_Seq/data/RAW/data_README.md
# 72 controls — 17658 features across 262411 samples
sc=ad.read_h5ad(data_dir2+"22_08_12_seurat_raw_FILTERED_neurons_only_controls_only_lognormalized.h5ad")
print(sc.shape)
print(len(np.unique(sc.obs['ID'])))
#(262411, 2000)
#72
sc = sc[sc.obs['ID'].isin(set(train_val_df['ID']))]
sc.var.index = sc.var['features']
print(sc.shape)
print(len(np.unique(sc.obs['ID'])))
#(224645, 2000)
#58
### Add nCount_RNA / ID
sc.obs['nCount_RNA_rc'] = np.reciprocal(sc.obs['nCount_RNA'])
sc.obs['nCount_RNA_ID'] = sc.obs['ID'].map(sc.obs.groupby('ID')['nCount_RNA'].mean())
sc.obs['nCount_RNA_ID_CV'] = sc.obs['ID'].map(sc.obs.groupby('ID')['nCount_RNA'].std() / sc.obs.groupby('ID')['nCount_RNA'].mean()*100)

###################################
###################################

# Informative feautres
geneset = pd.read_csv(out_dir+"24_10_02_model_M24_coef_gene_filtering_mean_age_corrected_pred_age.csv", index_col = 0)
pos = pd.DataFrame(geneset[geneset['M_24']> 0].index)
neg = pd.DataFrame(geneset[geneset['M_24']< 0].index)
features = pd.concat([pos, neg], ignore_index = True)[0].tolist()
features.remove('nCount_RNA_ID')
len(features)
#sc = sc[:,features]

# Random feautres
#random = np.random.choice(sc.var['features'],len(features), replace = False)
#sc1 = sc[:,random]
#sc1


###################################
###################################

## Forloop in ROUNDS -- n iterations

df_coef = pd.DataFrame()
df_mae = pd.DataFrame()
df_r2 = pd.DataFrame()

n_models= 100
for i in range(n_models):
    print("Model "+str(i))
        
    ## ROUND 1 - Age
    ######################################################################
    ###### Preprocessing. 

    if i == 0 : 
        sc1 = sc[:,features]
    if i > 0 : 
        # Random feautres
        random = np.random.choice(sc.var['features'],len(features), replace = False)
        sc1 = sc[:,random]
        print(str(sc1.var_names.shape[0])+" features used for fitting")
     
    ###### Dataset   
    ###### separate into train/val
    train = sc1[sc1.obs['ID'].isin(set(train_df['ID']))]
    val = sc1[sc1.obs['ID'].isin(set(val_df['ID']))]
    
    ###### perform feature normalisation over the training data, 
    #  Then perform normalisation on testing instances as well, 
    # but this time using the mean and variance of training explanatory variables.
    scaler = StandardScaler(with_mean=False)
    scaler.fit(train.X)  # it should be fit to the train data, and used on the test.
    cts = scaler.transform(train.X)
    cts_val = scaler.transform(val.X)
    print(cts.shape)
    print(cts_val.shape)
    
    #### scale nCount_RNA_ID
    nscaler = StandardScaler()
    nscaler.fit(train.obs['nCount_RNA_ID'].values.reshape(-1, 1))
    train.obs['nCount_RNA_ID']=nscaler.transform(train.obs['nCount_RNA_ID'].values.reshape(-1, 1))
    val.obs['nCount_RNA_ID']=nscaler.transform(val.obs['nCount_RNA_ID'].values.reshape(-1, 1))
    #### scale weights
    wscaler = MinMaxScaler()
    wscaler.fit(train.obs["nCount_RNA_rc"].values.reshape(-1, 1))
    train.obs["nCount_RNA_rc"]=wscaler.transform(train.obs["nCount_RNA_rc"].values.reshape(-1, 1))
    val.obs["nCount_RNA_rc"]=wscaler.transform(val.obs["nCount_RNA_rc"].values.reshape(-1, 1))
    #### scale the age
    age_scaler = StandardScaler()
    age_scaler.fit(train.obs["Age"].values.reshape(-1, 1))
    train.obs["Age"]=age_scaler.transform(train.obs["Age"].values.reshape(-1, 1))
    val.obs["Age"]=age_scaler.transform(val.obs["Age"].values.reshape(-1, 1))  
    
    r2_res = pd.DataFrame(columns=["alpha", "l1_ratio", "train", "val"], index=range(0))
    mse_res = pd.DataFrame(columns=["train", "val"], index=range(0))
    mae_res = pd.DataFrame(columns=["train", "val"], index=range(0))
    pearson_res = pd.DataFrame(columns=["train", "val"], index=range(0))
    
    ###### prediction task 
    # fit the elastic net.
    
    print("fit started on model"+str(i))
    
    gc.collect()
    datestring="24_05_01"
    
    # make a model with distinct random state.
    m = ElasticNetCV(cv = 5)
                     #random_state=i, 
                     #normalize=False,
                     #cv = cv_iterable,
                     #selection="random")
  
    m.fit(X =hstack((cts,train.obs[["nCount_RNA_ID"]].values)).astype(np.float64),
          y =np.array(train.obs["Age"].values.tolist()),
          sample_weight = np.array(train.obs["nCount_RNA_rc"].values.tolist())) 
        
    ######## record the model details.
    r2_res.loc[i,"l1_ratio"] = m.l1_ratio_
    r2_res.loc[i,"alpha"] = m.alpha_
        
    ######## record the results on train.
    train_res=m.predict(X =hstack((cts,train.obs[["nCount_RNA_ID"]].values)).astype(np.float64))

    mse_res.loc[i,"train"] = mean_squared_error(age_scaler.inverse_transform(train.obs["Age"].values.reshape(-1, 1)), 
                                                     age_scaler.inverse_transform(train_res.reshape(-1, 1)))
    mae_res.loc[i,"train"] = mean_absolute_error(age_scaler.inverse_transform(train.obs["Age"].values.reshape(-1, 1)), 
                                                     age_scaler.inverse_transform(train_res.reshape(-1, 1)))
     
       
    r2_res.loc[i,"train"] = m.score(X =hstack((cts,train.obs[["nCount_RNA_ID"]].values)).astype(np.float64),
                                    y =np.array(train.obs["Age"].values.tolist()),
                                    sample_weight = np.array(train.obs["nCount_RNA_rc"].values.tolist()))
        
    ######## record the results on val.
    print("VAL FIT STARTED")
    val_res=m.predict(X =hstack((cts_val,val.obs[["nCount_RNA_ID"]].values)).astype(np.float64))

    mse_res.loc[i,"val"] = mean_squared_error(age_scaler.inverse_transform(val.obs["Age"].values.reshape(-1, 1)),
                                                age_scaler.inverse_transform(val_res.reshape(-1, 1)))
    mae_res.loc[i,"val"] = mean_absolute_error(age_scaler.inverse_transform(val.obs["Age"].values.reshape(-1, 1)),
                                                age_scaler.inverse_transform(val_res.reshape(-1, 1)))
          
    if i == 0 : 
        df_mse = mse_res
        df_mae = mae_res
        df_r2 = r2_res
        print("fit made on model"+str(i))         
        print("mae res")
        print(mae_res)
        print("r2_res")
        print(r2_res)

    if i > 0 : 
        df_mse = pd.concat([df_mse, mse_res])
        df_mae = pd.concat([df_mae, mae_res])
        df_r2 = pd.concat([df_r2, r2_res])
        print("fit made on model"+str(i))         
        print("mae res")
        print(mae_res)
        print("r2_res")
        print(r2_res)

print('writing output')
df_mse.to_csv(out_dir+"model_M"+str(i)+"_feature_sampling_mse_pa.csv")
df_mae.to_csv(out_dir+"model_M"+str(i)+"_feature_sampling_mae_pa.csv")
df_r2.to_csv(out_dir+"model_M"+str(i)+"_feature_sampling_r2_pa.csv")