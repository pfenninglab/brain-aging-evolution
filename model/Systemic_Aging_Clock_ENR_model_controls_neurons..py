#!/usr/bin/env python3
 
# =============================================================================
# Single-round Elastic Net "aging clock" (Systemic Aging Model) for
# single-nucleus RNA-seq (neurons, controls)
# -----------------------------------------------------------------------------
# Goal: predict chronological age from per-cell gene expression using one
# cross-validated Elastic Net fit. Unlike the iterative "intrinsic" version,
# this script runs a SINGLE round (n_models = 1) with NO gene pruning:
#
#   Round 0 (the only round):
#       Fit a cross-validated Elastic Net to predict scaled Age from all genes
#       plus a per-donor sequencing-depth covariate (nCount_RNA_ID). Cells are
#       weighted by the reciprocal of their total counts. Record coefficients,
#       train/val errors (MAE, MSE) and R2, and save the predicted age.
#
# Results (coefficients, MSE, MAE, R2) and the full obs table with predictions
# are written to CSV at the end.
# =============================================================================

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
# Donor-level split tables (one row per sample/ID). Train + val are used here;
# the held-out test set is left commented out so it never leaks into fitting.
train_df = pd.read_csv(data_dir2+ "22_11_02_df_train.csv")
val_df = pd.read_csv(data_dir2+ "22_11_02_df_val.csv")
train_val_df =pd.concat( [train_df, val_df], axis=0)
#test_df = pd.read_csv(data_dir+ "22_11_02_df_test.csv")

#### load corresponding data
## Data Documentation -- /singleCell/Ruzicka_snRNA_Seq/data/RAW/data_README.md
# Log-normalized, neuron-only, control-only expression matrix (cells x genes).
#* (17658 genes x 183056 cells, 69 controls)
sc=ad.read_h5ad(data_dir2+"22_08_12_seurat_raw_FILTERED_neurons_only_controls_only_lognormalized.h5ad")
print(sc.shape)
print(len(np.unique(sc.obs['ID'])))

# Keep only cells belonging to donors in the train+val split, and use gene
# symbols as the var index so we can subset by gene name later.
sc = sc[sc.obs['ID'].isin(set(train_val_df['ID']))]
sc.var.index = sc.var['features']
print(sc.shape)
print(len(np.unique(sc.obs['ID'])))

### Add nCount_RNA / ID
sc.obs['nCount_RNA_rc'] = np.reciprocal(sc.obs['nCount_RNA'])
sc.obs['nCount_RNA_ID'] = sc.obs['ID'].map(sc.obs.groupby('ID')['nCount_RNA'].mean())
sc.obs['nCount_RNA_ID_CV'] = sc.obs['ID'].map(sc.obs.groupby('ID')['nCount_RNA'].std() / sc.obs.groupby('ID')['nCount_RNA'].mean()*100)

###################################
###################################


###################################
###################################

## Forloop in ROUNDS -- n iterations

df_coef = pd.DataFrame()
df_mae = pd.DataFrame()
df_r2 = pd.DataFrame()

# ENR - Systemic Aging Model
n_models= 1
for i in range(n_models):
    
    print("Model "+str(i))
    
    if i == 0 :
        
        ## ROUND 1 - Age
        ######################################################################
        ###### Preprocessing. 
        
        print(str(sc.var_names.shape[0])+" features used for fitting")
         
        ###### Dataset   
        ###### separate into train/val
        train = sc[sc.obs['ID'].isin(set(train_df['ID']))]
        val = sc[sc.obs['ID'].isin(set(val_df['ID']))]
        
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
        datestring="24_10_02"
        
        # make a model with distinct random state.
        m = ElasticNetCV(random_state=i, 
                         #normalize=False,
                         #cv = cv_iterable,
                         cv = 5,
                         selection="random")
      
        m.fit(X =hstack((cts,train.obs[["nCount_RNA_ID"]].values)).astype(np.float64),
              y =np.array(train.obs["Age"].values.tolist()),
              sample_weight = np.array(train.obs["nCount_RNA_rc"].values.tolist()))
            
        # LOAD THE MODEL...
            
        # write the coefficients.
        ftlist = sc.var['features'].tolist()
        ftlist.append('nCount_RNA_ID')
        #ftlist.append('Age')        
        df = pd.DataFrame()
        df.index = ftlist
        coef = deepcopy(m.coef_)
        #coef = np.append(coef,'NA')
        df["M"+"_"+str(i)] = coef
        df_coef = df
        #df.to_csv(out_dir+datestring+"_model_M"+str(i)+"_coef_optim_mean_age_corrected.csv",index=True)
            
        ######## record the model details.
        r2_res.loc[i,"l1_ratio"] = m.l1_ratio_
        r2_res.loc[i,"alpha"] = m.alpha_
            
        ######## record the results on train.
        train_res=m.predict(X =hstack((cts,train.obs[["nCount_RNA_ID"]].values)).astype(np.float64))

        mse_res.loc[i,"train"] = mean_squared_error(age_scaler.inverse_transform(train.obs["Age"].values.reshape(-1, 1)), 
                                                         age_scaler.inverse_transform(train_res.reshape(-1, 1)))
        mae_res.loc[i,"train"] = mean_absolute_error(age_scaler.inverse_transform(train.obs["Age"].values.reshape(-1, 1)), 
                                                         age_scaler.inverse_transform(train_res.reshape(-1, 1)))

        print(age_scaler.inverse_transform(train.obs["Age"].values.reshape(-1, 1)))
        print(age_scaler.inverse_transform(train_res.reshape(-1, 1)))      

            
        temp = pd.DataFrame(train.obs)
        temp["predicted_age"] =age_scaler.inverse_transform(train_res.reshape(-1, 1))
        #temp["Age"]= age_scaler.inverse_transform(train.obs["Age"].values.reshape(-1, 1))
        #temp.to_csv(out_dir+datestring+"_model_M"+str(i)+"_values_TRAIN_optim_mean_age_corrected.csv")
        
        r2_res.loc[i,"train"] = m.score(X =hstack((cts,train.obs[["nCount_RNA_ID"]].values)).astype(np.float64),
                                        y =np.array(train.obs["Age"].values.tolist()),
                                        sample_weight = np.array(train.obs["nCount_RNA_rc"].values.tolist()))
        print(train.obs["Age"])
            
        ######## record the results on val.
        print("VAL FIT STARTED")
        val_res=m.predict(X =hstack((cts_val,val.obs[["nCount_RNA_ID"]].values)).astype(np.float64))

        mse_res.loc[i,"val"] = mean_squared_error(age_scaler.inverse_transform(val.obs["Age"].values.reshape(-1, 1)),
                                                    age_scaler.inverse_transform(val_res.reshape(-1, 1)))
        mae_res.loc[i,"val"] = mean_absolute_error(age_scaler.inverse_transform(val.obs["Age"].values.reshape(-1, 1)),
                                                    age_scaler.inverse_transform(val_res.reshape(-1, 1)))
              
        print(age_scaler.inverse_transform(val.obs["Age"].values.reshape(-1, 1)))
        print(age_scaler.inverse_transform(val_res.reshape(-1, 1)))      
        
        temp2 = pd.DataFrame(val.obs)
        temp2["predicted_age"] =age_scaler.inverse_transform(val_res.reshape(-1, 1))
        #temp2["Age"]= age_scaler.inverse_transform(val.obs["Age"].values.reshape(-1, 1))
        #temp2.to_csv(out_dir+datestring+"_model_M"+str(i)+"_values_VAL_optim_mean_age_corrected.csv")
        
        r2_res.loc[i,"val"] = m.score(X = hstack((cts_val,val.obs[["nCount_RNA_ID"]].values)).astype(np.float64),
                                      y =np.array(val.obs["Age"].values.tolist()),
                                      sample_weight = np.array(val.obs["nCount_RNA_rc"].values.tolist()))
        df_mse = mse_res
        df_mae = mae_res
        df_r2 = r2_res
    
        print(val.obs["Age"])
        print("fit made on model"+str(i))         
        print("mae res")
        print(mae_res)
        print("r2_res")
        print(r2_res)

        # Merge predicted age from current fit to metadata 
        col_name = "predicted_age"+"_"+str(i)
        print(col_name)
        pred_age = pd.concat([temp[["predicted_age"]],temp2[["predicted_age"]]])
        pred_age = pred_age.rename(columns={"predicted_age":col_name})
        sc.obs = sc.obs.merge(pred_age, left_index = True, right_index = True, how = "left")
