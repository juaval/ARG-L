import os
import pandas as pd
import numpy as np
import warnings
import scipy.stats as stat

from multiprocessing import Pool
from tqdm import tqdm
warnings.filterwarnings("ignore", category=RuntimeWarning) #this will otherwise clutter the screen when calculating corrs

os.chdir("../results/Corr_res")
long_zap = pd.read_csv("long_zap.csv", index_col = 0)

def statistic(x, y):
    rs = stat.spearmanr(x, y).statistic  # ignore pvalue
    dof = len(x) - 2 # will only work for cases where x and y are equal
    transformed = rs * np.sqrt(dof / ((rs+1.0)*(1.0-rs)))
    return transformed

def permuted_spearman2(item):
    group = item[1]
    corr = stat.permutation_test((group["var1_vals"], group["var2_vals"]), statistic, n_resamples = 999, alternative="two-sided", permutation_type="pairings")
    return [item[0][0], item[0][1], corr.statistic, corr.pvalue]

corr_res = []
items = [(name, group) for name, group in long_zap.groupby(["var1", "var2"])]

print("There's ", len(items), " groups to calculate")
items_tenth = int(len(items)/10)
#items = items[items_tenth:items_tenth*2]
#items = items[:int(len(items)/2)]
#items = items[:int(len(items)/2)]
#items = items[:int(len(items)/2)]
#items_n = items[:int(len(items)/2)]

items_n = items[:items_tenth]
for item in tqdm(items_n):
    corr_res.append(permuted_spearman2(item))
print("1/10 done")
corr_res = pd.DataFrame(corr_res, columns = ["var1", "var2", "statistic", "p-val"])
corr_res = corr_res.loc[corr_res["p-val"] <= 0.05]
for name, group in tqdm(corr_res.groupby(["var1", "var2"])):
    num_data = long_zap.loc[(long_zap["var1"] == name[0]) & (long_zap["var2"] == name[1])]
    coef = stat.spearmanr(num_data["var1_vals"], num_data["var2_vals"]).statistic
    corr_res.loc[(corr_res["var1"] == name[0]) & (corr_res["var2"] == name[1]), "coef"] = coef
    corr_res.loc[corr_res["var1"]== name[0], "var1"] = name[0]
    corr_res.loc[corr_res["var2"]== name[1], "var2"] = name[1]


corr_res.to_csv("ZAP_corrs_1.csv")

corr_res = []
items_n = items[items_tenth*2:items_tenth*3]
for item in tqdm(items_n):
    corr_res.append(permuted_spearman2(item))
print("10/10 done")
corr_res = pd.DataFrame(corr_res, columns = ["var1", "var2", "statistic", "p-val"])
corr_res = corr_res.loc[corr_res["p-val"] <= 0.05]
for name, group in corr_res.groupby(["var1", "var2"]):
    num_data = long_zap.loc[(long_zap["var1"] == name[0]) & (long_zap["var2"] == name[1])]
    coef = stat.spearmanr(num_data["var1_vals"], num_data["var2_vals"]).statistic
    corr_res.loc[(corr_res["var1"] == name[0]) & (corr_res["var2"] == name[1]), "coef"] = coef
corr_res.to_csv("ZAP_corrs_3.csv")

del corr_res

print("DONE!")

