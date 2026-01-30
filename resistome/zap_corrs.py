import os
import pandas as pd
import numpy as np
import scipy.stats as stat
from multiprocessing import Pool
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore", category=RuntimeWarning) #this will otherwise clutter the screen when calculating corrs

os.chdir("../results/Corr_res")
#long_zap = pd.read_csv("long_zap.csv", index_col = 0)

# this is the task I want to parallelize
def permuted_spearman(item):
    # I need each processor to have access to the functions to use, otherwise it will result in a mutex
    import pandas as pd
    import scipy.stats as stat

    def statistic(x, y):  # explore all possible pairings by permuting `x`
        rs = stat.spearmanr(x, y).statistic  # ignore pvalue
        dof = len(x) - 2 # will only work for cases where x and y are equal
        transformed = rs * np.sqrt(dof / ((rs+1.0)*(1.0-rs)))
        return transformed
    group = item[1]
    corr = stat.permutation_test((group["var1_vals"], group["var2_vals"]), statistic, n_resamples = 999, alternative="two-sided", permutation_type="pairings")
    return [item[0][0], item[0][1], corr.statistic, corr.pvalue]

### ZAP RESULTS ####
#corr_res = []
#with Pool() as pool:
#    items = [(name, group) for name, group in long_zap.groupby(["var1", "var2"])]
#    print("There's ", len(items), " groups to calculate")
#    for result in tqdm(pool.imap(permuted_spearman, items)):
#        corr_res.append(result)
#corr_res = pd.DataFrame(corr_res, columns = ["var1", "var2", "statistic", "p-val"])
#corr_res = corr_res.loc[corr_res["p-val"] <= 0.05]
#print("permutations done, recalculating coef")
#for name, group in tqdm(corr_res.groupby(["var1", "var2"])):
#    num_data = long_zap.loc[(long_zap["var1"] == name[0]) & (long_zap["var2"] == name[1])]
#    coef = stat.spearmanr(num_data["var1_vals"], num_data["var2_vals"]).statistic
#    corr_res.loc[(corr_res["var1"] == name[0]) & (corr_res["var2"] == name[1]), "coef"] = coef
#corr_res.to_csv("ZAP_corrs.csv")

#del corr_res
#del long_zap

long_zac = pd.read_csv("long_zac.csv", index_col = 0)
corr_res = []
with Pool() as pool:
    print("making the item list")
    items = [(name, group) for name, group in long_zac.groupby(["var1", "var2"])]
    print("There's ", len(items), " groups to calculate")
    for result in tqdm(pool.imap(permuted_spearman, items)):
        corr_res.append(result)
corr_res = pd.DataFrame(corr_res, columns = ["var1", "var2", "statistic", "p-val"])
corr_res = corr_res.loc[corr_res["p-val"] <= 0.05]
print("permutations done, recalculating correlation coef")
for name, group in tqdm(corr_res.groupby(["var1", "var2"])):
    num_data = long_zac.loc[(long_zac["var1"] == name[0]) & (long_zac["var2"] == name[1])]
    coef = stat.spearmanr(num_data["var1_vals"], num_data["var2_vals"]).statistic
    corr_res.loc[(corr_res["var1"] == name[0]) & (corr_res["var2"] == name[1]), "coef"] = coef
corr_res.to_csv("ZAC_corrs.csv")
