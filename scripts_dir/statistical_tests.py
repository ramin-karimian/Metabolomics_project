import os

import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu, pearsonr, spearmanr, ttest_rel,ttest_ind,probplot, shapiro
import statsmodels.api
from conf_file import *
from scripts_dir.utils import  *
import pickle


def apply_test(df,conf,cond1,cond2,t1,t2,var,test):

    df1 = df[df["Condition"]==cond1]
    df2 = df[df["Condition"]==cond2]

    df3 = df1[df1["Time"]==t1]
    df4 = df2[df2["Time"]==t2]

    df3.index = df3['ID'].values
    df4.index = df4['ID'].values

    # if cond1!=cond2:
    #     for i in df3.index:
    #        if df3['ID'][i] not in conf['ptlist1']:
    #             df3.drop(index=[i],inplace = True )
    #             df4.drop(index=i,inplace = True )

    if len(df3)>len(df4):df3 = df3.loc[df4.index]
    else:df4 = df4.loc[df3.index]
    if (df3.loc[:,var].values - df4.loc[:,var].values ==[0]*len(df4.loc[:,var].values)).all():
        print(var,"\nbe careful of normalization for this pvalue ==10 and rank 1000")
        return 1000,10

    rnk,pvalue = test(df3.loc[:,var].values,df4.loc[:,var].values)

    return rnk,pvalue

def hypothesis_testing_phase(writer,df,fea_dict,alpha = 0.05,consider_normality=False):
    pvals = []
    if alpha != 0.05: alpha = conf['alpha']
    ct = conf["keep_time_condition"]
    cond1,cond2,t1,t2 =  ct[0][0], ct[1][0],int(ct[0][1:]), int(ct[1][1:])
    print(cond1,cond2,t1,t2)
    for var in df.columns:
        if var not in fea_dict.keys(): continue
        _,norm = probplot(df.loc[:,var].values)
        _,pvl = shapiro(df.loc[:,var].values)

        if consider_normality==False:
            test = wilcoxon
            fea_dict[var]['normal'] = None
        elif pvl <conf['alpha']:
            test = wilcoxon
            fea_dict[var]['normal'] = "NO"
        else:
            test = ttest_rel
            fea_dict[var]['normal'] = "YES"

        rnk,pvalue = apply_test(df,conf,cond1,cond2,t1,t2,var,test)

        fea_dict[var]['metabolite'] = var
        fea_dict[var][f"raw_pvalue({cond1}{t1}_{cond2}{t2})"] = pvalue
        pvals.append((var,pvalue))

    FDRres = statsmodels.stats.multitest.multipletests([pvl for v,pvl in pvals], alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)


    for i in range(len(pvals)):
        fea_dict[pvals[i][0]][f"pvalue({cond1}{t1}_{cond2}{t2})"] = FDRres[1][i]
        if FDRres[1][i] < alpha:
            fea_dict[pvals[i][0]][f"significance({cond1}{t1}_{cond2}{t2})"] = 1
        else:
            fea_dict[pvals[i][0]][f"significance({cond1}{t1}_{cond2}{t2})"] = 0


    res = pd.DataFrame(fea_dict).transpose()
    res.to_excel(writer,sheet_name= f"{cond1}{t1}_{cond2}{t2}")
    return res , fea_dict,FDRres

def combine_imp_fea_from_bothgroups(conf):
    df1 = pd.read_excel(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_ptlist1{'_sgnf' if conf['remove_nonsign'] else ''}/top{conf['topn_imp_features']}.xlsx",index_col=0)
    df2 = pd.read_excel(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_ptlist2{'_sgnf' if conf['remove_nonsign'] else ''}/top{conf['topn_imp_features']}.xlsx",index_col=0)
    common_fea = list(set(list(df1.index)).intersection(set(list(df2.index))))
    return common_fea

def statistical_test_pca_features(writer,df,fea_dict,conf,alpha = 0.05):
    common_fea = combine_imp_fea_from_bothgroups(conf)
    pvals = []
    if alpha != 0.05: alpha = conf['alpha']
    ct = conf["keep_time_condition"]
    cond1,cond2,t1,t2 =  ct[0][0], ct[1][0],int(ct[0][1:]), int(ct[1][1:])
    print(cond1,cond2,t1,t2)
    for var in common_fea:
        test = wilcoxon
        fea_dict[var]['normal'] = None
        rnk,pvalue = apply_test(df,conf,cond1,cond2,t1,t2,var,test)

        fea_dict[var]['metabolite'] = var
        fea_dict[var][f"raw_pvalue({cond1}{t1}_{cond2}{t2}_combined)"] = pvalue
        pvals.append((var,pvalue))

    for var in df.columns:
        if var not in common_fea and var not in ['ID','Time','Condition']:
            fea_dict[var][f"raw_pvalue({cond1}{t1}_{cond2}{t2}_combined)"] = None
            fea_dict[var][f"pvalue({cond1}{t1}_{cond2}{t2}_combined)"] = None
            fea_dict[var][f"significance({cond1}{t1}_{cond2}{t2}_combined)"] = None

    FDRres = statsmodels.stats.multitest.multipletests([pvl for v,pvl in pvals], alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)


    for i in range(len(pvals)):
        fea_dict[pvals[i][0]][f"pvalue({cond1}{t1}_{cond2}{t2}_combined)"] = FDRres[1][i]
        if FDRres[1][i] < alpha:
            fea_dict[pvals[i][0]][f"significance({cond1}{t1}_{cond2}{t2}_combined)"] = 1
        else:
            fea_dict[pvals[i][0]][f"significance({cond1}{t1}_{cond2}{t2}_combined)"] = 0


    res = pd.DataFrame(fea_dict).transpose()
    res.to_excel(writer,sheet_name= f"{cond1}{t1}_{cond2}{t2}_combined")
    return res , fea_dict,common_fea


if __name__=="__main__":
    conf['normal'] = False
    # conf['alpha'] = 0.05
    conf['topn_imp_features'] = 100
    conf['fig_size']= (6, 4)
    conf['fea_color'] = 'black'
    conf['keep_time_condition'] = [['A0','A90']][0]
    conf['remove_nonsign'] = [True,False][1]
    if conf['keep_time_condition'][0][0]!=conf['keep_time_condition'][1][0]:
        conf['factors'] = ['individual','condition']
    else:
        conf['factors'] = ['individual','time']
    writer = pd.ExcelWriter(f"../outputs/statistical_test_metabolites.xlsx")
    ## load and preprocess data
    df = prepare_file(data_file,rmv_samples,rmv_features,normal = conf['normal'], remove_metadata = False)
    conf['num_of_patients'] = max(set([int(x.split("_")[0]) for x in df.index ]))

    ## create feature dictionary
    if "features_dictionary.pkl" not in os.listdir("../outputs") :
        fea_dict = create_feature_dict(df, conf)
    else:
        with open(f"../outputs/features_dictionary.pkl","rb" ) as f :
            fea_dict = pickle.load(f)

    ## this is for all features together
    # res , fea_dict,FDRres =  hypothesis_testing_phase(writer,df,fea_dict,conf,consider_normality=False)

    # this is for pca extacted features for both groups combined features
    res , fea_dict,common_fea = statistical_test_pca_features(writer,df,fea_dict,conf,alpha = 0.05)

    # pvals = statistical_test_pca_features(writer,df,fea_dict,conf,alpha = 0.05)


    with open(f"../outputs/features_dictionary.pkl","wb") as f :
        pickle.dump(fea_dict,f)
    pd.DataFrame(fea_dict).T.to_csv(f"../outputs/features_dictionary.csv")
