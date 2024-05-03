import operator
import os
import pickle
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc, RocCurveDisplay, confusion_matrix, precision_recall_curve
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict, cross_val_score, LeaveOneOut, StratifiedKFold

import statsmodels.api
from scripts_dir.utils import *
from scripts_dir.conf_file import *
from scripts_dir.pca import separate_samples


def run_model(x,y):
    m = np.zeros(x.shape[1]).tolist()
    for i in range(100):
        print(i)
        rf = RandomForestClassifier(n_jobs=-1, n_estimators=200)
        model = rf.fit(x, y)
        m1 =  model.feature_importances_.tolist()
        m = list(map(operator.add, m,m1))
    # to make a table for important features
    n = x.columns.tolist()
    np.append([m],[n], axis=0)
    fe_imp= pd.DataFrame(np.append([m],[n], axis=0)).T
    fe_imp[0] = fe_imp[0].astype(float)#**6
    fe_imprt = fe_imp.sort_values(by=[0], ascending=False)#.iloc[0:20]

    return fe_imprt

def plot_model_accuracy_prcurve(x, y):
    # cv= 4
    cv= LeaveOneOut() # changed
    # y_tests = []
    # y_predicts = []
    # y_predictions = []
    # thresholds_track = [] # mine
    rf_confs = [] # mine
    pvals = []
    tprs = []
    aucs = []
    mean_fpr = np.linspace (0, 1, 100)
    plt.figure(figsize=(10, 7.5))
    for i in range(10):
        rf = RandomForestClassifier(n_jobs=-1, n_estimators=200, max_features = 'auto')
        rf_proba = cross_val_predict(rf, x,y, cv=cv, method='predict_proba')
        rf_confusion = cross_val_predict(rf, x, y, cv=cv)
        # rf_confusion = [-1*(x-1) for x in rf_confusion]
        # y = [-1*(x-1) for x in y]

        rf_conf = confusion_matrix(y,rf_confusion)
        print("Hi\n",rf_proba,"\n",y,"\n",rf_confusion,"\n",rf_conf,"\n\n")

        _, p_value = fisher_exact(rf_conf)
        pvals.append(p_value)
        rf_confs.append(rf_conf)

        rf_scores = rf_proba[:, 1] #Get the probability of the positive class
        fpr, tpr, thresholds = roc_curve(y,rf_scores)
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.5,c='steelblue',linestyle='--')
        # y_tests.append(y)
        # y_predicts.append(rf_scores)
        # y_predictions.append(rf_confusion)
    plt.plot([0, 1], [0, 1], linestyle='-.', lw=2, color='black',label='Luck', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)

    rf_confs_mean = np.mean(rf_confs,axis=0)
    _, p_value = fisher_exact(rf_confs_mean)

    plt.plot(mean_fpr, mean_tpr, color='teal', label=r'Mean ROC' ,lw=3, alpha=1)

    # print(p_value)
    # print(pvals)
    # print(rf_confs)

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(r'ROC Curve -- Classifier: Random Forests -- AUC = %0.2f (Pvalue = %0.4f)' % (mean_auc,p_value))
    plt.legend(loc="lower right")
    plt.savefig(f"../outputs/rf_model_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][0][1:]}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][1][1:]}.png",dpi = conf['dpi'],bbox_inches = 'tight')
    plt.show()

    f = open("test.txt",'w')
    [f.write(f"[[{x[0][0]},{x[0][1]}],[{x[1][0]},{x[1][1]}]]\n") for x in  rf_confs]
    f.close()
    return aucs, rf_proba, rf_confs, pvals, p_value
    # return scores, y_tests, y_predicts, y_predictions, aucs, rf_proba, thresholds_track, rf_confs, pvals, p_value



if __name__=="__main__":
    conf['group'] = ['ptlist1','ptlist2','all'][2]
    conf['normal'] = True
    conf['topn_imp_features'] = 20
    conf['keep_time_condition'] = [['A0','A90']][0]
    convertor={"A":1,"B":0,"0":0,"90":1}
    if conf['keep_time_condition'][0][0]!=conf['keep_time_condition'][1][0]:
        col2 = 'condition'
        conf['factors'] = ['individual','condition']
    else:
        col2 = 'time'
        conf['factors'] = ['individual','time']
    ## load and preprocess data
    df = prepare_file(data_file,rmv_samples,rmv_features,normal = conf['normal'], remove_metadata = True)

    ## create feature dictionary
    if "features_dictionary.pkl" not in os.listdir("../outputs") :
        fea_dict = create_feature_dict(df, conf)
    else:
        with open(f"../outputs/features_dictionary.pkl","rb" ) as f :
            fea_dict = pickle.load(f)

    ## keep only significant features
    df.drop(columns = [x for x in df.columns if fea_dict[x]['pvalue(A0_A90)'] >= 0.05 ],inplace=True)

    ## seprate samples of interest for analysis
    df , _ = separate_samples(df.copy(),conf)
    conf['num_of_patients'] = max(set([int(x.split("_")[0]) for x in df.index ]))


    dfx = df.copy()
    y = [convertor[x.split("_")[1 if col2=='time' else 2]] for x in df.index ]
    dfx['target']= y
    dfy = dfx.drop(columns=[x for x in dfx.columns if x!="target"])
    dfx.drop(columns=["target"],inplace = True)

    fe_imprt = run_model(dfx,y)
    aucs, rf_proba, rf_confs, pvals, p_value = plot_model_accuracy_prcurve(df, y)
    # run_model_accuracy(dfx.values,dfy.values)


    for i in range(len(fe_imprt)):
        importance , fea = fe_imprt.iloc[i,:]
        fea_dict[fea]['rf_importance'] = importance
    for fea in fea_dict:
        if fea not in fe_imprt.iloc[:,1].values:
            fea_dict[fea]['rf_importance'] = 0

    with open(f"../outputs/features_dictionary.pkl","wb" ) as f :
        pickle.dump(fea_dict,f)
        pd.DataFrame(fea_dict).T.to_csv(f"../outputs/features_dictionary.csv")
        fe_imprt_df = pd.DataFrame(fe_imprt)
        fe_imprt_df.columns=['rf_importance','fea']
        fe_imprt_df.to_csv("../outputs/features_imp.csv")

