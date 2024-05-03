import os
import pickle

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist
from scipy.spatial import distance_matrix
from scripts_dir.utils import normalize_list
import statsmodels.api
from scripts_dir.utils import *
from scripts_dir.conf_file import *
# from skbio.stats.distance  import permanova as skbio_permanova
from skbio.stats.distance  import DistanceMatrix
from scipy.spatial.distance import squareform
from scripts_dir.boxplots import boxplot
import permanova
import os

def execute_cmd(cmd):
    print("Running:\n",cmd)
    return_code = os.system(cmd)
    print(f"Done! (return_code:{return_code})\n")
    return return_code

def permanova_test(df,conf, method='euclidean'):
    col1 , col2 = conf['factors']
    c1,t1,c2,t2 = conf['keep_time_condition'][0][0],conf['keep_time_condition'][0][1:],conf['keep_time_condition'][1][0],conf['keep_time_condition'][1][1:]
    if col2 == 'time':
        df[col1] = [x.split("_")[0] for x in df.index]
        df[col2] = [x.split("_")[1] for x in df.index]
    elif col2 == 'condition':
        df[col1] = [x.split("_")[0] for x in df.index]
        df[col2] = [x.split("_")[2] for x in df.index]
    elif col2 == 'both':
        col1, col2 = col1, 'time'
        df[col1] = [x.split("_")[0] for x in df.index]
        df[col2] = [x.split("_")[2] for x in df.index]
    df.to_csv(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}/temp/df.csv")
    df.drop(columns=[col1,col2],inplace=True)
    cmd = f"Rscript --vanilla {conf['rscript_path']} {c1} {t1} {c2} {t2} {col1} {col2} {conf['group']}"
    execute_cmd(cmd)
    print(f"Attention: ",f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}/Rpermanova_results.csv")
    pvalue_df = pd.read_csv(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}/Rpermanova_results.csv",index_col=0)
    return pvalue_df


def extract_important_features(data_reduced,features,fea_dict,group='total',topn= 10):
    dis_from_center = pd.Series(index = range(len(data_reduced)))
    for i in range(len(data_reduced)):
        dis_from_center[i] = np.sqrt(np.sum(data_reduced[i]**2))
    sorted_dis_from_center = dis_from_center.sort_values()[-1:0:-1]
    sorted_indexes = sorted_dis_from_center.index[:]
    sorted_dis_from_center.index = features[sorted_indexes]
    # ss = pd.Series(index= features[dis_from_center.index[:]])
    # ss[:] = normalize_list(dis_from_center.values)

    for f in sorted_dis_from_center.index.values:
        fea_dict[f][f'pca_value_{group}'] = sorted_dis_from_center[f]
    ## Note : dis_from_center is soretd based on features distance from center ( high to low) with index being the feature name
    return sorted_dis_from_center ,sorted_indexes, fea_dict


def pca_plot_data(df,model,colors_list,conf):
    # plt.clf()
    df_reduced = model.transform(df)

    x_ulimit, y_ulimit = conf['x_ulimit'], conf['y_ulimit']
    x_llimit, y_llimit = conf['x_llimit'], conf['y_llimit']

    labels=[]
    s= conf['marker_size']
    # plt.figure(figsize=conf['fig_size'], dpi=conf['dpi'])


    for i in range(len(df_reduced)):
        if int(df.index[i].split("_")[0]) in conf['ptlist1']:m="^"
        else: m = "o"
        clr =  colors_list[int(df.index[i].split("_")[0])-1]
        if int(df.index[i].split("_")[0]) not in labels:
            labels.append(int(df.index[i].split("_")[0]))
            plt.scatter(df_reduced[i,0],df_reduced[i,1],c = clr, s = s,
                        marker = m, label = int(df.index[i].split("_")[0]) )
        else :
            plt.scatter(df_reduced[i,0],df_reduced[i,1],c = clr , s=s,marker = m)

        plt.text(df_reduced[i,0],df_reduced[i,1],"_".join(df.index[i].split("_")[1:]),fontsize = conf['font_size']-5,c = 'black')
    txt = f"PCA plot for {conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][0][1:]}_{conf['keep_time_condition'][1][0]}{conf['keep_time_condition'][1][1:]}"
    plt.title(txt)
    plt.xlabel(f"PC1 ({np.around(model.explained_variance_ratio_[0],decimals = 2)})",fontsize = conf['font_size'])
    plt.ylabel(f"PC2 ({np.around(model.explained_variance_ratio_[1],decimals = 2)})",fontsize = conf['font_size'])
    plt.legend(loc='center left',bbox_to_anchor=(1,0.5),fontsize = conf['font_size']-2)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    x1,x2 = ax.get_xlim()
    step_2 = (x2-x1)/5
    plt.xticks([ x1 + i*step_2 for i in range(6)], [ np.around(x1 + i*step_2,decimals=2) for i in range(6)],fontsize = conf['font_size']-2)
    ax.spines['bottom'].set_bounds(x1,x2)

    plt.savefig(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}/{txt}.png",dpi=conf['dpi'],bbox_inches = 'tight')
    # plt.close()


def plot_pca_features_old(model,features,fea_dict,conf):

    data_reduced = model.components_.T
    plt.scatter(data_reduced[:,0],data_reduced[:,1],s = 2,c= 'black')

    x_ulimit,y_ulimit= 0.2 , 0.2
    x_llimit , y_llimit =- 0.1 , -0.2

    # plt.xlim(x_llimit,x_ulimit)
    # plt.ylim(y_llimit,y_ulimit)
    # # plt.axis('scaled')

    dis_from_center ,sorted_indexes, fea_dict = extract_important_features(data_reduced,features,fea_dict,topn = conf['topn_imp_features']) # ss (sorted)==> fea. : value, sorted_s ==> index:value

    # dist_mat = dist_matrix(data_reduced)

    for i in range(len(data_reduced)):
        if features[i] in dis_from_center[:conf['topn_imp_features']]: font,c = conf['font_size']-7, 'red'
        else : font,c = conf['font_size']-7, 'black'
        plt.text(data_reduced[i,0],data_reduced[i,1],features[i],fontsize = font,c = c)

    plt.plot([0,0],[y_ulimit,y_llimit],linestyle='dashed',c = 'gray')
    plt.plot([x_ulimit,x_llimit],[0,0],linestyle='dashed',c = 'gray')

    plt.xlabel(f"PC1({np.around(model.explained_variance_ratio_[0],decimals = 2)})", labelpad = 0.1,fontsize = conf['font_size'])
    plt.ylabel(f"PC2({np.around(model.explained_variance_ratio_[1],decimals = 2)})", labelpad = 0.1,fontsize = conf['font_size'])
    # plt.title(f"Metabolites",fontsize = conf['font_size'])

    plt.savefig(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}/pca_features.png",dpi = conf['dpi'],bbox_inches = 'tight')
    plt.close()
    return dis_from_center, fea_dict


def plot_pca_features(model,features,fea_dict,conf):

    sorted_dis_from_center ,sorted_indexes, fea_dict = extract_important_features(model.components_.T,features,fea_dict,group=f"{conf['group']}",topn = conf['topn_imp_features']) # ss (sorted)==> fea. : value, sorted_s ==> index:value

    ratio = conf['ratio']
    for i in range(len(model.components_.T)):
        if features[i] in sorted_dis_from_center[:conf['topn_imp_features_ongraph']]:
            font,c = conf['font_size']-7, 'darkcyan'
            plt.text(model.components_.T[i,0]*ratio,model.components_.T[i,1]*ratio,features[i],fontsize = font + 2,c = c)
            plt.arrow(0, 0, model.components_.T[i,0]*ratio,model.components_.T[i,1]*ratio ,width=0.00005,color='gray',alpha=0.6,overhang=0.5,head_width=0.05, head_length=0.01)
            print(features[i])
        # else : font,c = conf['font_size']-7, 'black'
        else: continue


    # plt.xlabel(f"PC1({np.around(model.explained_variance_ratio_[0],decimals = 2)})", labelpad = 0.1,fontsize = conf['font_size'])
    # plt.ylabel(f"PC2({np.around(model.explained_variance_ratio_[1],decimals = 2)})", labelpad = 0.1,fontsize = conf['font_size'])

    plt.savefig(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}/pca_features.png",dpi = conf['dpi'],bbox_inches = 'tight')
    # plt.close()
    return sorted_dis_from_center, fea_dict


def pca_model(df,fea_dict, conf):

    model = PCA(n_components=conf['num_dim'])

    X = model.fit_transform(df)
    features = df.columns

    dis_from_center , fea_dict = plot_pca_features(model,features,fea_dict,conf)
    # clrs = ["red","orange","blue","green","yellow","brown","purple","pink","gray"]
    colors_list = create_distinct_colors(conf['num_of_patients'])
    pca_plot_data(df,model,colors_list, conf)
    return fea_dict


def separate_samples(df,conf):
    inds = []
    if conf['group']!="all":
        df.drop(index = [y for y in df.index if int(y.split("_")[0]) not in conf[conf['group']]], inplace=True)
    for x in conf['keep_time_condition']:
        df_temp = df.drop(index = [y for y in df.index if y.split("_")[-1] != x[0]])
        df_temp = df_temp.drop(index = [y for y in df_temp.index if y.split("_")[1] !=x[1:]])
        inds.extend(list(df_temp.index))

    df.drop(index = [x for x in df.index if x not in inds ],inplace = True)
    return df, inds

if __name__ == "__main__":

    conf['normal'] = True
    conf['ratio']= 3  # 10 for all non-sign. ones
    conf['topn_imp_features'] = 100
    conf['topn_imp_features_ongraph'] = 10 ## just to show on PCA graphs
    conf['num_dim'] = 2
    conf['marker_size'] = 50
    conf['x_ulimit'], conf['y_ulimit'] , conf['x_llimit'] ,conf['y_llimit']  = 4.1, 3, -4.1, -2.5
    conf['x-y_llimit'] = -4.1
    conf['fig_size']= (6, 4)
    conf['fea_color'] = 'black'
    conf['keep_time_condition'] = [['A0','A90'],['B0','B90'],['A90','B0'],['B90','A0'],['A90','B90'],['A0','B0']][5]
    conf['group'] = ['ptlist1','ptlist2','all'][0]
    conf['remove_nonsign'] = [True,False][1]
    conf['rscript_path'] = "permanova_script.R"
    # conf['factors'] = [['individual','time'],['individual','condition']][0]
    if conf['keep_time_condition'][0][0]!=conf['keep_time_condition'][1][0]:
        if conf['keep_time_condition'][0][1:]==conf['keep_time_condition'][1][1:]:
            conf['factors'] = ['individual','condition']
        else:
            conf['factors'] = ['individual','time']

    else:
        conf['factors'] = ['individual','time']

    if f"pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}" not in os.listdir("../outputs"):
        os.mkdir(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}")
    if "temp" not in os.listdir(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}"):
        os.mkdir(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}/temp")

    ## load and preprocess data
    df = prepare_file(data_file,rmv_samples,rmv_features,normal = conf['normal'], remove_metadata = False)
    conf['num_of_patients'] = max(set([int(x.split("_")[0]) for x in df.index ]))

    ## create feature dictionary
    if "features_dictionary.pkl" not in os.listdir("../outputs") :
        fea_dict = create_feature_dict(df, conf)
    else:
        with open(f"../outputs/features_dictionary.pkl","rb" ) as f :
            fea_dict = pickle.load(f)

    ## removed non-significant features
    if conf['remove_nonsign'] :
        df.drop(columns=[x for x in df.columns[3:] if fea_dict[x][f'significance(A0_A90)']!=1],inplace=True)

    ## preprocess data again
    print(f"df shape Before preprocessing (dropna) : {df.shape}")
    df.dropna(inplace =True)
    print(f"df shape After preprocessing (dropna) : {df.shape}")
    for cl in df.columns:
        if cl not in fea_dict.keys():
            print(f"{pca_model.__name__}_warning {cl} not in keys")
            df.drop(columns=[cl],inplace= True)

    ## seprate samples of interest for analysis
    df , _ = separate_samples(df.copy(),conf)


    ## run PCA model
    fea_dict = pca_model(df,fea_dict, conf)

    ## statistical tests
    pvalue_df = permanova_test(df,conf, method = 'euclidean')
    fea_df = pd.DataFrame(fea_dict).T[f'pca_value_{conf["group"]}']
    fea_df = fea_df.dropna()
    top_fea_df = fea_df.sort_values()[::-1][:conf['topn_imp_features']]
    top_fea_df.to_excel(f"../outputs/pca_{'-'.join(conf['keep_time_condition'])}_{conf['group']}{'_sgnf' if conf['remove_nonsign'] else ''}/top{conf['topn_imp_features']}.xlsx")


    # dff = prepare_file(data_file,rmv_samples,rmv_features,normal = conf['normal'], remove_metadata = False)
    # boxplot(dff,top_fea_df)

    if conf['keep_time_condition'] == ['A0','A90'] :
        with open(f"../outputs/features_dictionary.pkl","wb" ) as f :
            pickle.dump(fea_dict,f)
        pd.DataFrame(fea_dict).T.to_csv(f"../outputs/features_dictionary.csv")

