import networkx as nx
import pandas as pd
import numpy as np
import os , pickle
import community as community_louvain
import statsmodels.api
from scipy.spatial.distance import pdist
from scipy import stats

from scripts_dir.utils import *
from conf_file import  *

def prepare_df(df,t1,t2,cond):
    xf = df[df['Condition']==cond]
    x1 = xf[xf['Time']==t1]
    x2 = xf[xf['Time']==t2]
    del xf ,df
    x1.index = x1.ID
    x2.index = x2.ID

    x1,x2  = x1.iloc[:,3:] , x2.iloc[:,3:]

    xf = x1 - x2
    xf.index = x2.index
    return xf

def Louvain_modularity_community(g):
    partition = community_louvain.best_partition(g)
    for n in g.nodes():
        g.nodes[n]["louvain_modularity"] = partition[n]
    return  g


def centrality_phase(g):
    res = normalize_dict(nx.betweenness_centrality(g))
    res1= normalize_dict(nx.closeness_centrality(g))
    res2 = normalize_dict(nx.degree_centrality(g))
    # res = nx.betweenness_centrality(g)
    # res1= nx.closeness_centrality(g)
    # res2 = nx.degree_centrality(g)
    for n in g.nodes:
        g.nodes[n]['betweenness_centrality'] = res[n]
        g.nodes[n]['closeness_centrality'] = res1[n]
        g.nodes[n]['degree_centrality'] = res2[n]
    return g

def calculate_distance(df):
    # Calculate the pairwise distances between samples using the euclidean distance metric
    distance_matrix = pdist(df.values, metric='euclidean')
    return distance_matrix

def correlation(df,dropNan = False):
    # xf = prepare_df(df_alpha,t1,t2,cond,group,fea_dict,dropNan)
    # xf = xf.transpose()
    fea_names= df.index.values
    if dropNan:
        df=df.dropna()
    corr,pvalue =stats.spearmanr(df.values,axis = 1) # If axis=0 (default), then each column represents a variable,
    # with observations in the rows. If axis=1, the relationship is transposed: each row represents
    # a variable, while the columns contain observations. If axis=None, then both arrays will be
    # raveled.
    corrected_pvalues = np.ones((len(corr),len(corr)))
    FDRres= statsmodels.stats.multitest.multipletests([pvalue[i,j] for i in range(len(pvalue)) for j in range(i+1,len(pvalue)) ], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    # FDR_pvalues = FDRres[1].reshape((len(pvalue),-1))
    c = 0
    for i in range(len(corrected_pvalues)):
        for j in range(i+1,len(corrected_pvalues)):
            corrected_pvalues[i][j] = FDRres[1][c]
            c+=1

    return corr,pvalue, corrected_pvalues, fea_names,FDRres[1]

def construct_network(dist,pvalue,fea_names,fea_dict,th_pvalue,th_corr):
    g_pos = nx.Graph()
    g_neg = nx.Graph()
    g = nx.Graph()

    for i in range(len(dist)):
        for j in range(i+1,len(pvalue)):
            if pvalue[i][j] <= th_pvalue and dist[i][j] >= th_corr:
                if fea_names[i] not in g.nodes:
                        # exec("%s = %d" % (f"{f}",5000))
                    # g.add_node(fea_names[i],label=fea_dict[fea_names[i]]['color'])
                    g.add_node(fea_names[i])
                    g_pos.add_node(fea_names[i])
                    g_neg.add_node(fea_names[i])

                    for k,v in fea_dict[fea_names[i]].items():
                        if k in ['normal','metabolite','pvalue(A0_A180)','raw_pvalue(A0_A180)','significance(A0_A180)','pvalue(A0_B0)','raw_pvalue(A0_B0)','significance(A0_B0)','pvalue(B0_B90)','raw_pvalue(B0_B90)','significance(B0_B90)']: continue
                        g.nodes[fea_names[i]][k]=v
                        g_pos.nodes[fea_names[i]][k]=v
                        g_neg.nodes[fea_names[i]][k]=v

                if fea_names[j] not in g.nodes:
                    # g.add_node(fea_names[j],label=fea_dict[fea_names[j]]['color'])
                    g.add_node(fea_names[j])
                    g_pos.add_node(fea_names[j])
                    g_neg.add_node(fea_names[j])
                    for k,v in fea_dict[fea_names[j]].items():
                        if k in ['normal','metabolite','pvalue(A0_A180)','raw_pvalue(A0_A180)','significance(A0_A180)','pvalue(A0_B0)','raw_pvalue(A0_B0)','significance(A0_B0)','pvalue(B0_B90)','raw_pvalue(B0_B90)','significance(B0_B90)']: continue
                        g.nodes[fea_names[j]][k]=v
                        g_pos.nodes[fea_names[j]][k]=v
                        g_neg.nodes[fea_names[j]][k]=v

                # g.add_edge(fea_names[i],fea_names[j],weight= 1 - pvalue[i][j])
                g.add_edge(fea_names[i],fea_names[j],weight= np.around(abs(dist[i][j]),decimals=2) , pvalue = pvalue[i][j],sign = np.sign(dist[i][j]))
                if dist[i][j]>=0:
                    g_pos.add_edge(fea_names[i],fea_names[j],weight= np.around(abs(dist[i][j]),decimals=2) , pvalue = pvalue[i][j],sign = np.sign(dist[i][j]) )
                    # g_pos.add_edge(g.nodes[fea_names[i]],g.nodes[fea_names[j]])
                else:
                    g_neg.add_edge(fea_names[i],fea_names[j],weight= np.around(abs(dist[i][j]),decimals=2), pvalue = pvalue[i][j], sign = np.sign(dist[i][j]))
                    # g_neg.add_edge(g.nodes[fea_names[i]],g.nodes[fea_names[j]])
    # g = centrality_phase(g)
    # g_neg = centrality_phase(g_neg)
    # g_pos = centrality_phase(g_pos)
    # g =Louvain_modularity_community(g)
    # g_neg =Louvain_modularity_community(g_neg)
    # g_pos =Louvain_modularity_community(g_pos)

    return g,g_neg,g_pos

if __name__=="__main__":

    conf['normal'] = True
    conf['topn_imp_features'] = 10
    conf['num_dim'] = 2
    conf['marker_size'] = 50
    conf['x_ulimit'], conf['y_ulimit'] , conf['x_llimit'] ,conf['y_llimit']  = 4.1, 3, -4.1, -2.5
    conf['x-y_llimit'] = -4.1
    conf['fig_size']= (6, 4)
    conf['fea_color'] = 'black'
    conf['keep_time_condition'] = [['A0','A90'],['B0','B90'],['A0','B0'],['A90','B90']][0]
    conf['group'] = ['ptlist1','ptlist2','all'][2]
    conf['remove_nonsign'] = [True,False][0]
    conf['approach'] = ['stat_rf','pca_stat'][1]
    conf['rscript_path'] = "permanova_script.R"
    conf['th_pvalue'] = 0.05
    conf['th_corr'] = 0.6
    conf['topn'] = 40  # from topn features from stat_rf or pca_stat approach
    # conf['factors'] = [['individual','time'],['individual','condition']][0]
    if conf['keep_time_condition'][0][0]!=conf['keep_time_condition'][1][0]:
        conf['factors'] = ['individual','condition']
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

    if conf['approach'] =='stat_rf':
        ## removed non-significant features
        if conf['remove_nonsign']:
            df.drop(columns=[x for x in df.columns[3:] if fea_dict[x][f'significance(A0_A90)']!=1],inplace=True)

        ## just keep the rf_model importatn features
        index_from_rf = pd.DataFrame(fea_dict).T["rf_importance"].sort_values()[::-1].iloc[:conf['topn']].index
        df.drop(columns=[x for x in df.columns[3:] if x not in index_from_rf],inplace = True)
    elif conf['approach'] =='pca_stat':
        ## just keep the combined (patient groups) pca imp features
        df.drop(columns=[x for x in df.columns[3:] if fea_dict[x][f'significance(A0_A90_combined)']!=1],inplace=True)



    ## get the A90-A0 df
    df_dif = prepare_df(df,t1=int(conf['keep_time_condition'][0][1:]),t2=int(conf['keep_time_condition'][1][1:]),cond=conf['keep_time_condition'][0][0])
    df_dif = df_dif.T

    ## preprocess data again
    print(f"df shape Before preprocessing (dropna) : {df.shape}")
    df.dropna(inplace =True)
    print(f"df shape After preprocessing (dropna) : {df.shape}")
    for cl in df.columns:
        if cl not in fea_dict.keys():
            df.drop(columns=[cl],inplace= True)


    # pd.DataFrame(df.columns).to_csv(f"../outputs/network/{conf['topn']}_{str(conf['th_pvalue']).split('.')[1]}_{str(conf['th_corr']).split('.')[1]}/{conf['approach']}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][0][1:]}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][1][1:]}_important_metabolites.csv")

    corr,pvalue, corrected_pvalues, fea_names,res = correlation(df_dif)
    # dist = calculate_distance(df_dif)
    g,g_neg,g_pos = construct_network(corr,corrected_pvalues,fea_names,fea_dict,conf['th_pvalue'],conf['th_corr'])

    if f"{conf['topn']}_{str(conf['th_pvalue']).split('.')[1]}_{str(conf['th_corr']).split('.')[1]}" not in os.listdir("../outputs/network"):
        os.mkdir(f"../outputs/network/{conf['topn']}_{str(conf['th_pvalue']).split('.')[1]}_{str(conf['th_corr']).split('.')[1]}")
    nx.write_gexf(g,f"../outputs/network/{conf['topn']}_{str(conf['th_pvalue']).split('.')[1]}_{str(conf['th_corr']).split('.')[1]}/{conf['approach']}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][0][1:]}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][1][1:]}_g.gexf")
    pd.DataFrame(corr,index = df_dif.index, columns = df_dif.index ).to_excel(f"../outputs/network/{conf['topn']}_{str(conf['th_pvalue']).split('.')[1]}_{str(conf['th_corr']).split('.')[1]}/{conf['approach']}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][0][1:]}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][1][1:]}_corrs.xlsx")
    pd.DataFrame(corrected_pvalues,index = df_dif.index, columns = df_dif.index ).to_excel(f"../outputs/network/{conf['topn']}_{str(conf['th_pvalue']).split('.')[1]}_{str(conf['th_corr']).split('.')[1]}/{conf['approach']}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][0][1:]}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][1][1:]}_corrected-pvalues.xlsx")
    pd.DataFrame(pvalue,index = df_dif.index, columns = df_dif.index ).to_excel(f"../outputs/network/{conf['topn']}_{str(conf['th_pvalue']).split('.')[1]}_{str(conf['th_corr']).split('.')[1]}/{conf['approach']}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][0][1:]}_{conf['keep_time_condition'][0][0]}{conf['keep_time_condition'][1][1:]}_pvalues.xlsx")
