from scripts_dir.utils import *
import matplotlib.pyplot as plt
import numpy as np
from scripts_dir.utils import *
from conf_file import  *
import pickle, os



def boxplot(dff):
    cond1 , cond2, t1, t2= 'A','B',0,90
    dff.dropna(inplace =True)
    print(len(dff))
    dff1 = dff[dff['Condition']==cond1]
    dff2 = dff1[dff1['Time']==t1]
    del dff1
    dff1 = dff[dff['Condition']==cond1]
    dff3 = dff1[dff1['Time']==t2]
    del dff1
    dff1 = dff[dff['Condition']==cond1]
    dff3_2 = dff1[dff1['Time']==180]
    del dff1
    dff1 = dff[dff['Condition']==cond2]
    dff4 = dff1[dff1['Time']==t1]
    del dff1
    dff1 = dff[dff['Condition']==cond2]
    dff5 = dff1[dff1['Time']==t2]
    del dff1
    dff1 = dff[dff['Condition']==cond2]
    dff5_2 = dff1[dff1['Time']==180]
    del dff1

    # del dff
    for fea in dff.columns:
        if fea in ['ID','Time','Condition']: continue
        # plt.subplot(100+i)
        # if group =="black":
        c1,c2 = 'black', 'white'
        # elif group == "green": c1,c2 = 'lime', 'white',
        plt.clf()
        x1= dff2[fea]
        x2= dff3[fea]
        x2_2= dff3_2[fea]

        x3= dff4[fea]
        x4= dff5[fea]
        x4_2= dff5_2[fea]

        # c1,c2 = 'blue', 'cyan'
        plt.boxplot(
            [x1,x2,x3,x4],
            # [x1,x2,x2_2,x3,x4,x4_2],

                    labels=[f"{cond1}{t1}",
                            f"{cond1}{t2}",
                            # f"{cond1}{180}",
                            f"{cond2}{t1}",
                            f"{cond2}{t2}",
                            # f"{cond2}{180}"
                            ],
                    meanline=True, showmeans=True, patch_artist=True,
                    boxprops=dict(facecolor=c2, color=c1),
                    capprops=dict(color=c1),
                    whiskerprops=dict(color=c1),
                    flierprops=dict(color=c1, markeredgecolor=c2),
                    medianprops=dict(color=c1)
                    )
        txt = f"{fea}"
        plt.title(txt, pad = 20,fontdict={"fontsize":12})
        print(fea)

        if "/" in fea : fea = fea.replace("/","_")
        if f"imp_metabolites_from_{conf['approach']}" not in os.listdir("../outputs"):
            os.mkdir(f"../outputs/imp_metabolites_from_{conf['approach']}")
        plt.savefig(f"../outputs/imp_metabolites_from_{conf['approach']}/{fea.replace(':','')}.png",bbox_inches = 'tight',dpi = 500)



def boxplot_old(dff,ss,group='black'):
    # cond1,cond2,t1,t2,ptlist1,ptlist2 = conf["cond1"],conf["cond2"],conf["t1"],conf["t2"], conf["ptlist1"], conf["ptlist2"]
    cond1 , cond2, t1, t2= 'A','B',0,90
    dff.dropna(inplace =True)
    # dff.index = dff['ID']
    # dff=dff.loc[conf['ptlist1']]
    print(len(dff))
    dff1 = dff[dff['Condition']==cond1]
    dff2 = dff1[dff1['Time']==t1]
    del dff1
    dff1 = dff[dff['Condition']==cond1]
    dff3 = dff1[dff1['Time']==t2]
    del dff1
    dff1 = dff[dff['Condition']==cond1]
    dff3_2 = dff1[dff1['Time']==180]
    del dff1
    dff1 = dff[dff['Condition']==cond2]
    dff4 = dff1[dff1['Time']==t1]
    del dff1
    dff1 = dff[dff['Condition']==cond2]
    dff5 = dff1[dff1['Time']==t2]
    del dff1
    dff1 = dff[dff['Condition']==cond2]
    dff5_2 = dff1[dff1['Time']==180]
    del dff1

    del dff
    for i in range(10):
        # plt.subplot(100+i)
        if group =="black": c1,c2 = 'black', 'white'
        elif group == "green": c1,c2 = 'lime', 'white',
        plt.clf()
        x1= dff2[ss.index[i]]
        x2= dff3[ss.index[i]]
        x2_2= dff3_2[ss.index[i]]

        x3= dff4[ss.index[i]]
        x4= dff5[ss.index[i]]
        x4_2= dff5_2[ss.index[i]]

        # c1,c2 = 'blue', 'cyan'
        plt.boxplot(
            [x1,x2,x3,x4],
            # [x1,x2,x2_2,x3,x4,x4_2],

                    labels=[f"{cond1}{t1}",
                            f"{cond1}{t2}",
                            # f"{cond1}{180}",
                            f"{cond2}{t1}",
                            f"{cond2}{t2}",
                            # f"{cond2}{180}"
                            ],
                    meanline=True, showmeans=True, patch_artist=True,
                    boxprops=dict(facecolor=c2, color=c1),
                    capprops=dict(color=c1),
                    whiskerprops=dict(color=c1),
                    flierprops=dict(color=c1, markeredgecolor=c2),
                    medianprops=dict(color=c1)
                    )
        txt = f"{ss.index[i]} _ Rank #{i+1}"
        plt.title(txt, pad = 20,fontdict={"fontsize":12})
        print(ss.index[i])
        s=ss.index[i]
        if "/" in ss.index[i] : s = ss.index[i].replace("/","_")
        plt.savefig(f"../outputs/pca/boxplots/{group}_Rank {i+1}_{s.replace(':','')}.png",bbox_inches = 'tight',dpi = 500)

def boxplot2(file,ss,rmvlist,colsrmv,conf):

    # cond1,cond2,t1,t2,ptlist1,ptlist2 = conf["cond1"],conf["cond2"],conf["t1"],conf["t2"], conf["ptlist1"], conf["ptlist2"]
    cond1,cond2,t1,t2 = 'A','B',0,90
    dff = prepare_file(file,rmvlist,colsrmv,normal =True,remove_metadata=False)
    dff.dropna(inplace =True)
    # dff.index = dff['ID']
    # dff=dff.loc[conf['ptlist1']]
    print(len(dff))
    dff1 = dff[dff['Condition']=="A"]
    dff2 = dff1[dff1['Time']==0]
    del dff1
    dff1 = dff[dff['Condition']=="A"]
    dff3 = dff1[dff1['Time']==90]
    del dff1
    dff1 = dff[dff['Condition']=="A"]
    dff3_2 = dff1[dff1['Time']==180]
    del dff1
    dff1 = dff[dff['Condition']=="B"]
    dff4 = dff1[dff1['Time']==0]
    del dff1
    dff1 = dff[dff['Condition']=="B"]
    dff5 = dff1[dff1['Time']==90]
    del dff1
    dff1 = dff[dff['Condition']=="B"]
    dff5_2 = dff1[dff1['Time']==180]
    del dff1
    del dff
    # for i in range(10):
        # plt.clf()
        # x1= dff2[ss.index[i]]
        # x2= dff3[ss.index[i]]
        # x3= dff4[ss.index[i]]
        # x4= dff5[ss.index[i]]
    for x in conf['temp']:
        plt.clf()
        x1= dff2[x]
        x2= dff3[x]
        x2_2= dff3_2[x]

        x3= dff4[x]
        x4= dff5[x]
        x4_2= dff5_2[x]

        c1,c2 = 'black', 'white'
        plt.boxplot([x1,x2,x2_2,x3,x4,x4_2],
                    labels=[f"{cond1}{t1}",
                            f"{cond1}{t2}",
                            f"{cond1}{180}",

                            f"{cond2}{t1}",
                            f"{cond2}{t2}",
                            f"{cond2}{180}",

                            ],
                    meanline=True, showmeans=True, patch_artist=True,
                    boxprops=dict(facecolor=c2, color=c1),
                    capprops=dict(color=c1),
                    whiskerprops=dict(color=c1),
                    flierprops=dict(color=c1, markeredgecolor=c2),
                    medianprops=dict(color=c1)
                    )
        txt = f"{x}"
        plt.title(txt, pad = 20,fontdict={"fontsize":12})
        plt.savefig(f"figures/boxplots/new/{x.replace(':','')}.png")


def boxplot3(file,cond1,cond2,t1,t2,fea_dict,conf,rmvlist,colsrmv,varlist,xf):
    dff = prepare_file(file,rmvlist,colsrmv,normal =True,remove_metadata=False)
    dff.dropna(inplace =True)
    dff1 = dff[dff['Condition']=="A"]
    dff2 = dff1[dff1['Time']==0]
    del dff1
    dff1 = dff[dff['Condition']=="A"]
    dff3 = dff1[dff1['Time']==90]
    del dff1
    dff1 = dff[dff['Condition']=="A"]
    dff4 = dff1[dff1['Time']==180]
    del dff1

    del dff

    for x in varlist:
        plt.clf()
        x1= dff2[x]
        x2= dff3[x]
        # x3= dff4[x]
        # x4= dff5[x]
        c1,c2 = 'black', 'white'
        plt.boxplot(
            # [x1,x2,x3,x4],
                    [x1,x2],
                    labels=['A0', 'A90'
                        # ,'A180'
                            # f"{x}\n{cond1}{t2}",
                            # f"{x}\n{cond2}{t1}",
                            # f"{cond2}{t2}"
                            ],
                    meanline=True, showmeans=True, patch_artist=True,
                    boxprops=dict(facecolor=c2, color=c1),
                    capprops=dict(color=c1),
                    whiskerprops=dict(color=c1),
                    flierprops=dict(color=c1, markeredgecolor=c2),
                    medianprops=dict(color=c1)
                    )
        txt = f"{x} ({conf[fea_dict[x]['color']]})\nFDR_pvalue = {np.around(xf[f'pvalue({cond1}{t1}_{cond2}{t2})'][x],decimals=3)}"
        plt.title(txt, pad = 20,fontdict={"fontsize":12})
        plt.savefig(f"figures/boxplots/significant/{cond1}{t1}_{cond2}{t2}/{x.replace(':','')}.png")




if __name__=="__main__":

    conf['normal'] = True
    conf['topn_imp_features'] = 10
    conf['num_dim'] = 2
    conf['fig_size']= (6, 4)
    conf['fea_color'] = 'black'
    conf['keep_time_condition'] = [['A0','A90'],['B0','B90'],['A0','B0'],['A90','B90']][0]
    conf['approach'] = ['stat_rf','pca_stat'][1]
    conf['group'] = ['ptlist1','ptlist2','all'][2]
    conf['remove_nonsign'] = [True,False][0]
    conf['rscript_path'] = "permanova_script.R"
    conf['topn'] = 20  # from topn features from rf model
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
        ## just keep the rf_model importatn features
        index_from_rf = pd.DataFrame(fea_dict).T["rf_importance"].sort_values()[::-1].iloc[:conf['topn']].index
        df.drop(columns=[x for x in df.columns[3:] if x not in index_from_rf],inplace = True)
    elif conf['approach'] =='pca_stat':
        ## just keep the combined (patient groups) pca imp features
        df.drop(columns=[x for x in df.columns[3:] if fea_dict[x][f'significance(A0_A90_combined)']!=1],inplace=True)



    # ## removed non-significant features
    # if conf['remove_nonsign'] :
    #     df.drop(columns=[x for x in df.columns[3:] if fea_dict[x][f'significance(A0_A90)']!=1],inplace=True)
    # index_from_rf = pd.DataFrame(fea_dict).T["rf_importance"].sort_values()[::-1].iloc[:conf['topn']].index
    #
    # ## just keep the rf_model importatn features
    # df.drop(columns=[x for x in df.columns[3:] if x not in index_from_rf],inplace = True)



    boxplot(df.copy())

