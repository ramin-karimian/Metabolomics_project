import matplotlib.pyplot as plt
import numpy as np
import os

import pandas as pd
from skbio.diversity.alpha import shannon
from scipy.stats import wilcoxon
import statsmodels.api
from conf_file import *
from utils import *

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)

    if (a_set & b_set):
        s = a_set & b_set
    else:
        print(a_set , b_set)
        print("No common elements")
    return s

def statistical_test(df1,df2):
    df1.index = [x.split("_")[0] for x in df1.index ]
    df2.index = [x.split("_")[0] for x in df2.index ]

    ind = common_member(list(df1.index),list(df2.index))
    df1.drop(index = [x for x in df1.index if x not in ind], inplace = True)
    df2.drop(index = [x for x in df2.index if x not in ind], inplace = True)

    rnk, pvalue = wilcoxon(df1.values,df2.values)

    return pvalue


def diversity(df,start_index = 3):
    divs = {}
    df['alpha'] = None
    for ind in range(len(df.index)):
        div = shannon(df.iloc[ind,start_index:-1].values.tolist())
        divs[df.index[ind]] = div
        df.iloc[ind,-1] = div

    return df

def alpha_diversity_plot(patients,divs_df,conf,timeslots,ylim,diversity=False):
    df = divs_df[divs_df['Time'] <= timeslots]
    samp = np.unique(df['ID'])

    if diversity: dive = f"_{diversity}"
    else: dive = ""

    if patients == "ptlist1": xticks = ['A-0', 'A-90', 'A-180', 'B-0', 'B-90', 'B-180']
    else: xticks = ['B-0', 'B-90', 'B-180','A-0', 'A-90', 'A-180']

    if conf['timeslots']==180:
        xvalues= np.array([0,1,2,4,5,6])/2
    elif conf['timeslots']==90:
        xvalues= np.array([0,1,3,4])/2
        xticks = [x  for x in xticks if '180' not in x ]

    d = {}
    for p,v in zip(xticks,xvalues):
        d[p]=v
    plt.clf()
    plt.xticks(xvalues,xticks)

    for s in samp:
        if s not in conf[patients]:
            # print(s)
            # print(conf[patients])
            continue
        ser = divs_df[divs_df['ID']==s]["alpha"]
        xs = [k for k in d.keys() if k in [f"{divs_df.loc[x, 'Condition']}-{divs_df.loc[x, 'Time']}" for x in ser.index] ]
        # print([d[x] for x in xs],ser[[str(s)+"-"+x for x in xs]].values)
        ser.index = [f"{str(s)}-{divs_df.loc[x, 'Condition']}-{divs_df.loc[x, 'Time']}" for x in ser.index]
        plt.plot([d[x] for x in xs],ser[[str(s)+"-"+x for x in xs]].values,marker="^",label = s)

    # plt.legend(loc="upper right")
    plt.legend(loc='center left',bbox_to_anchor=(1,0.5),fontsize = conf['font_size']-2)

    plt.title(f"Shannon diversity for condition { 'A => B'  if patients == 'ptlist1' else  'B => A'}", fontsize=conf['font_size'])
    plt.xlabel("Condition-Time",labelpad = 0.1, fontsize=conf['font_size'])
    plt.ylabel("Shannon diversity", fontsize=conf['font_size'] )
    # plt.ylim(ylim[0],ylim[1])
    txt=f"Alpha diversities for patients {conf[patients]} (A --> B)" if patients=="ptlist1" else f"Alpha diversities for patients {conf[patients]} (B --> A)"
    # plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=conf['font_size'])
    plt.show()
    if timeslots <=90:
        plt.savefig(f'../outputs/diversities/90/{patients}{dive}.png', dpi = conf['dpi'] , bbox_inches = 'tight')
    else:
        plt.savefig(f'../outputs/diversities/180/{patients}{dive}.png', dpi = conf['dpi'] , bbox_inches = 'tight')


def alpha_diversity_boxplot(patients,divs_df,conf,timeslots,diversity=False):


    cond1 , cond2, t1, t2= 'A','B',0,90
    divs_df = divs_df[divs_df['Time'] <= timeslots]

    divs_df.drop(index = [s for s in divs_df.index if int(s.split("_")[0]) not in conf[patients]],inplace = True)

    # divs_df = divs_df['alpha']
    divs_df.dropna(inplace =True)
    # print(len(divs_df))
    divs_df1 = divs_df[divs_df['Condition']==cond1]
    divs_df2 = divs_df1[divs_df1['Time']==t1]

    divs_df1 = divs_df[divs_df['Condition']==cond1]
    divs_df3 = divs_df1[divs_df1['Time']==t2]

    divs_df1 = divs_df[divs_df['Condition']==cond2]
    divs_df4 = divs_df1[divs_df1['Time']==t1]

    divs_df1 = divs_df[divs_df['Condition']==cond2]
    divs_df5 = divs_df1[divs_df1['Time']==t2]


    if diversity: dive = f"_{diversity}"
    else: dive = ""


    if patients == "ptlist1": xticks = ['A-0', 'A-90', 'B-0', 'B-90']
    else: xticks = ['B-0', 'B-90','A-0', 'A-90']


    xvalues= np.array([1,2,3,4])
    xticks = [x  for x in xticks if '180' not in x ]

    d = {}
    for p,v in zip(xticks,xvalues):
        d[p]=v
    plt.clf()
    plt.xticks(xvalues,xticks)
    x1,x2,x3,x4 = divs_df2['alpha'].values , divs_df3['alpha'].values, divs_df4['alpha'].values, divs_df5['alpha'].values

    if patients == "ptlist1":
        c1 , c2 = 'black', 'white'
        plt.boxplot(
            [x1,x2,x3,x4],
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
        print(f"XX1: \n {divs_df2['alpha'],divs_df3['alpha'],divs_df4['alpha'],divs_df5['alpha']}")
        print(f"XX11: \n {divs_df2['alpha'].mean(),divs_df3['alpha'].mean(),divs_df4['alpha'].mean(),divs_df5['alpha'].mean()}")
        pval1 = statistical_test(divs_df2['alpha'],divs_df3['alpha'])
        pval2 = statistical_test(divs_df3['alpha'],divs_df4['alpha'])
        pval3 = statistical_test(divs_df4['alpha'],divs_df5['alpha'])
        pvalues = [pval1,pval2,pval3]
        FDRres = statsmodels.stats.multitest.multipletests(pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        corrected_pvalues = FDRres[1]
        pvalues_df = pd.DataFrame([pvalues,corrected_pvalues],index= ['pvalues','corrected_pvalues'], columns = ['1st-2nd','2nd-3rd','3rd-4th'])

    else:
        c1,c2 = 'black', 'white'
        plt.boxplot(
            [x3,x4,x1,x2],
                    labels=[f"{cond2}{t1}",
                            f"{cond2}{t2}",
                            f"{cond1}{t1}",
                            f"{cond1}{t2}",
                            ],
                    meanline=True, showmeans=True, patch_artist=True,
                    boxprops=dict(facecolor=c2, color=c1),
                    capprops=dict(color=c1),
                    whiskerprops=dict(color=c1),
                    flierprops=dict(color=c1, markeredgecolor=c2),
                    medianprops=dict(color=c1)
                    )
        print(f"XX2: \n {divs_df4['alpha'],divs_df5['alpha'],divs_df2['alpha'],divs_df3['alpha']}")
        print(f"XX22: \n {divs_df4['alpha'].mean(),divs_df5['alpha'].mean(),divs_df2['alpha'].mean(),divs_df3['alpha'].mean()}")

        pval1 = statistical_test(divs_df4['alpha'],divs_df5['alpha'])
        pval2 = statistical_test(divs_df5['alpha'],divs_df2['alpha'])
        pval3 = statistical_test(divs_df2['alpha'],divs_df3['alpha'])
        pvalues = [pval1,pval2,pval3]
        FDRres = statsmodels.stats.multitest.multipletests(pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        corrected_pvalues = FDRres[1]
        pvalues_df = pd.DataFrame([pvalues,corrected_pvalues],index= ['pvalues','corrected_pvalues'], columns = ['1st-2nd','2nd-3rd','3rd-4th'])



    plt.title(f"Shannon diversity for condition { 'A => B'  if patients == 'ptlist1' else  'B => A'}", fontsize=conf['font_size'])
    plt.xlabel("Condition-Time",labelpad = 0.1, fontsize=conf['font_size'])
    plt.ylabel("Shannon diversity", fontsize=conf['font_size'] )
    txt=f"Alpha diversities for patients {conf[patients]} (A --> B)" if patients=="ptlist1" else f"Alpha diversities for patients {conf[patients]} (B --> A)"
    plt.show()
    if timeslots <=90:
        plt.savefig(f'../outputs/diversities/90/{patients}{dive}.png', dpi = conf['dpi'] , bbox_inches = 'tight')
    else:
        plt.savefig(f'../outputs/diversities/180/{patients}{dive}.png', dpi = conf['dpi'] , bbox_inches = 'tight')


    return  pvalues, corrected_pvalues , pvalues_df

if __name__ == "__main__":

    ## override initial conf file
    conf['normal'] = False ## for shannon we do not normalize the data
    conf['timeslots'] = 90 ## if 90, ignores 180 timepoint for all samples
    conf['font_size'] = 10
    ## load and preprocess data
    df = prepare_file(data_file,rmv_samples,rmv_features,normal = conf['normal'], remove_metadata = False)

    ## calculate diversities
    divs_df = diversity(df.copy(), start_index = 3)

    ## plot diversities

    # alpha_diversity_plot('ptlist1',divs_df,conf,conf['timeslots'],ylim=(3,4.25))
    # alpha_diversity_plot('ptlist2',divs_df,conf,conf['timeslots'],ylim=(3,4.25))

    ## plot diversity boxplots
    pvalues1, corrected_pvalues1 , pvalues_df1 = alpha_diversity_boxplot('ptlist1',divs_df.copy(),conf,conf['timeslots'])
    pvalues2, corrected_pvalues2 , pvalues_df2 = alpha_diversity_boxplot('ptlist2',divs_df.copy(),conf,conf['timeslots'])

    pvalues_df1.T.to_csv(f"../outputs/diversities/90/pvalues1.csv")
    pvalues_df2.T.to_csv(f"../outputs/diversities/90/pvalues2.csv")
