import pandas as pd
import numpy as np
from distinctipy import distinctipy

def normalize_dict(d):
    values = list(d.values())
    keys = list(d.keys())
    values = normalize_list(values)
    d1 = {}
    for i in range(len(keys)):
        d1[keys[i]]=values[i]
    return d1

def normalize_list(l):
    l2 = []
    max = np.max(l)
    min = np.min(l)
    for i in range(len(l)):
        l2.append((l[i]-min)/(max-min))
    return l2

def normalize(df):
    for c in df.columns:
        if c not in ["ID","Time","Condition"]:
            # print(c)
            df[c] = df[c].apply(lambda x: (x-min(df[c]))/(max(df[c])-min(df[c])) )
    return df

def prepare_file(file,rmvlist,colsrmv,normal =True,remove_metadata=False):
    # df = pd.read_excel(file)
    df1 = pd.read_excel(file,sheet_name='data_1')
    df2 = pd.read_excel(file,sheet_name='data_2')
    df = pd.concat([df1,df2])
    df.index = range(len(df))

    with pd.ExcelWriter(file, engine='openpyxl', mode='a',  if_sheet_exists='replace') as writer:
        df.to_excel(writer, sheet_name='data_concat')

    rmv_index = []
    for i in rmvlist:
        rmv_index.extend(list(df[df['ID']==i].index.values))
    df.drop(index = rmv_index, inplace =True)

    df.drop(columns=colsrmv,inplace = True)
    df.replace("< 0",0,inplace = True)
    df.replace("< LOD",0,inplace = True)

    df.index = range(len(df))

    with pd.ExcelWriter(file, engine='openpyxl', mode='a',  if_sheet_exists='replace') as writer:
        df.to_excel(writer, sheet_name='data_concat_preprocessed')

    if normal :
        df = normalize(df)
        with pd.ExcelWriter(file, engine='openpyxl', mode='a',  if_sheet_exists='replace') as writer:
            df.to_excel(writer, sheet_name='data_concat_preprocessed_normal')

    ## fix index
    inds = []
    for i in range(len(df)):
        inds.append(f"{df.iloc[i,0]}_{df.iloc[i,1]}_{df.iloc[i,2]}")
    df.index = inds

    with pd.ExcelWriter(file, engine='openpyxl', mode='a',  if_sheet_exists='replace') as writer:
        df.to_excel(writer, sheet_name='data_final')

    if remove_metadata:
        df.drop(columns=["ID","Condition","Time"],inplace = True)


    return df


## this is changed in the new project directory since all the features are now Plasma,and WBC features are excluded from the analysis

def create_feature_dict(df, conf):
    d = {}
    for x in df.columns:
        if x in ['ID','Time', 'Condition']: continue
        if x in d.keys():print(x, "Repeated feature")
        d[x] = {}
        # d[x]['color'] = conf['fea_color']
    return d



def create_distinct_colors(n_clrs):
    colors_list = distinctipy.get_colors(n_clrs) ## returns 3 digit tupe for each color
    return colors_list
