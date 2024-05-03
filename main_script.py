
from scipy.stats import wilcoxon, mannwhitneyu, pearsonr, spearmanr, ttest_rel,ttest_ind
from scripts_dir.utils import *
from scripts_dir.conf_file import  *


if __name__=="__main__":
    conf = {
        'fea_color':'black',

    }

    df = prepare_file(data_file,rmv_samples,rmv_features,normal = conf['normal'])
    fea_dict = create_feature_dict(df,conf)
    group = ["1","2","12"][1]
    conf = {
        "cond1":"B","cond2": "B",
        # "cond1":"A","cond2": "A",
        "t1":0,"t2":90,
        "ptlist1": [2,5,7,11,13,14,17,18,20] , #[2,5,7],
        "ptlist2": [1,4,6,9,12,15,16], #[1,4,6,9],
        "G":"WBC",
        "B":"Plasma",
        "green":"WBC",
        "black":"Plasma",
        # "total":"total",
        'alpha':0.05,
        'test': wilcoxon,
        'temp':['FA(18:2)',"FA(16:0)","FA(18:1)",'C3','Tyrosine']
            # "Acetyl-CoA",
                # "ADPR",'Malic acid',
        # "beta-Hydroxybutyric acid","C6:1","PC aa C40:6","UDP-glucuronic acid"]
    }
    # cond1, cond2 = "A", "A"
    # t1,t2=0,90
    num_dim = 2
    group ="green"


