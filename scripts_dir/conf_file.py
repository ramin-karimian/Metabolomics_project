
data_file = "../data/dataset.xlsx" ## initial data file after doing manual modifications
# rmvlist = [3,5,6,8]
rmv_samples = [3,8] ## since they have missing data for whole condition B and A, respectively
rmv_features = ["FA(12:1)", ## all zero
                'Histamine', 'Phenylethylamine', 'cis-Hydroxyproline', 'Dopamine', 'DOPA',  ## they are not important
                'Carnosine' ,'Nitro-Tyrosine' ,'Diacetylspermine' , 'Tyramine', 'Phosphocreatine',   ## they are not important
                'FA(16:2)', 'FA(20:1)' ,'FA(20:2)', 'FA(20:3)' ,'FA(22:2)' ,'FA(22:3)', ## uncommon (in data2)
                'Methylgluaryl carnitine' ,'C14:1 carnitine' ,'C16:1 carnitine', 'C19 carnitine', 'C20:0 carnitine' ,'C20:3 carnitine',  ## uncommon (in data2)
                'FA(08:0-DC)' ,'FA(10:0-DC)' ,'FA(12:1-DC)' ,'FA(14:1)' ,'FA(16:0-OH)' ,'FA(18:0-OH)', 'Glutaryl carnitine', ## uncommon (in data1)
                # 'FA(16:1, 9Z)', ## uncommon (in data1)
                # 'FA(16:1)' ## uncommon (in data2)
]

conf = {

    "normal": True, # feature normalization
    "cond1":"B","cond2": "B",
    # "cond1":"A","cond2": "A",
    "t1":0,"t2":90,
    "ptlist1": [2,5,7,11,13,14,17,18,20] , #[2,5,7],
    "ptlist2": [1,4,6,9,12,15,16], #[1,4,6,9],
    "G":"WBC",
    "B":"Plasma",
    "green":"WBC",
    "black":"Plasma",
    'alpha':0.05,
    # 'test': wilcoxon,
    'temp':['FA(18:2)',"FA(16:0)","FA(18:1)",'C3','Tyrosine'],
        # "Acetyl-CoA",
            # "ADPR",'Malic acid',
    # "beta-Hydroxybutyric acid","C6:1","PC aa C40:6","UDP-glucuronic acid"]
    "font_size": 10,
    "dpi":500
}
