import pandas as pd
from tqdm import tqdm

def seq_str_cleaned(x):
    symbols = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'I', 'N', 'D', 'Q', 'K', 'E', 'M',
     'H', 'F', 'U', 'R', 'Y', 'W', 'S[79.96]', 'T[79.96]',
     'Y[79.96]', 'M[15.99]', 'C[57.02]', 'K[42.01]', 'K[229.16]']
    e_str = str(x) if str(x)!='nan' else ""
    e_str_cleaned = "".join([i for i in e_str if i in symbols+['']])
    return e_str_cleaned 

def seqlist(pre, seq, post, n_m, n_results):
    aa_mass = {
        "G": 57.021463735,
        "A": 71.037113805,
        "S": 87.032028435,
        "P": 97.052763875,
        "V": 99.068413945,
        "T": 101.047678505,
        "C": 103.009184505,
        "L": 113.084064015,
        "I": 113.084064015,
        "N": 114.042927470,
        "D": 115.026943065,
        "Q": 128.058577540,
        "K": 128.094963050,
        "E": 129.042593135,
        "M": 131.040484645,
        "H": 137.058911875,
        "F": 147.068413945,
        "U": 150.953633405,
        "R": 156.101111050,
        "Y": 163.063328575,
        "W": 186.079312980,
        "S[79.96]": 79.966331 + 87.032028435,
        "T[79.96]": 79.966331 + 101.047678505,
        "Y[79.96]": 79.966331 + 163.063328575,
        "M[15.99]": 15.994915 + 131.040484645,
        "C[57.02]": 57.021464 + 103.009184505,
        "K[42.01]": 42.010565 + 128.094963050,
        "K[229.16]": 229.16 + 128.094963050,
        "HOH": 1.00794 + 1.00794 + 15.9994
    }
    mod = {
        "S": "S[79.96]",
        "T": "T[79.96]",
        "Y": "Y[79.96]",
        "M": "M[15.99]",
        "C": "C[57.02]",
        "K": "K[42.01]"
        }
    
    pre_list, seq_list, post_list = (list(seq_str_cleaned(pre))[::-1], 
                                     list(seq_str_cleaned(seq)),
                                     list(seq_str_cleaned(post)))
    pre_vars = [[""]] + [pre_list[0:i+1][::-1] for i in range(len(pre_list))]
    post_vars = [[""]] + [post_list[0:i+1] for i in range(len(post_list))]
    pre_seq_vars = [pre_var + seq_list for pre_var in pre_vars]
    full_seq_vars = [pre_seq_var + post_var for pre_seq_var in pre_seq_vars for post_var in post_vars]
    full_seq_vars_mods = []
    for var in full_seq_vars:
        while '' in var:
            try:
                var.remove('')
            except:
                pass
        mods = [[]]
        for e in var:
            if e in mod.keys():
                mods_yes = [x + [mod[e]] for x in mods]
                mods_no = [x + [e] for x in mods]
                mods = mods_yes+mods_no
            else:
                modss = [x + [e] for x in mods]
                mods = modss
        for e in mods:
            if e[-1] == "K[42.01]":
                e[-1] = "K[229.16]"
        full_seq_vars_mods = full_seq_vars_mods+mods
    full_seq_vars_mods = [[seq_mod, sum([aa_mass[x] for x in seq_mod]) + aa_mass['HOH']] 
                                for seq_mod in full_seq_vars_mods]
    [x.append(x[1]-n_m) for x in full_seq_vars_mods]
    
    full_seq_vars_mods = [x for x in full_seq_vars_mods if abs(x[2]) <= (aa_mass['G'] + aa_mass['W'])]
    out = []
    all_aa  = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'I', 'N', 'D', 'Q', 'K', 'E', 'M',
     'H', 'F', 'U', 'R', 'Y', 'W', 'S[79.96]', 'T[79.96]',
     'Y[79.96]', 'M[15.99]', 'C[57.02]', 'K[42.01]', 'K[229.16]']
    change_table = [[x, y, aa_mass[y] - aa_mass[x]] for x in all_aa for y in all_aa]
    out = []
    for variant in full_seq_vars_mods:
        var = variant[0]
        var_change = sorted([x for x in change_table if x[0] in var], key = lambda f: abs(variant[2] + f[2]))
        result = [variant]
        for change in var_change:
            for i in range(len(var)):
                if var[i] == change[0]:
                    if not ((change[1] == "K[229.16]" and i != len(var)-1) or (change[1] == "K[42.01]" and i == len(var)-1)):
                        result +=[[var[0:i] + [change[1]]+var[i+1::], variant[1]+change[2], variant[2] + change[2]]]
            if len(result)>=n_results:
                break
            out+=result
    out =  [[''.join(x[0])] + x[1::] for x in out]         
    out = [list(x) for x in list(dict.fromkeys(tuple(i) for i in out))]      
    out = sorted(out, key = lambda f: abs(f[2]))[0:n_results]
    
    return out
    
df = pd.read_csv('data.csv') # insert your data
n_results = 30
columns_names = df.columns.to_list()
pre_seq = df['Pre_Seq'].tolist()
seq = df['Sequence'].tolist()
post_seq = df['Post_Seq'].tolist()
n_mass = df["spectrum neutral mass"].tolist()
result = []
for i in tqdm(range(len(df))):
    row = df.iloc[i].to_list()
    row_updated_set = [row + x for x in seqlist(pre_seq[i], seq[i], post_seq[i], n_mass[i], n_results)]
    result += row_updated_set
  
new_df = pd.DataFrame(data = result, columns = columns_names + ['Prod_sequence', 'Prod_seq_mass', 'Prod_seq_delta_mass'])
new_df.drop_duplicates(subset = ['Prod_sequence'], inplace = True)
new_df.to_csv('result.csv')