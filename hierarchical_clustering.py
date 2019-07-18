import pandas as pd
import numpy as np
#TODO: implement hiearchical clustering necessary functions

def compute_dist_matrix(data_path, gene_path, verbose=0):
    p_data = pd.read_csv(data_path, delimiter='\t')
    gene_list = pd.read_csv(gene_path, delimiter='\t')
    col = list(gene_list['Gene'])
    df = pd.DataFrame(columns=col, index=col)
    print(df.shape)
    count = 0
    #print(gene_list['Gene'])

    for index, row in gene_list.iterrows():
        data1 = p_data[p_data['Gene'] == row['Gene']]
        data1 = data1[p_data.columns[2:20]].mean()

        for j, r in gene_list.iterrows():
            if r['Gene'] == row['Gene']:
                corr = 0
                df.loc[row['Gene'], r['Gene']] = corr

            else:
                data2 = p_data[p_data['Gene'] == r['Gene']]
                data2 = data2[p_data.columns[2:20]].mean()
                corr = np.square(np.triu(np.corrcoef(data1, data2)))[0, 1]
                df.loc[row['Gene'], r['Gene']] = corr

    df = df.multiply(1000)
    df = df.astype(int)

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    df.to_csv('data/matched_phosphorylation_dist_matrix.csv', sep=',', index=True)

    return df