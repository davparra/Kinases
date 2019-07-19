import pandas as pd
import numpy as np
#TODO:Fix column names
def compute_dist_matrix_from_file(input_folder, verbose=0):
    p_data = pd.read_csv(input_folder + '/matched_sub_data.csv', delimiter='\t')
    gene_list = pd.read_csv(input_folder + '/matched_gene_list.csv', delimiter='\t')
    col = list(gene_list['Kinase'])
    df = pd.DataFrame(columns=col, index=col)

    for index, row in gene_list.iterrows():
        data1 = p_data[p_data['Kinase'] == row['Kinase']]
        data1 = data1[p_data.columns[2:]].mean()

        for j, r in gene_list.iterrows():
            if r['Kinase'] == row['Kinase']:
                corr = 0
                df.loc[row['Kinase'], r['Kinase']] = corr

            else:
                data2 = p_data[p_data['Kinase'] == r['Kinase']]
                data2 = data2[p_data.columns[2:]].mean()
                corr = np.square(np.triu(np.corrcoef(data1, data2)))[0, 1]
                df.loc[row['Kinase'], r['Kinase']] = corr

    print(df)
    df = df.multiply(1000)
    df = df.astype(int)

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    df.to_csv(input_folder + '/matched_phosphorylation_dist_matrix.csv', sep=',', index=True)

    return df