import pandas as pd
from table_methods import merge_dataframes, split_column, merge_columns, strip_last_char, split_float_column
import numpy as np

def merge_sub_site(ksa_path, output_folder, delimiter):
    #merges the substrate and the site to generate the
    #phosphosite name from KSA_human.txt file
    df = pd.read_csv(ksa_path, delimiter=delimiter)
    tmp = merge_columns(df, ['Substrate', 'Site'], 'Phosphosite', '-')
    tmp[['Kinase', 'Phosphosite']].to_csv(output_folder + '/Phosphosite.csv', sep='\t', index=False)

def output_breastc_to_csv(phospho_path, output_folder, delimiter):
    df = pd.read_csv(phospho_path, delimiter=delimiter)
    tmp = merge_columns(df, ['Gene', 'Site'], 'Phosphosite', '-')
    tmp.dataframe[2:].to_csv(output_folder + '/Phosphorylation.csv', sep='\t', index=False)

def remove_rows_with_multisites(dataframe, sites_column, delimiter):
    df = dataframe.copy()
    split_column(df, sites_column, delimiter)
    df = df[df[sites_column].map(len) == 1]

    return df

def remove_rows_with_multiposition(dataframe, sites_column, delimiter):
    df = dataframe.copy()
    split_float_column(df, sites_column, delimiter)

    for index, row in df.iterrows():
        if len(row[sites_column]) != 1:
            df.drop(index, inplace=True)

    df[sites_column] = df[sites_column].apply(lambda row: str(row[0]))
    df = df.reset_index(drop=True)

    return df

def list_length(list):
    return len(list)

def remove_rows_with_nan_on_col(dataframe, column):
    df = dataframe.copy()
    df.dropna(subset=column, how='all', inplace = True)
    return df

def make_phosphosite_column(dataframe, gene_column, site_column):
    df = dataframe.copy()
    df = merge_columns(df, [gene_column, site_column], 'Phosphosite', '-')
    return df

def delete_rows_with_fewer_samples(dataframe):
    df = dataframe.copy()
    df = df.reset_index(drop=True)
    n_col = len(df.columns) - 1
    Nan_list = df.isnull().sum(axis=1).divide(n_col)

    for index, row in df.iterrows():
        if Nan_list[index] > .5:
            df.drop(index, inplace=True)

    return df

def average_nans(dataframe):
    df = dataframe.copy()
    df = df.reset_index(drop=True)

    for index, row in df.iterrows():
        avg = row[2:].mean()
        row = row.fillna(avg)
        df.iloc[index] = row

    return df

def remove_k_less_than_2(phospho_data, gene_list):
    df = phospho_data.copy()
    df = df.reset_index(drop=True)

    for index, row in gene_list.iterrows():
        tmp = df[df['Kinase'] == row['Kinase']]

        if len(df[df['Kinase'] == row['Kinase']].index) < 2:
            gene_list.drop(index, inplace=True)
            df.drop(tmp.index, inplace=True)

    return df, gene_list

def get_substrate_mean(phospho_data, gene_list):
    df = pd.DataFrame(columns=phospho_data.columns.values)

    for index, row in gene_list.iterrows():
        tmp = phospho_data[phospho_data['Kinase'] == row['Kinase']]
        new_row = tmp.groupby('Kinase').mean()
        new_row.reset_index(level=0, inplace=True)
        df = df.append(new_row)

    df = df.reset_index(drop=True)

    return df
