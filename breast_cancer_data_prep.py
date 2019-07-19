from table import table
from aliases import rename_gene_to_original
from table_methods import merge_tables, split_column, strip_column, merge_columns
from distance_matrix import compute_dist_matrix_from_file
from sequence_match import match_seq_to_gene_list
from excel import export_xl_page_to_csv
import re
import pandas as pd

#TODO: add notes about the steps of preparing the data
#TODO: remove unnecessary functions

def main():
    print('Starto')
    sub_path = 'breast_cancer_data/Phosphosite.csv'
    info_path = 'data/kinase_info.csv'
    bc_path = 'breast_cancer_data/BreastCancerData.csv'
    phospho_path = 'breast_cancer_data/Phosphorylation_data.csv'
    p_path = 'breast_cancer_data/Phosphorylation.csv'
    mp_path = 'breast_cancer_data/matched_phosphorylation_data.csv'
    mps_path = 'breast_cancer_data/matched_phosphorylation_samples.csv'
    rnan_path = 'breast_cancer_data/matched_phosphorylation_rnan.csv'
    delimiter = '\t'

    #export_xl_page_to_csv()
    #output_phoshosite_csv(sub_path, 'breast_cancer_data', delimiter)
    #remove_rows_with_multisites(bc_path, 'breast_cancer_data', delimiter)
    #clean_sites_col(phospho_path, 'breast_cancer_data', delimiter)
    #join_gene_site(p_path, 'breast_cancer_data', delimiter)
    #phospho_table = table()
    #phospho_table.initialize_from_file(path=phospho_path, delimiter=delimiter)
    #print(phospho_table.dataframe.iloc[:3])
    #sub_table = table()
    #sub_table.initialize_from_file(path=sub_path, delimiter=delimiter)
    #print(sub_table.dataframe.iloc[:3])
    #joined_table = merge_tables('Phosphosite', sub_table, phospho_table, 'breast_cancer_data')
    #joined_table.sort('Kinase', True)
    #print(joined_table.dataframe)
    #Nan_list = joined_table.dataframe.isnull().sum(axis=1).divide(27)
    #print(Nan_list)
    #delete_rows_with_fewer_samples(mp_path, 'breast_cancer_data', delimiter)
    #average_nans(mps_path, 'breast_cancer_data', delimiter)
    #tmp1 = pd.read_csv(rnan_path, delimiter=delimiter)
    #merged_table = table()
    #merged_table.initialize_from_dataframe(tmp1)
    #tmp2 = pd.read_csv(info_path, delimiter=',')
    #info_table = table()
    #info_table.initialize_from_dataframe(tmp2)
    #df, gene_list = rename_gene_to_original(merged_table, info_table, 'breast_cancer_data', 0)
    #gene_list = pd.read_csv('breast_cancer_data/matched_gene_list.csv')
    #match_seq_to_gene_list('data/human_kinase_domain.fasta', gene_list, 'breast_cancer_data')
    d = compute_dist_matrix_from_file('breast_cancer_data')

def output_phoshosite_csv(sub_path, output_folder, delimiter):
    df = pd.read_csv(sub_path, delimiter=delimiter)
    tmp = merge_columns(df, ['Substrate', 'Site'], 'Phosphosite', '-')
    tmp.dataframe[['Kinase', 'Phosphosite']].to_csv(output_folder + '/Phosphosite.csv', sep='\t', index=False)

def output_breastc_to_csv(phospho_path, output_folder, delimiter):
    df = pd.read_csv(phospho_path, delimiter=delimiter)
    tmp = merge_columns(df, ['Gene', 'Site'], 'Phosphosite', '-')
    tmp.dataframe[2:].to_csv(output_folder + '/Phosphorylation.csv', sep='\t', index=False)

def remove_rows_with_multisites(phospho_path, output_folder, delimiter):
    df = pd.read_csv(phospho_path, delimiter=delimiter)
    df2 = df.copy()
    split_column(df2, 'variableSites', ' ')

    for index, row in df2.iterrows():
        tmp = [x for x in row['variableSites'] if x != '']
        #print('{}\n{}'.format(df.iloc[index]['variableSites'], row['variableSites']))

        if len(tmp) != 1:
            #print('{}\n{}\n{}'.format(index, row['variableSites'], len(row['variableSites'])))
            df.drop(index, inplace=True)

    df.to_csv(output_folder + '/Phosphorylation.csv', sep='\t', index=False)

def clean_sites_col(phospho_path, output_folder, delimiter):
    #to be called after remove_rows_with_multisites function
    df = pd.read_csv(phospho_path, delimiter=delimiter)
    df['variableSites'] = df['variableSites'].apply(remove_char_list)
    df.to_csv(output_folder + '/Phosphorylation.csv', sep='\t', index=False)

def remove_char_list(s):
    char_list = ['\'', '[', ']', ',', ' ']
    return re.sub("|".join(char_list), "", s)

def to_str(list):
    return ''.join(list)

def join_gene_site(phospho_path, output_folder, delimiter):
    df = pd.read_csv(phospho_path, delimiter=delimiter)
    strip_column(df, 'variableSites', ' ')
    strip_column(df, 'variableSites', 's')
    strip_column(df, 'variableSites', 't')
    strip_column(df, 'variableSites', 'y')
    strip_column(df, 'variableSites', 'q')
    strip_column(df, 'variableSites', 'm')
    print(df['variableSites'])
    tmp = merge_columns(df, ['geneSymbol','variableSites'], 'Phosphosite', '-')
    print(tmp.dataframe['Phosphosite'])
    colnames = tmp.dataframe.columns.tolist()
    print(colnames)
    colnames = colnames[-1:] + colnames[:-1]
    print(colnames)
    tmp.dataframe = tmp.dataframe[colnames]
    tmp.dataframe = tmp.dataframe.drop(['geneSymbol', 'variableSites'], axis=1)
    print(tmp.dataframe)
    tmp.dataframe.to_csv(output_folder + '/Phosphorylation_data.csv', sep='\t', index=False)

def delete_rows_with_fewer_samples(phospho_path, output_folder, delimiter):
    df = pd.read_csv(phospho_path, delimiter=delimiter)
    Nan_list = df.isnull().sum(axis=1).divide(27)

    for index, row in df.iterrows():
        if Nan_list[index] > .5:
            #print('{}\n{}'.format(df.iloc[index]['Kinase'], row))
            df.drop(index, inplace=True)

    df.to_csv(output_folder + '/matched_phosphorylation_samples.csv', sep='\t', index=False)

def average_nans(phospho_path, output_folder, delimiter):
    df = pd.read_csv(phospho_path, delimiter=delimiter)

    for index, row in df.iterrows():
        avg = row[2:].mean()
        row = row.fillna(avg)
        df.iloc[index] = row

    print(df)

    df.to_csv(output_folder + '/matched_phosphorylation_rnan.csv', sep='\t', index=False)


if __name__ == "__main__":
    main()

