from aliases import rename_gene_to_original
from table_methods import merge_dataframes, strip_last_char
from distance_matrix import compute_dist_matrix_from_file
from sequence_match import match_seq_to_gene_list
from excel import export_xl_page_to_csv
import pandas as pd
import data_prep as dp

'''
1.  Export page from Excel file
2.  Read csv file
3.  Remove unnecessary columns
4.  Rename site to phosphosite
5.  Load phosphosite data from ksa
6.  Remove rows with nan values on the phosphosite column
7.  Strip last character from sites column
8.  match phosphosites to substrate table
9.  Delete rows with fewer samples
10. Average nan values on each row
11. Check for aliases
11. Remove genes with less than 2 substrates
12. Match to sequence data
13. Compute distance matrix
14. Prepare table for Kmeans
'''

def main():
    print('Starto')
    folder = 'ovarian_cancer_data'
    xl_path = 'ddata/ovarian_cancer.xlsx'
    xl_output = 'Phosphorylation'
    xl_page = 'FullOvarianPhosphoproteomeData'
    ksa_path = 'ddata/KSA_human.txt'
    sub_path = folder + '/Phosphosite.csv'
    info_path = 'data/kinase_info.csv'
    phos_path = folder + '/Phosphorylation.csv'
    delimiter = '\t'

    '''
    print('preparing data from xl file...')
    export_xl_page_to_csv(xl_path=xl_path, output_name=xl_output, output_folder=folder, page_name=xl_page)
    '''

    #TODO: add drop columns method
    print('Read csv file')
    phos_data = pd.read_csv(filepath_or_buffer=phos_path, delimiter=delimiter)

    print('droping unnecessary columns from phosphorylation data...')
    rcol = ['Gene', 'RefSeq', 'Peptide', 'AScores']
    phos_data = phos_data.drop(rcol, axis=1)

    print('Renaming site column...')
    phos_data = phos_data.rename(columns={'Site': 'Phosphosite'})

    print('loading phosphosite data from ksa...')
    sub_data = pd.read_csv(filepath_or_buffer=sub_path, delimiter=delimiter)

    print('removing rows with nan values on the phosphosite column...')
    phos_data = dp.remove_rows_with_nan_on_col(dataframe=phos_data, column=['Phosphosite'])

    print('striping last character from sites column...')
    phos_data = strip_last_char(dataframe=phos_data, column='Phosphosite')

    print('match phosphosites to substrate table...')
    merged_data = merge_dataframes(columns='Phosphosite', df1=phos_data, df2=sub_data)

    print('deleting rows with fewer samples...')
    merged_data = dp.delete_rows_with_fewer_samples(merged_data)

    print('averaging nan values...')
    merged_data = dp.average_nans(merged_data)

    print('checking on info table for kinase aliases...')
    alias_data = pd.read_csv(info_path, delimiter=',')
    alias_data = alias_data.loc[:, ~alias_data.columns.str.contains('^Unnamed')]
    df, gene_list = rename_gene_to_original(merged_data, alias_data, folder, 0)

    print('removing kinases with less than 2 substrates...')
    df, gene_list = dp.remove_k_less_than_2(df, gene_list)
    df.to_csv(folder + '/matched_phosphorylation.csv', sep='\t', index=False)
    gene_list.to_csv(folder + '/matched_gene_list.csv', sep='\t', index=False)

    print('matching to sequencing data...')
    match_seq_to_gene_list('data/human_kinase_domain.fasta', gene_list, folder)

    print('computing distance matrix...')
    d = compute_dist_matrix_from_file(df, gene_list, folder)

    print('ready to compute hierarchical clustering from distance matrix...')

    print('preparing data for kmeans...')
    k_phospho_data = df.drop(['Phosphosite'], axis=1)
    k_data = dp.get_substrate_mean(phospho_data=k_phospho_data, gene_list=gene_list)
    k_data.to_csv(folder + '/oc_kmeans_table.csv', sep=delimiter, index=False)

if __name__ == "__main__":
    main()