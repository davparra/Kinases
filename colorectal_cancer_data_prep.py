from aliases import rename_gene_to_original
from table_methods import merge_dataframes, strip_last_char, lowercase_column, merge_columns
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
    folder = 'colorectal_cancer_data'
    xl_path = 'ddata/colorectal_cancer.xlsx'
    xl_output = 'Phosphorylation'
    xl_page = 'Class1 pSTY site (13411)'
    ksa_path = 'ddata/KSA_human.txt'
    sub_path = folder + '/Phosphosite.csv'
    info_path = 'data/kinase_info.csv'
    phos_path = folder + '/Phosphorylation.csv'
    delimiter = '\t'
    '''
    print('preparing data from xl file...')
    export_xl_page_to_csv(xl_path=xl_path, output_name=xl_output, output_folder=folder, page_name=xl_page)
    '''
    print('Read csv file')
    phos_data = pd.read_csv(filepath_or_buffer=phos_path, delimiter=delimiter)
    #print(phos_data.columns)

    print('droping unnecessary columns from phosphorylation data...')
    rcol = ['Protein names', 'Leading proteins (Uniprot ID)', 'Delta score',
       'Localization prob E1', 'Localization prob E2', 'Localization prob E3',
       'Assigned in PhosphositePlus database']

    phos_data = phos_data.drop(rcol, axis=1)

    print('Switch aminoacid column to lowercase')
    #print(phos_data['Amino acid'])
    phos_data = lowercase_column(phos_data, ['Amino acid'])

    print('Removing rows with nan on gene')
    phos_data = dp.remove_rows_with_nan_on_col(dataframe=phos_data, column=['Gene names'])

    print('Removing rows with nan on position')
    phos_data = dp.remove_rows_with_nan_on_col(dataframe=phos_data, column=['Positions within proteins'])

    print('Removing rows with multiple positions')
    phos_data = dp.remove_rows_with_multiposition(phos_data, 'Positions within proteins', ';')

    print('Merging position with aminoacid columns')
    phos_data = merge_columns(phos_data, ['Positions within proteins', 'Amino acid'], 'p+a','')

    print('Mergeing p+a column with gene column')
    phos_data = dp.make_phosphosite_column(phos_data, 'Gene names', 'p+a')

    print(phos_data['Phosphosite'])
    print(phos_data.columns.values)
    print()

    #TODO: Match data to ksa and sequences


if __name__ == "__main__":
    main()