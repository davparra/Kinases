from aliases import rename_gene_to_original
from table_methods import merge_dataframes, strip_last_char
from distance_matrix import compute_dist_matrix_from_file
from sequence_match import match_seq_to_gene_list
from excel import export_xl_page_to_csv
import pandas as pd
import data_prep as dp

def main():
    print('Starto')
    folder = 'breast_cancer_data'

    phos_path = folder + '/matched_phosphorylation.csv'
    gene_list_path = folder + '/matched_gene_list.csv'
    delimiter = '\t'

    l = [.1298]
    print(type(str(l[0])))

if __name__ == "__main__":
    main()