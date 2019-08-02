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

    print('loading phosphorylation data to prepare for kmeans...')
    phospho_data = pd.read_csv(filepath_or_buffer=phos_path, delimiter=delimiter)
    phospho_data = phospho_data.drop(['Phosphosite'], axis=1)
    gene_list = pd.read_csv(filepath_or_buffer=gene_list_path, delimiter=delimiter)

    km_data = dp.get_substrate_mean(phospho_data=phospho_data, gene_list=gene_list)


if __name__ == "__main__":
    main()