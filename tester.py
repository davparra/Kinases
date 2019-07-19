from table import table
from aliases import rename_gene_to_original
from table_methods import merge_tables, split_column, strip_column
from distance_matrix import compute_dist_matrix_from_file
from sequence_match import match_seq_to_gene_list
import pandas as pd
#algorithm
#1. Read phospho data and

def main():
    print('Starto')
    sub_path = 'data/Kinase_substrates.txt'
    info_path = 'data/kinase_info.csv'
    phospho_path = 'data/phosphorylation_data.txt'
    delimiter = '\t'
    #info_table = dp.get_gene_alias('data/kinase_info.csv', 0)
    #sub_table = dp.get_substrate_table('data/Kinase_Substrates.txt', 0)
    #phospho_table = dp.get_phospho_table('data/phosphorylation_data.txt', 0)
    #t = dp.match_p_data_to_substrates(sub_table, phospho_table, 0)
    #n, gene_list = dp.rename_gene_to_original(t, info_table, 1)
    #dp.match_gene_list_to_seq('data/human_kinase_domain.fasta', gene_list, 0)
    #dp.compute_dist_matrix('data/matched_phosphorylation_data.csv', 1)
    #substrate_table = table()
    #print('Dataframe is empty: {}'.format(substrate_table.empty()))
    #print('Columns: {}'.format(substrate_table.columns()))
    #substrate_table.initialize_from_file(sub_path, delimiter=delimiter)
    #print('Dataframe is empty: {}'.format(substrate_table.empty()))
    #substrate_table.sort(['Gene'], True)
    #phospho_table = table().initialize_from_file(phospho_path, delimiter=delimiter)
    #strip_column(phospho_table.dataframe, 'Substrate', '\'')
    #print(phospho_table.length)
    #info_table = table().initialize_from_file(path=info_path, delimiter=',', columns=['Gene', 'Alias'])
    #split_column(info_table.dataframe, 'Alias', ',')
    #info_table.sort(['Gene'], True)

    #print(info_table)
    #print(info_table.columns())

    #merged_table = merge_tables('Substrate', substrate_table, phospho_table)
    #merged_table.sort(['Gene'], True)
    #print(merged_table)

    #df, gene_list = rename_gene_to_original(merged_table, info_table, 'phos_data', 0)

    #m = compute_dist_matrix_from_file('phos_data', 0)

    #df = pd.read_excel('xlsx/BreastCancerData.xlsx', sheet_name='data')
    #print(df.columns.values)
    gene_list = pd.read_csv('phos_data/matched_gene_list.csv', delimiter='\t')
    print(gene_list)
    tmp = match_seq_to_gene_list('data/human_kinase_domain.fasta', gene_list, 'phos_data')

if __name__ == "__main__":
    main()