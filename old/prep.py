#import data_prep as dp
from table import table
import pandas as pd

def main():
    print('Starto')
    infile = 'data/Kinase_substrates.txt'
    delimiter = '\t'
    #info_table = dp.get_gene_alias('data/kinase_info.csv', 0)
    #sub_table = dp.get_substrate_table('data/Kinase_Substrates.txt', 0)
    #phospho_table = dp.get_phospho_table('data/phosphorylation_data.txt', 0)
    #t = dp.match_p_data_to_substrates(sub_table, phospho_table, 0)
    #n, gene_list = dp.rename_gene_to_original(t, info_table, 1)
    #dp.match_gene_list_to_seq('data/human_kinase_domain.fasta', gene_list, 0)
    #dp.compute_dist_matrix('data/matched_phosphorylation_data.csv', 1)
    substrate_table = table()
    print(substrate_table.columns())
    substrate_table.initialize_from_file(infile, '\t')
    substrate_table.sort(['Gene', 'Substrate'], True)
    #print(substrate_table)

    df = pd.ExcelFile('xlsx/BreastCancerData.xlsx')
    print(df)


if __name__ == "__main__":
    main()