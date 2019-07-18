import os
import csv
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet

def get_gene_alias(infile, verbose=0):
    #read the kinase_info.csv file to get the aliases
    df = pd.read_csv(infile, usecols=['Gene', 'Alias'])
    df['Alias'] = df['Alias'].str.split(',')
    df = df.sort_values(by=['Gene'])

    if verbose == 1:
        print(df.keys())
        print('k = {}'.format(df.shape))
        print(df.head(10))

    return df

def get_phospho_table(infile, verbose=0):
    #read the phosphorylation_data.txt file to get the phosphorylation data
    df = pd.read_csv(infile, delimiter='\t')
    #renamed to make table operations easier
    df = df.rename(columns={'Phosphosite': 'Substrate'})
    df['Substrate'] = df['Substrate'].str.strip('\'')
    df = df.sort_values(by=['Substrate'])

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    return df

def get_substrate_table(infile, verbose=0):
    #read the kinase_substrates data.txt file to get the phosphorylation data
    df = pd.read_csv(infile, delimiter='\t', usecols=['Gene', 'Substrate'])
    df = df.sort_values(by=['Gene'])

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    return df

def rename_gene_to_original(sub_table, info_table, verbose=0):
    #for each gene in the sub table, check if the name matches the genes in info table
    #if found then go to the next one
    #if not then look into the aliases of each gene
    #if found then rename gene
    #if not then go to the next one
    COL = list(sub_table)

    df = pd.DataFrame(columns=COL)
    gene_list = pd.DataFrame(columns=['Gene'])
    gene_count = 0

    for index, row in sub_table.iterrows():
        gene = row['Gene']
        if info_table[info_table['Gene']==gene].empty:
            #look into the aliases
            for index, item in info_table.iterrows():
                in_alias = gene in item['Alias']
                if in_alias:
                    if verbose==1:
                        print('Found alias:{} in gene:{}'.format(gene, item['Gene']))

                    row['Gene'] = item['Gene']
                    df.loc[index] = row


                    if gene_list[gene_list['Gene'] == item['Gene']].empty:
                    #if item['Gene'] not in gene_list['Gene']:
                        gene_list.loc[gene_count] = item['Gene']
                        gene_count += 1
                    break
        else:
            if verbose==1:
                print('Found gene:{}'.format(row['Gene']))

            df.loc[index] = row
            if gene_list[gene_list['Gene'] == row['Gene']].empty:
                gene_list.loc[gene_count] = row['Gene']
                gene_count += 1

    if verbose==1:
        print(df.sort_values(by='Gene'))
        print(len(df))
        print(gene_list.shape)
        print(len(gene_list))

    #output table to text file for easier use
    df.to_csv('data/matched_phosphorylation_data.csv', sep='\t', index=False)
    gene_list.to_csv('data/matched_gene_list.csv', sep='\t', index=False)
    #print('donezo')
    x = pd.read_csv('data/matched_gene_list.csv', delimiter='\t')
    print(x)
    #print('donezo again')
    return df, gene_list

def match_p_data_to_substrates(sub_table, phospho_table, verbose=0):
    #merges both tables on the substrate column
    df = pd.merge(sub_table, phospho_table, on='Substrate', how='inner').sort_values(by=['Gene'])

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    return df

def match_gene_list_to_seq(infile, gene_list, verbose=0):
    seq_list = SeqIO.index(infile, 'fasta')
    output = []

    gene = 'Hs_ACTR2.kin_dom'
    seq_rec = seq_list[gene]
    print(seq_rec.seq)

    for gene in gene_list:
        tmp_name = 'Hs_' + gene + '.kin_dom'
        seq_rec = seq_list[tmp_name]
        tmp_rec = SeqRecord(seq_rec.seq,
                            id=gene,
                            name=gene,
                            description='')
        output.append(tmp_rec)

    SeqIO.write(output, 'data/matched_gene_sequences.fasta', 'fasta')

def compute_dist_matrix(infile, verbose=0):
    p_data = pd.read_csv(infile, delimiter='\t')
    gene_list = pd.read_csv('data/matched_gene_list.csv', delimiter='\t')
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

def split_on_equal(str):
    return str.split('=')[1]

def strip_name_on_fasta(str):
    return str.split('_', 1)[-1].split('.', 1)[0]


'''
def get_unique_k_from_substrate_list(infile, verbose=0):
    df = pd.read_csv(infile, delimiter='\t', usecols=['Gene'])
    df = pd.DataFrame.drop_duplicates(df).reset_index().drop('index', axis=1)
    df = df.sort_values(by=['Gene'])

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    return df

def merge_on_gene(df1, df2, verbose):
    df = pd.merge(df1, df2, on='Gene', how='inner')

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    return df

def match_p_data_to_substrates(sub_table, phospho_table, verbose=0):
    count = 0
    #for each phosphosite, look it up in the substrate table
    #if found, return row wich contains 'Gene' 'Substrate' 'Phosphorylation Data'
    
    for index, row in phospho_table.iterrows():
        phosphosite = row['Phosphosite']

        tmp = sub_table[sub_table['Substrate'] == phosphosite]
        if not tmp.empty:
            count+= count
        else:
            print('{} not found in table'.format(phosphosite))
        # for each phosphosite, look it up in the substrate table
        # if found, return row wich contains 'Gene' 'Substrate' 'Phosphorylation Data'
    print(count)

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    return df 
    
def match_sub_data_to_sequence(alias_table, sub_table, verbose=0):
    print('starto')
    count = 0
    df = pd.DataFrame(columns=['Gene', 'Substrate'])
    gene = 'MGC99656'  # row['Gene']
    ne = alias_table[alias_table['Alias'] == gene]
    for i1, row in gene_table.iterrows():
        gene = row['Gene']
        #print(gene)
        a_search = alias_table[alias_table['Gene']==gene]
        s_search = sub_table[sub_table['Gene']==gene]

        if not a_search.empty:
            #merge the rows from substrate table that contain the gene
            print('k found as gene in the table')
            tmp = s_search.groupby('Gene')['Substrate'].apply(' '.join).reset_index()
            df.loc[count] = [tmp['Gene'][0], tmp['Substrate'][0]]
            count += 1
        else:
            #a_search is empty. Look into all the aliases to see if we find matching kinase
            print('k not found as gene, searching in aliases...')
            for i2, entry in alias_table.iterrows():
                if gene in entry['Alias']:
                    s_search = sub_table[sub_table['Gene'] == gene]
                    print('k found as alias in the table')
                    tmp = s_search.groupby('Gene')['Substrate'].apply(' '.join).reset_index()
                    df.loc[count] = [entry['Gene'], tmp['Substrate'][0]]
                    print('name:{} alias:{}'.format(entry['Gene'],tmp['Gene'][0]))
                    count+=1

    df = df.sort_values(by=['Gene']).reset_index().drop(['index'], axis=1)

    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)
        print('{}'.format(count))

    return df.values  
    
def merge_rows_and_columns(p_sub_table, gene_list, verbose=0):
    #for each gene in the gene list, find all ocurrences in p_sub_table and merge
    #the rows
    print(p_sub_table)
    for row in gene_list:
        s_search = p_sub_table[p_sub_table['Gene'] == gene_list]
        tmp = s_search.groupby('Gene')['Substrate'].apply(' '.join)
        print(tmp)
    print(df)
    if verbose == 1:
        print(df.head(10))
        print(df.keys())
        print(df.shape)

    return df
'''


