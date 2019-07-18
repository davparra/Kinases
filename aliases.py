#TODO: optimize code
import pandas as pd

def rename_gene_to_original(substrates, aliases, verbose=0):
    #for each gene in the sub table, check if the name matches the genes in info table
    #if found then go to the next one
    #if not then look into the aliases of each gene
    # if found then rename gene
    #if not then go to the next one
    COL = list(substrates.dataframe)

    df = pd.DataFrame(columns=COL)
    gene_list = pd.DataFrame(columns=['Gene'])
    gene_count = 0

    for index, row in substrates.dataframe.iterrows():
        gene = row['Gene']
        if aliases.dataframe[aliases.dataframe['Gene'] == gene].empty:
            #look into the aliases
            for index, item in aliases.dataframe.iterrows():
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