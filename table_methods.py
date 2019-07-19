from table import table
import pandas as pd

#TODO: Think about what other functionality would be useful
#TODO: Clean up code
def merge_tables(columns, table1, table2, output_folder, join='inner'):
    #merge 2 dataframes based on a column or columns
    #join method can be specified: {‘left’, ‘right’, ‘outer’, ‘inner’}
    #default join value is ‘inner’
    #check that both tables are not empty
    tmp = pd.DataFrame()
    table1.dataframe[columns] = table1.dataframe[columns].astype(str)
    table2.dataframe[columns] = table2.dataframe[columns].astype(str)

    print(table2.dataframe[columns].dtype)

    if not table1.dataframe.empty and not table2.dataframe.empty:
        try:
            tmp = pd.merge(table1.dataframe, table2.dataframe, on=columns, how=join)
        except:
            print('An exception ocurred.')
    else:
        print('A table is empty')

    new_table = table().initialize_from_dataframe(tmp)

    new_table.dataframe.to_csv(output_folder + '/matched_phosphorylation_data.csv', sep='\t', index=False)


    return new_table

def merge_columns(df, columns, col_name, sep):
    #returns a new table that is a copy of the old table plus a column
    #produced by the columns merged
    print(df[columns].values)
    df[col_name] = df[columns].apply(lambda x: sep.join(x), axis=1)
    print(df[col_name])
    return table().initialize_from_dataframe(df)

'''
def merge_rows():
    pass

def match_tables(columns):
    #returns a table containing the rows of that matched on the column
    pass
'''



def strip_column(dataframe, column, char):
    #removes character char from all entries in column
    dataframe[column] = dataframe[column].str.strip(char)

def split_column(dataframe, column, char):
    #splits the entries of a column if the entry is a list divided by a character char
    dataframe[column] = dataframe[column].str.split(char)
    return dataframe