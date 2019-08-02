import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

#TODO: Think about what other functionality would be useful
#TODO: Clean up code
#TODO: drop_columns function
def merge_dataframes(columns, df1, df2, join='inner'):
    #merge 2 dataframes based on a column or columns
    #join method can be specified: {‘left’, ‘right’, ‘outer’, ‘inner’}
    #default join value is ‘inner’
    #check that both tables are not empty
    new_dataframe = pd.DataFrame()
    df1[columns] = df1[columns].astype(str)
    df2[columns] = df2[columns].astype(str)

    if not df1.empty and not df2.empty:
        try:
            new_dataframe = pd.merge(df1, df2, on=columns, how=join)
        except:
            print('An exception ocurred.')
    else:
        print('A table is empty')

    colnames = new_dataframe.columns.tolist()
    colnames = colnames[-1:] + colnames[:-1]
    new_dataframe = new_dataframe[colnames]

    return new_dataframe
    #new_dataframe.to_csv(output_folder + '/matched_phosphorylation_data.csv', sep='\t', index=False)

def merge_columns(dataframe, columns, col_name, sep):
    #returns a new table that is a copy of the old table plus a column
    #produced by the columns merged
    dataframe[col_name] = dataframe[columns].apply(lambda x: sep.join(x), axis=1)
    dataframe = dataframe.drop(columns, axis=1)
    colnames = dataframe.columns.tolist()
    colnames = colnames[-1:] + colnames[:-1]
    dataframe = dataframe[colnames]

    return dataframe

'''
def merge_rows_by_avg(dataframe, oncolumn, withcolumn):
    df = dataframe.groupby(oncolumn)[withcolumn].apply(' '.join).reset_index()
    print(df.columns)
    pass
'''

def strip_column(dataframe, column, char):
    #removes character char from all entries in column
    dataframe[column] = dataframe[column].str.strip(char)

def strip_last_char(dataframe, column):
    #TODO: fix parsing of col: variableSites
    #removes last character from all entries in a column
    dataframe[column] = dataframe[column].apply(lambda row: row[:-1])

    return dataframe

def split_column(dataframe, column, char):
    #splits the entries of a column if the entry is a list divided by a character char
    dataframe[column] = dataframe[column].str.split(char)
    #remove empty items from list
    dataframe[column] = dataframe[column].apply(lambda row: [x for x in row if x != ''])

    return dataframe