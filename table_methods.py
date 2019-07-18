from table import table
import pandas as pd

#TODO: Think about what other functionality would be useful
def merge_tables(columns, table1, table2, join='inner'):
    #merge 2 dataframes based on a column or columns
    #join method can be specified: {‘left’, ‘right’, ‘outer’, ‘inner’}
    #default join value is ‘inner’
    #check that both tables are not empty
    tmp = pd.DataFrame()

    if not table1.dataframe.empty and not table2.dataframe.empty:
        try:
            tmp = pd.merge(table1.dataframe, table2.dataframe, on=columns, how=join)
        except:
            print('An exception ocurred.')
    else:
        print('A table is empty')

    new_table = table().initialize_from_dataframe(tmp)

    return new_table


def merge_rows():
    pass
'''
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