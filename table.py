import pandas as pd
import os

class table:
    def __init__(self):
        self.dataframe = pd.DataFrame()
        self.length = 0
        self.shape = None

    def __repr__(self):
        return self.dataframe.to_string()

    def initialize_from_file(self, path, delimiter, columns=None):
        #initializes the table from an input file
        if os.path.exists(path):
            try:
                self.dataframe = pd.read_csv(path, delimiter=delimiter)
                self.length = len(self.dataframe)
                self.shape = self.dataframe.shape

                if columns:
                    self.dataframe = self.dataframe[columns]

                return self
            except:
                print('An exception ocurred.')
        else:
            print('Check that the path to the file is correct.')


    def initialize_from_dataframe(self, dataframe):
        #initializes the table from a given dataframe object
        if not dataframe.empty:
            try:
                self.dataframe = dataframe
                self.length = len(dataframe)
                self.shape = dataframe.shape

                return self
            except:
                print('An exception ocurred.')
        else:
            print('Dataframe input is empty.')

    def sort(self, column, order):
        #sorts the table using a given column or columns
        #order is bool, True for ascending order, False for descending order
        self.dataframe = self.dataframe.sort_values(by=column, ascending=order)

    def columns(self):
        #returns an array of the columns in the table
        return self.dataframe.columns.values

    def empty(self):
        #returns wether the dataframe of the table is empty
        return self.dataframe.empty