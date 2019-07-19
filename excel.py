import pandas as pd

def export_xl_page_to_csv(file_path, output_name, output_folder, page_name, col_converter=None, converter=None):
    df = pd.read_excel(file_path, sheet_name=page_name, parse_dates=False)
    print(df)
    df.to_csv(output_folder + '/' + output_name + '.csv', sep='\t', index=False)