import os
import pandas

this_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
csv_path = os.path.join(this_dir, "enzymes_data", "enzymes_infos.csv")
dataframe = pandas.read_csv(csv_path, index_col='enzyme', sep=';')
enzymes_infos = dataframe.T.to_dict()