import glob
import pickle
import openpyxl
import argparse
import pandas as pd

class Generate_Preprocessing_Table(object):

    def __init__(self, summary_pickle):
        self.output_dir = summary_pickle
        self.pickle_dict = {}
        for pickle_file in glob.glob(summary_pickle+'/*.pkl'):
            with open(pickle_file, 'rb') as handle:
                data = pickle.load(handle)
                self.pickle_dict[pickle_file.split('/')[-1].split('.')[0]] = data

    def generate_data_frame(self):
        columns = set()
        index = []

        for seq_run, batches in self.pickle_dict.items():
            for batch_key in batches.keys():
                if len(batch_key.split('_')) < 2:
                    index.append(batch_key)
                    for data_key in batches[batch_key].keys():
                        columns.add(data_key)
        return pd.DataFrame(index=index, columns=list(columns))

    def populate_data_frame(self, report_df):
        for batch, data in self.pickle_dict.items():
            for items in data.keys():
                if len(items.split('_')) < 2:
                    for keys in data[items].keys():
                        report_df[keys][items] = self.pickle_dict[batch][items][keys]
        return report_df

    def export_excel_report(self):
        empty_report_df = self.generate_data_frame()
        report_df = self.populate_data_frame(empty_report_df)
        for batch, data in self.pickle_dict.items():
            for items in data.keys():
                if len(items.split('_')) < 2:
                    for keys in data[items].keys():
                        report_df[keys][items] = self.pickle_dict[batch][items][keys]
        output_prefix = self.output_dir.split('/')[0]
        report_df.to_excel(self.output_dir + output_prefix + "_preprocessing_summary_table.xlsx")


def main():
    # create a parser object
    parser = argparse.ArgumentParser(description = "Generate preprocessing summary table.")

    parser.add_argument("--summary_pickle", type = str, nargs = 1,
                        help = "Path to the pre-processed summary output pickle file.")

    args = parser.parse_args()
    if args.summary_pickle[0] != None:
        execute = Generate_Preprocessing_Table(args.summary_pickle[0])
        execute.export_excel_report()

if __name__ == '__main__':
    main()
