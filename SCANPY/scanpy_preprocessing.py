
# importing required modules
import matplotlib
import pylab as plt
matplotlib.pyplot.switch_backend('agg')
import os
import glob
import copy
import argparse
import numpy as np
import pandas as pd
import scanpy.api as sc
from itertools import chain

class Single_Cell_Data_Wrangling(object):

    def __init__(self, cell_path_dict, ensembl2symbol):

        self.cell_path_dict = cell_path_dict
        self.ensembl2symbol = ensembl2symbol
        self.cell_batch_names = [batch for batch in self.cell_path_dict.keys()]
        self.cell_dict = {cell_name.split('_')[0]: [] for cell_name in cell_path_dict.keys()}

        for keys, values in self.cell_path_dict.items():
            adata = sc.read(values['filename_data'][0]).transpose()
            adata.var.index = ["-".join(str(hugo).replace("'","-").split("-")[1:-1]) for hugo in np.loadtxt(values['filename_genes'][0], dtype='S')[:, 1]]
            adata.var = adata.var.rename(index=self.ensembl2symbol)
            adata.obs_names = [keys.split('_')[-1] + "_" + "-".join(str(seq).replace("'","-").split("-")[1:-1]) for seq in np.loadtxt(values['filename_barcodes'][0], dtype='S')]
            adata.obs['batch_name'] = keys
            self.cell_dict[keys.split('_')[0]].append({keys.split('_')[-1]:adata})

        # Concatenate each cell batch data set
        self.concatenated_cell_dict = {key: None for key in self.cell_dict.keys()}
        for keys, values in self.cell_dict.items():
                self.concatenated_cell_dict[keys] = values[0][list(values[0].keys())[-1]].concatenate([list(items.values())[0] for items in values[1:]], join='outer')

        # Create a mapping dictionary
        self.batch_mapping_dict = {}
        inc = 0
        for keys in self.concatenated_cell_dict.keys():
            for k in self.concatenated_cell_dict[keys].obs['batch'].unique():
                self.batch_mapping_dict[k] = inc
                inc += 1

        for keys in self.concatenated_cell_dict.keys():
            self.concatenated_cell_dict[keys].obs["batch"].replace(self.batch_mapping_dict, inplace=True)

        self.minimum_cells = 3
        self.minimum_genes = 200

        self.thrsh_mito = 0.2
        self.up_thrsh_genes = 5000
        self.low_thrsh_genes = 50

        self.cell_counts = 1e4


    def filter_cells_and_genes(self):
        cell_and_gene_filtered_dict = copy.deepcopy(self.concatenated_cell_dict)
        for key in cell_and_gene_filtered_dict.keys():
            sc.pp.filter_cells(cell_and_gene_filtered_dict[key], min_genes = self.minimum_cells)
            sc.pp.filter_genes(cell_and_gene_filtered_dict[key], min_cells = self.minimum_genes)
        return cell_and_gene_filtered_dict

    def mitochondria_statistics(self, cell_and_gene_filtered_dict):
        mito_stats_dict = copy.deepcopy(cell_and_gene_filtered_dict)
        for key in mito_stats_dict.keys():
            mito_genes = [name for name in mito_stats_dict[key].var_names if name.startswith('MT.') or name.startswith('MT-')]
            mito_stats_dict[key].obs['percent_mito'] = np.sum(mito_stats_dict[key][:, mito_genes].X, axis=1) / np.sum(mito_stats_dict[key].X, axis=1)
            mito_stats_dict[key].obs['n_counts'] = np.sum(mito_stats_dict[key].X, axis=1)
            sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'percent_mito', title=key, save=key+"_percent_mito_vs_n_counts")
            sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'n_genes', title=key, save=key+"_n_genes_vs_n_count")
        return mito_stats_dict

    def compute_upper_lower_gene_thresholds(self, x, y):
        n_counts_vs_n_genes_pd = pd.DataFrame({'x':x, 'y':y})
        count, bins = np.histogram(n_counts_vs_n_genes_pd)
        up_thrsh_genes = [(x,y) for x,y in zip(count,bins) if x>= 5][-1][-1]
        low_thrsh_genes = [(x,y) for x,y in zip(count,bins)][0][1]
        return up_thrsh_genes, low_thrsh_genes

    def compute_mitochondria_threshold(self, x_, y_):
        n_counts_vs_mito_pct_pd = pd.DataFrame({'x': x_, 'y': y_})
        count, bins = np.histogram(n_counts_vs_mito_pct_pd)
        thrsh_mito = [(x,(y/100000)) for x,y in zip(count, bins) if x>=5][-1][-1]
        return thrsh_mito

    def mitochondria_filtering(self, mito_stats_dict):
        mito_filtered_cell_dict = copy.deepcopy(mito_stats_dict)
        for key in mito_filtered_cell_dict.keys():
            up_thrsh_genes, low_thrsh_genes = self.compute_upper_lower_gene_thresholds(np.array(mito_filtered_cell_dict[key].obs['n_counts']), np.array(mito_filtered_cell_dict[key].obs['n_genes']))
            thrsh_mito = self.compute_mitochondria_threshold(np.array(mito_filtered_cell_dict[key].obs['n_counts']), np.array(mito_filtered_cell_dict[key].obs['percent_mito']))
            mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']< up_thrsh_genes, :]
            mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']> low_thrsh_genes, :]
            mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['percent_mito']< thrsh_mito, :]
        return mito_filtered_cell_dict

    def cell_normalization(self, mito_filtered_cell_dict):
        normalized_mito_filtered_cell_dict = copy.deepcopy(mito_filtered_cell_dict)
        for key in normalized_mito_filtered_cell_dict.keys():
            sc.pp.normalize_per_cell(normalized_mito_filtered_cell_dict[key], counts_per_cell_after=self.cell_counts)
        return normalized_mito_filtered_cell_dict

    def variable_gene_filtering(self, normalized_mito_filtered_cell_dict):
        variable_gene_filtered_cell_dict = copy.deepcopy(normalized_mito_filtered_cell_dict)
        gene_dispersion_dict = {}
        for key in variable_gene_filtered_cell_dict.keys():
            gene_dispersion_dict[key] = sc.pp.filter_genes_dispersion(variable_gene_filtered_cell_dict[key].X, min_mean=0.0125, max_mean=3, min_disp=0.5)
            variable_gene_filtered_cell_dict[key] = variable_gene_filtered_cell_dict[key][:, gene_dispersion_dict[key].gene_subset]
            print('Number of variable genes identified in ' + key + ': ', sum(gene_dispersion_dict[key].gene_subset))
            sc.pl.filter_genes_dispersion(gene_dispersion_dict[key], save=key+"_gene_dispersion_vs mean_expression")
        return variable_gene_filtered_cell_dict, gene_dispersion_dict

    def log_transformation(self, variable_gene_filtered_cell_dict):
        log_filtered_cell_dict = copy.deepcopy(variable_gene_filtered_cell_dict)
        for key in log_filtered_cell_dict.keys():
            sc.pp.log1p(variable_gene_filtered_cell_dict[key].X)
        return log_filtered_cell_dict

    def regress_out(self, log_filtered_cell_dict):
        regress_out_filtered_cell_dict = copy.deepcopy(log_filtered_cell_dict)
        for key in regress_out_filtered_cell_dict.keys():
            sc.pp.regress_out(regress_out_filtered_cell_dict[key], ['n_counts', 'percent_mito'])
        return regress_out_filtered_cell_dict

    def handler(self):
        cell_and_gene_filtered_dict = self.filter_cells_and_genes()
        mito_stats_dict = self.mitochondria_statistics(cell_and_gene_filtered_dict)
        mito_filtered_cell_dict = self.mitochondria_filtering(mito_stats_dict)
        normalized_mito_filtered_cell_dict = self.cell_normalization(mito_filtered_cell_dict)
        variable_gene_filtered_cell_dict, gene_dispersion_dict = self.variable_gene_filtering(normalized_mito_filtered_cell_dict)
        log_filtered_cell_dict = self.log_transformation(variable_gene_filtered_cell_dict)
        regress_out_filtered_cell_dict = self.regress_out(log_filtered_cell_dict)


def generate_gene_mapping_dictionary(args):
    valid_path(args.gene_id_conversion_file[0])
    symbol2ensemble = pd.read_csv(args.gene_id_conversion_file[0], sep='\t', index_col=1)\
['ensembl']
    symbol2ensemble = symbol2ensemble[symbol2ensemble.notnull()]
    ensembl2symbol = dict(zip(symbol2ensemble.values,symbol2ensemble.index.values))
    return ensembl2symbol

def generate_cell_batch_dictionary(args):
    path_prefix = "/*/outs/filtered_gene_bc_matrices"
    if glob.glob(args.matrix_file[0]+path_prefix):
        cell_path_dict = {}
        for matrix_path in glob.glob(args.matrix_file[0]+path_prefix):
            cell_path_dict[matrix_path.split("/")[-3]] = {'filename_data': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1].split('.')[-1] == 'mtx'], \
            'filename_genes': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1] == 'genes.tsv'], \
            'filename_barcodes': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1] == 'barcodes.tsv']}
        return cell_path_dict
    else:
        print(INVALID_DIRECTORY_MSG%(args.matrix_file[0]))


# error messages
INVALID_PATH_MSG = "Error: Invalid file path/name. Path %s does not exist."
INVALID_DIRECTORY_MSG = "Error: Invalid prefix path/name. Path %s does not contain cellRanger output directories containing a metrics_summary.csv file."

def validate_file(file_name):
    '''
    validate file name and path.
    '''
    if not valid_path(file_name):
        print(INVALID_PATH_MSG%(file_name))
        quit()
    return


def valid_path(path):
    # validate file path
    return os.path.exists(path)


def main():
    # create a parser object
    parser = argparse.ArgumentParser(description = "Preprocess a sc-RNASeq run using the cellRanger filtered gene-barcode matrices containing only cellular barcodes in MEX format.")

    # defining arguments for parser object
    parser.add_argument("-m", "--matrix_file", type = str, nargs = 1, default = None,
                        help = "Path to the cellRanger output directory/directories.")

    parser.add_argument("-c", "--gene_id_conversion_file", type = str, nargs = 1, default = None,
                        help = "Path to the Ensemble to Hugo gene ID conversion file.")

    # parse the arguments from standard input
    args = parser.parse_args()

    if args.matrix_file != None and args.gene_id_conversion_file != None:
        execute = Single_Cell_Data_Wrangling(generate_cell_batch_dictionary(args), generate_gene_mapping_dictionary(args))
        execute.handler()


if __name__ == '__main__':
    main()
