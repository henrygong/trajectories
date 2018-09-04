
# importing required modules
import matplotlib
import pylab as plt
matplotlib.pyplot.switch_backend('agg')
import os
import glob
import copy
import argparse
import pickle
import numpy as np
import pandas as pd
import scanpy.api as sc
# from itertools import chain

class Single_Cell_Data_Wrangling(object):
    def __init__(self, cell_path_dict, ensembl2symbol,\
                minimum_cells, minimum_genes, counts_per_cell_after, \
                thrsh_mito, up_thrsh_genes, low_thrsh_genes, output_dir,\
                load_h5ad, output_unprocessed_h5ad, only_export_unprocessed_h5ad,\
                marker_genes):

        self.output_dir = output_dir
        try:
            if not os.path.exists(self.output_dir[0]):
                os.makedirs(self.output_dir[0])
        except OSError:
            print("Error creating directory.")

        self.load_h5ad = load_h5ad

        if self.load_h5ad:
            # Update
            print("1. Loading single cell data from h5 file.")
            print("\n")
            print("2. Skipping concatenation step.")
            print("\n")
            if glob.glob(self.load_h5ad[0]+"/*unprocessed.h5ad"):
                self.concatenated_cell_dict = {}
                if glob.glob(self.load_h5ad[0]+"/cache"):
                    for items in glob.glob(self.load_h5ad[0]+"/*unprocessed.h5ad"):
                        print(items.split("/")[-1])
                        print(items.split("/")[-1].split('.')[0])

                        self.concatenated_cell_dict[items.split("/")[-1].split('.')[0]\
                        .split("_")[0]] = sc.read(self.load_h5ad[0]\
                        +"/cache/"+items.split("/")[-1].split('.')[0], cache=True)
                else:
                    for items in glob.glob(self.load_h5ad[0]+"/*unprocessed.h5ad"):
                        self.concatenated_cell_dict[items.split("/")[-1]\
                        .split('.')[0].split("_")[0]] = sc.read(items, cache=True)

        else:
            self.marker_genes = marker_genes
            self.cell_path_dict = cell_path_dict
            self.ensembl2symbol = ensembl2symbol
            self.cell_batch_names = [batch for batch in self.cell_path_dict.keys()]
            self.cell_dict = {cell_name: [] for cell_name in cell_path_dict.keys()}
            self.marker_gene_dict = {cell_name: {gene: {'sum': None, 'full_set': None, 'total_cells': None} for gene in self.marker_genes} for cell_name in cell_path_dict.keys()}

            print("1. Loading single cell data.")
            for keys, values in self.cell_path_dict.items():
                adata = sc.read(values['filename_data'][0]).transpose()

                adata.var.index = ["-".join(str(hugo).replace("'","-").\
                split("-")[1:-1]) \
                for hugo in np.loadtxt(values['filename_genes'][0], dtype='S')[:, 1]]

                adata.var = adata.var.rename(index=self.ensembl2symbol)
                adata.obs_names = [keys.split('_')[-1] + "_" + "-".\
                join(str(seq).replace("'","-").split("-")[1:-1])\
                 for seq in np.loadtxt(values['filename_barcodes'][0], dtype='S')]

                adata.obs['batch_name'] = keys

                print("Batch: ", keys)
                print(adata)
                print("Total cells: ", adata.shape[0])

                try:
                    if self.output_dir:
                        if not os.path.exists(self.output_dir[0]+"/"+keys):
                            os.makedirs(self.output_dir[0]+"/"+keys+"/preprocessing_figures")
                        if not os.path.exists(self.output_dir[0]+"/preprocessing_summary"):
                            os.makedirs(self.output_dir[0]+"/preprocessing_summary")
                    else:
                        if not os.path.exists(keys):
                            os.makedirs(keys+"/preprocessing_figures")
                except OSError:
                    print ('Error: Creating directory.')

                self.cell_dict[keys].append({keys:adata})
                gene_temp_list = []
                for gene in adata.var_names.tolist():
                    if gene in self.marker_genes:
                        print(gene + ": " + str(sum(adata[:,[gene]].X)))
                        self.marker_gene_dict[keys][gene]['sum'] = sum(adata[:,[gene]].X)
                        self.marker_gene_dict[keys][gene]['full_set'] = adata[:, [gene]].X
                        self.marker_gene_dict[keys][gene]['total_cells'] = adata.shape[0]
                        gene_temp_list.append(gene)
                self.missing_marker_gene(gene_temp_list)
                print("\n")

            # Concatenate each cell batch data set
            self.concatenated_cell_dict = {key: [] for key in self.cell_dict.keys()}
            # for keys, values in self.cell_dict.items():
            #     cell_dict[keys.split('_')[0]].append(values[0])

            print("2. Concatenating single cell data. ")
            for keys, values in self.cell_dict.items():
                    print("Batch:", keys)

                    self.concatenated_cell_dict[keys] = values[0][list(values[0].\
                    keys())[-1]].concatenate([list(items.values())[0]\
                     for items in values[1:]], join='outer')

                    print(self.concatenated_cell_dict[keys])
                    self.check_marker_gene(self.concatenated_cell_dict[keys])
                    print("\n")

            self.output_unprocessed_h5ad = output_unprocessed_h5ad
            if self.output_unprocessed_h5ad:
                for output_name in self.concatenated_cell_dict.keys():
                            if self.output_dir:
                                self.concatenated_cell_dict[output_name].\
                                write(self.output_dir[0]  +\
                                 "/" + output_name+ "/gene_matrices/" + output_name +\
                                  "_unprocessed.h5ad")
                            else:
                                self.concatenated_cell_dict[output_name].write(output_name+\
                                 "/gene_matrices/" + output_name + "_unprocessed.h5ad")

            # Update
            self.only_export_unprocessed_h5ad = only_export_unprocessed_h5ad
            if self.only_export_unprocessed_h5ad:
                for output_name in self.concatenated_cell_dict.keys():
                    if self.output_dir:
                        self.concatenated_cell_dict[output_name]\
                        .write(self.output_dir[0]  + output_name+"_unprocessed.h5ad")
                    else:
                        self.concatenated_cell_dict[output_name].write(output_name+ "_unprocessed.h5ad")
                exit()

        # Create a mapping dictionary
        self.batch_mapping_dict = {}
        inc = 0
        for keys in self.concatenated_cell_dict.keys():
            for k in self.concatenated_cell_dict[keys].obs['batch'].unique():
                self.batch_mapping_dict[k] = inc
                inc += 1

        for keys in self.concatenated_cell_dict.keys():
            self.concatenated_cell_dict[keys].obs["batch"]\
            .replace(self.batch_mapping_dict, inplace=True)

        if minimum_cells:
            self.minimum_cells = minimum_cells
        else:
            self.minimum_cells = 3

        if minimum_genes:
            self.minimum_genes = minimum_genes
        else:
            self.minimum_genes = 200

        self.cell_counts = counts_per_cell_after

        self.thrsh_mito = thrsh_mito
        self.up_thrsh_genes = up_thrsh_genes
        self.low_thrsh_genes = low_thrsh_genes

        self.output_summary_json_dict = {key: {"minimum_cells": self.minimum_cells,\
         "minimum_genes": self.minimum_genes, "counts_per_cell_after": None,\
          "up_thrsh_genes": None, "low_thrsh_genes": None, "thrsh_mito": None,\
           "number_of_variable_genes_found": None}\
            for key in self.concatenated_cell_dict.keys()}

        self.load_h5ad = load_h5ad

    def check_marker_gene(self, input_dict):
        gene_temp_list = []
        for gene in input_dict.var_names.tolist():
            if gene in self.marker_genes:
                print(gene + ": " + str(sum(input_dict[:,[gene]].X)))
                gene_temp_list.append(gene)
        self.missing_marker_gene(gene_temp_list)


    def missing_marker_gene(self, gene_temp_list):
        for gene in self.marker_genes:
            if gene not in gene_temp_list:
                print("Missing: ", gene)

    def filter_cells_and_genes(self):
        cell_and_gene_filtered_dict = copy.deepcopy(self.concatenated_cell_dict)
        print("3. Starting cell and gene filter.")
        for key in cell_and_gene_filtered_dict.keys():
            print("Batch: ", key)
            print(" - Minimum cell count: " + str(self.minimum_cells))
            print(" - Minimum gene count: " + str(self.minimum_genes))

            sc.pp.filter_cells(cell_and_gene_filtered_dict[key],\
             min_genes = self.minimum_cells)

            sc.pp.filter_genes(cell_and_gene_filtered_dict[key], \
            min_cells = self.minimum_genes)

            print(cell_and_gene_filtered_dict[key])
            self.check_marker_gene(cell_and_gene_filtered_dict[key])
            print("\n")
        return cell_and_gene_filtered_dict

    def mitochondria_statistics(self, cell_and_gene_filtered_dict):
        mito_stats_dict = copy.deepcopy(cell_and_gene_filtered_dict)
        print("4. Adding mitochondria percentage and gene count metadata to AnnData structure.")
        for key in mito_stats_dict.keys():
            print("Batch: ", key)
            mito_genes = [name for name in mito_stats_dict[key].var_names if name.startswith('MT.') or name.startswith('MT-')]
            mito_stats_dict[key].obs['percent_mito'] = np.sum(mito_stats_dict[key][:, mito_genes].X, axis=1) / np.sum(mito_stats_dict[key].X, axis=1)
            mito_stats_dict[key].obs['n_counts'] = np.sum(mito_stats_dict[key].X, axis=1)
            print(mito_stats_dict[key])
            if self.output_dir:
                sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'percent_mito', title=key)
                plt.savefig(self.output_dir[0] + "/"+key +"/preprocessing_figures/" +key+"_percent_mito_vs_n_counts")
                plt.close('all')

                sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'n_genes', title=key)
                plt.savefig(self.output_dir[0]+"/"+key +"/preprocessing_figures/" +key+"_n_genes_vs_n_count")
                plt.close('all')
            else:
                sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'percent_mito', title=key)
                plt.savefig(key + "/preprocessing_figures/" +key+"_percent_mito_vs_n_counts")
                plt.close('all')

                sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'n_genes', title=key)
                plt.savefig(key + "/preprocessing_figures/" +key+"_n_genes_vs_n_count")
                plt.close('all')
            print("\n")
        return mito_stats_dict

    def compute_upper_lower_gene_thresholds(self, title, x, y):
        n_counts_vs_n_genes_pd = pd.DataFrame({'x':x, 'y':y})
        n_counts_vs_n_genes_pd.plot(x='x', y='y', kind='hist')
        if self.output_dir:
            plt.savefig(self.output_dir[0]  + '/' + title  +'/preprocessing_figures/' +title+"_n_counts_vs_genes_hist.pdf")
            plt.close('all')
        else:
            plt.savefig(title + "/preprocessing_figures/" + title + "_n_counts_vs_genes_hist.pdf")
            plt.close('all')
        count, bins = np.histogram(y)
        up_thrsh_genes = [(x,y) for x,y in zip(count,bins) if x>10][-1][1]
        low_thrsh_genes = [(x,y) for x,y in zip(count,bins)][0][1]
        return up_thrsh_genes, low_thrsh_genes

    def compute_mitochondria_threshold(self, title, x_, y_):
        n_counts_vs_mito_pct_pd = pd.DataFrame({'x': x_, 'y': y_})
        n_counts_vs_mito_pct_pd.plot(x="x", y="y", kind='hist')
        if self.output_dir:
            plt.savefig(self.output_dir[0] +"/"+title+ "/preprocessing_figures/" + title + '_n_counts_vs_mito_pct.pdf')
            plt.close('all')
        else:
            plt.savefig(title+ "/preprocessing_figures/" + title + '_n_counts_vs_mito_pct.pdf')
            plt.close('all')
        count, bins = np.histogram(y_)
        thrsh_mito = [(x,y) for x,y in zip(count, bins) if x<60 and x>10][-1][-1]
        return thrsh_mito

    def mitochondria_filtering(self, mito_stats_dict):
        mito_filtered_cell_dict = copy.deepcopy(mito_stats_dict)

        print("5. Starting filter for upper gene threshold, lower gene threshold, and mitochondria threshold. ")
        for key in mito_filtered_cell_dict.keys():
            print("Batch: ", key)
            up_thrsh_genes, low_thrsh_genes = self.compute_upper_lower_gene_thresholds(key, np.array(mito_filtered_cell_dict[key].obs['n_counts']), np.array(mito_filtered_cell_dict[key].obs['n_genes']))
            thrsh_mito = self.compute_mitochondria_threshold(key, np.array(mito_filtered_cell_dict[key].obs['n_counts']), np.array(mito_filtered_cell_dict[key].obs['percent_mito']))


            if self.up_thrsh_genes:
                self.output_summary_json_dict[key]['up_thrsh_genes'] = self.up_thrsh_genes
                print(" - Upper gene threshold for " + key + " was set to: " + str(self.up_thrsh_genes))
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']< self.up_thrsh_genes[0], :]
            else:
                print(" - Upper gene threshold for " + key + " was set to: " + str(up_thrsh_genes))
                self.output_summary_json_dict[key]['up_thrsh_genes'] = up_thrsh_genes
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']< up_thrsh_genes, :]

            if self.low_thrsh_genes:
                print(" - Lower gene threshold for " + key + " was set to: " + str(self.low_thrsh_genes))
                self.output_summary_json_dict[key]['low_thrsh_genes'] = self.low_thrsh_genes
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']> self.low_thrsh_genes[0], :]
            else:
                print(" - Lower gene threshold for " + key + " was set to: " + str(low_thrsh_genes))
                self.output_summary_json_dict[key]['low_thrsh_genes'] = low_thrsh_genes
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']> low_thrsh_genes, :]

            if self.thrsh_mito:
                print(" - Mitochondria gene cutoff percentage for " + key + " was set to: " + str(self.thrsh_mito))
                self.output_summary_json_dict[key]['thrsh_mito'] = self.thrsh_mito
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['percent_mito']< self.thrsh_mito[0], :]
            else:
                print(" - Mitochondria gene cutoff percentage for " + key + " was set to: " + str(thrsh_mito))
                self.output_summary_json_dict[key]['thrsh_mito'] = thrsh_mito
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['percent_mito']< thrsh_mito, :]
            self.check_marker_gene(mito_filtered_cell_dict[key])
            print("\n")
        return mito_filtered_cell_dict

    def cell_normalization(self, mito_filtered_cell_dict):
        normalized_mito_filtered_cell_dict = copy.deepcopy(mito_filtered_cell_dict)
        print("6. Normalize single cell data.")
        for key in normalized_mito_filtered_cell_dict.keys():
            print("Batch: ", key)
            if self.cell_counts:
                self.output_summary_json_dict[key]['counts_per_cell_after'] = self.cell_counts[0]
                sc.pp.normalize_per_cell(normalized_mito_filtered_cell_dict[key], counts_per_cell_after=self.cell_counts[0])
            else:
                sc.pp.normalize_per_cell(normalized_mito_filtered_cell_dict[key])
            print(normalized_mito_filtered_cell_dict[key])
            self.check_marker_gene(normalized_mito_filtered_cell_dict[key])
            print("\n")
        return normalized_mito_filtered_cell_dict

    def variable_gene_filtering(self, normalized_mito_filtered_cell_dict):
        variable_gene_filtered_cell_dict = copy.deepcopy(normalized_mito_filtered_cell_dict)
        gene_dispersion_dict = {}
        print("7. Searching for highly variable genes in single cell data.")
        for key in variable_gene_filtered_cell_dict.keys():
            print("Batch: " + key)
            gene_dispersion_dict[key] = sc.pp.filter_genes_dispersion(variable_gene_filtered_cell_dict[key].X, min_mean=0.0125, max_mean=3, min_disp=0.5)
            variable_gene_filtered_cell_dict[key] = variable_gene_filtered_cell_dict[key][:, gene_dispersion_dict[key].gene_subset]
            self.output_summary_json_dict[key]['number_of_varible_genes_found'] = sum(gene_dispersion_dict[key].gene_subset)
            print('Number of variable genes identified in ' + key + ': ', sum(gene_dispersion_dict[key].gene_subset))
            print(variable_gene_filtered_cell_dict[key])
            if self.output_dir:
                sc.pl.filter_genes_dispersion(gene_dispersion_dict[key])
                plt.savefig(self.output_dir[0] + "/" + key + "/preprocessing_figures/" + key +"_gene_dispersion_vs mean_expression")
                plt.close('all')
            else:
                sc.pl.filter_genes_dispersion(gene_dispersion_dict[key])
                plt.savefig(key + "/preprocessing_figures/" + key +"_gene_dispersion_vs mean_expression")
                plt.close('all')
            self.check_marker_gene(variable_gene_filtered_cell_dict[key])
            print("\n")
        return variable_gene_filtered_cell_dict, gene_dispersion_dict

    def log_transformation(self, variable_gene_filtered_cell_dict):
        log_filtered_cell_dict = copy.deepcopy(variable_gene_filtered_cell_dict)
        print("8. Log transform single cell data.")
        for key in log_filtered_cell_dict.keys():
            print("Batch: ", key)
            sc.pp.log1p(log_filtered_cell_dict[key].X)
            print(log_filtered_cell_dict[key])
            self.check_marker_gene(log_filtered_cell_dict[key])
            print("\n")
        return log_filtered_cell_dict

    def regress_out(self, log_filtered_cell_dict):
        np.warnings.filterwarnings("ignore")
        from pandas.core import datetools

        regress_out_filtered_cell_dict = copy.deepcopy(log_filtered_cell_dict)
        print("9. Regress out variables.")
        for key in regress_out_filtered_cell_dict.keys():
            print("Batch: ", key)
            sc.pp.regress_out(regress_out_filtered_cell_dict[key], ['n_counts', 'percent_mito'])
            print(regress_out_filtered_cell_dict[key])
            self.check_marker_gene(regress_out_filtered_cell_dict[key])
            print("\n")
        np.warnings.resetwarnings()
        return regress_out_filtered_cell_dict

    def output_h5_file(self, output_dict, suffix):
        print("Exporting h5ad file.")
        for output_name in output_dict.keys():
            print("Batch: ", output_name)
            if not suffix:
                if self.output_dir:
                    output_dict[output_name].write(self.output_dir[0] + "/" + output_name + "/gene_matrices/" + output_name+".h5ad")
                else:
                    output_dict[output_name].write("/" + output_name + "/gene_matrices/" + output_name+".h5ad")
            else:
                if self.output_dir:
                        output_dict[output_name].write(self.output_dir[0] + "/" + output_name + "/gene_matrices/" + output_name+ suffix + ".h5ad")
                else:
                    output_dict[output_name].write("/gene_matrices/" + output_name + suffix + ".h5ad")
        print("\n")

    def output_summary_json(self):
        print("Exporting summary json dictionary.")
        print(self.output_summary_json_dict)
        for keys in self.marker_gene_dict.keys():
            self.output_summary_json_dict[keys].update(self.marker_gene_dict[keys])

        if self.output_dir:
            with open(self.output_dir[0] + '/preprocessing_summary/in_vitro_HoC_output_summary.pkl', 'wb') as outfile:
                pickle.dump(self.output_summary_json_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open('/preprocessing_summary/in_vitro_HoC_output_summary.pkl', 'w') as outfile:
                pickle.dump(self.output_summary_json_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    def handler(self):
        cell_and_gene_filtered_dict = self.filter_cells_and_genes()
        mito_stats_dict = self.mitochondria_statistics(cell_and_gene_filtered_dict)
        mito_filtered_cell_dict = self.mitochondria_filtering(mito_stats_dict)
        normalized_mito_filtered_cell_dict = self.cell_normalization(mito_filtered_cell_dict)

        # Full preprocessing
        variable_gene_filtered_cell_dict, gene_dispersion_dict = self.variable_gene_filtering(normalized_mito_filtered_cell_dict)
        log_filtered_cell_dict = self.log_transformation(variable_gene_filtered_cell_dict)
        regress_out_filtered_cell_dict = self.regress_out(log_filtered_cell_dict)
        self.output_h5_file(regress_out_filtered_cell_dict, False)
        self.output_summary_json()
        print("\n")

        # Skip variable gene filtering
        print("Rerunning pipeline at step 7 but skipping variable gene filtering.")
        log_filtered_cell_dict_v2 = self.log_transformation(normalized_mito_filtered_cell_dict)
        regress_out_filtered_cell_dict_v2 = self.regress_out(log_filtered_cell_dict_v2)
        self.output_h5_file(regress_out_filtered_cell_dict_v2, '_log_regressed')
        print("\n")

        # Skip variable gene filtering and regression
        print("Rerunning pipeline at step 7 but skipping variable gene filtering and regression.")
        log_filtered_cell_dict_v3 = self.log_transformation(normalized_mito_filtered_cell_dict)
        self.output_h5_file(log_filtered_cell_dict_v3, '_log')


def generate_gene_mapping_dictionary(args):
    print("0. Building ensembl to hugo gene mapping dictionary.")
    print("\n")
    valid_path(args.gene_id_conversion_file[0])
    symbol2ensemble = pd.read_csv(args.gene_id_conversion_file[0], sep='\t', index_col=1)\
['ensembl']
    symbol2ensemble = symbol2ensemble[symbol2ensemble.notnull()]
    ensembl2symbol = dict(zip(symbol2ensemble.values,symbol2ensemble.index.values))
    return ensembl2symbol

def generate_cell_batch_dictionary(args):
    # path_prefix = "/*/outs/filtered_gene_bc_matrices"
    path_prefix = "/*"
    if glob.glob(args.clr_out[0]+path_prefix):
        cell_path_dict = {}
        for matrix_path in glob.glob(args.clr_out[0]+path_prefix):
            cell_path_dict[matrix_path.split("/")[-1]] = {'filename_data': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1].split('.')[-1] == 'mtx'], \
            'filename_genes': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1] == 'genes.tsv'], \
            'filename_barcodes': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1] == 'barcodes.tsv']}
        return cell_path_dict
    else:
        print(INVALID_DIRECTORY_MSG%(args.clr_out[0]))


# error messagesEF09W6D
INVALID_PATH_MSG = "Error: Invalid file path/name. Path %s does not exist."
INVALID_DIRECTORY_MSG = "Error: Invalid prefix path/name. Path %s does not contain cellRanger output directories."

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
    parser.add_argument("--clr_out", type = str, nargs = 1,
                        help = "Path to the cellRanger output directory/directories.")

    parser.add_argument("--gene_id_conversion_file", type = str, nargs = 1,
                        help = "Path to the Ensemble to Hugo gene ID conversion file.")

    parser.add_argument("--only_output_unprocessed_h5ad", type = str, nargs=1, help="Output a h5ad file of the unprocessed data.")

    parser.add_argument("--load_h5ad", type = str, nargs = 1,
                        help = "Path to a filtered gene matrix h5 file.")

    parser.add_argument("--output_unprocessed_h5ad", type = str, nargs = 1, help = "Only create a annData h5 file from the filtered matrices output directory.")

    parser.add_argument("--output_dir", type = str, nargs = 1,
                        help = "Output directory for processed single cell data in h5 file format and process summary json file.")

    parser.add_argument("--min_cells", type = float, nargs = 1, default = 3,
                        help = "Mininmum number of cells. Default is set at 3 cells.")

    parser.add_argument("--min_genes", type = float, nargs = 1, default = 200,
                        help = "Minimum number of genes expressed. Default is set at 200 genes.")

    parser.add_argument("--counts_per_cell_after", type = float, nargs = 1, default = None,
                        help = "Normalize each cell by total counts over all genes, so that every cell has the same total count after normalization. If None, after normalization, each cell has total count equal to the median of the counts_per_cell before normalization. Recommended to be set at 1e4.")

    parser.add_argument("--thrsh_mito", type = float, nargs = 1, default = None,
                        help = "Percentage cutoff of mitochondria genes in a single cell.")

    parser.add_argument("--up_thrsh_genes", type = float, nargs = 1, default = None,
                        help = "Upper limit cutoff for number of genes in a single cell.")

    parser.add_argument("--low_thrsh_genes", type = float, nargs = 1, default = None,
                        help = "Lower limit cutoff for number of genes in a single cell.")

    parser.add_argument("--marker_genes", default=[], nargs='*')

    # parse the arguments from standard input
    args = parser.parse_args()
    if args.clr_out != None and args.gene_id_conversion_file != None:
        print('Running Scanpy version', sc.__version__)
        sc.logging.print_memory_usage()
        print("\n")

        if args.marker_genes:
            print("Marker gene input: ", args.marker_genes)
            print("\n")

        execute = Single_Cell_Data_Wrangling(generate_cell_batch_dictionary(args), generate_gene_mapping_dictionary(args), \
                                            args.min_cells, args.min_genes, args.counts_per_cell_after,
                                            args.thrsh_mito, args.up_thrsh_genes,\
                                            args.low_thrsh_genes, args.output_dir,\
                                             args.load_h5ad, args.output_unprocessed_h5ad,\
                                              args.only_output_unprocessed_h5ad, args.marker_genes)
        execute.handler()

    else:
        print("You are required to provide a path to the matrix file(s) and gene conver id file.")


if __name__ == '__main__':
    main()
