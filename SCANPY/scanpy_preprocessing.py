
# importing required modules
import matplotlib
import pylab as plt
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
matplotlib.pyplot.switch_backend('agg')
import os
import glob
import copy
import json
import argparse
import pickle
import numpy as np
import pandas as pd
import openpyxl
import scanpy.api as sc

class Single_Cell_Data_Wrangling(object):
    def __init__(self, clr_out, cell_path_dict, ensembl2symbol,\
                minimum_cells, minimum_genes, counts_per_cell_after, \
                thrsh_mito, up_thrsh_genes, low_thrsh_genes, output_dir,\
                load_h5ad, output_unprocessed_h5ad, only_export_unprocessed_h5ad,\
                marker_genes):

        if marker_genes:
            self.marker_genes = list(map(lambda x: x.strip(), marker_genes.readlines()))
            print(self.marker_genes)
            print("\n")
        self.clr_og = clr_out
        self.clr_out = clr_out[0].split("/")[-3]
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

            self.cell_path_dict = cell_path_dict
            self.ensembl2symbol = ensembl2symbol
            self.cell_batch_names = [batch for batch in self.cell_path_dict.keys()]
            self.cell_dict = {cell_name.split('_')[0]: [] for cell_name in cell_path_dict.keys()}
            # self.marker_gene_dict = {cell_name: {gene: [] for gene in self.marker_genes} for cell_name in cell_path_dict.keys()}
            self.marker_gene_dict = {cell_name: {gene: {'sum': None, 'full_set': None, 'total_cells': None} for gene in self.marker_genes} for cell_name in cell_path_dict.keys()}
            print("1. Loading single cell data.")
            # self.total_cell_count = {}
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
                # self.total_cell_count[keys] = adata.shape[0]

                try:
                    if self.output_dir:
                        if not os.path.exists(self.output_dir[0]+"/"+keys.split('_')[0]):
                            os.makedirs(self.output_dir[0]+"/"+keys.split('_')[0]+"/preprocessing_figures")
                        if not os.path.exists(self.output_dir[0]+"/preprocessing_summary"):
                            os.makedirs(self.output_dir[0]+"/preprocessing_summary")
                    else:
                        if not os.path.exists(keys.split('_')[0]):
                            os.makedirs(keys.split('_')[0]+"/preprocessing_figures")
                except OSError:
                    print ('Error: Creating directory.')


                self.cell_dict[keys.split('_')[0]].append({keys.split('_')[-1]:adata})
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
            self.concatenated_cell_dict = {key: None for key in self.cell_dict.keys()}
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

        self.output_summary = {key: {"minimum_cells": self.minimum_cells,\
         "minimum_genes": self.minimum_genes, "counts_per_cell_after": None,\
          "up_thrsh_genes": None, "low_thrsh_genes": None, "thrsh_mito": None,\
           "number_of_variable_genes_found": None}\
            for key in self.concatenated_cell_dict.keys()}

        self.load_h5ad = load_h5ad

    def generate_metric_summary(self):
        path_prefix = "/*/outs/metrics_summary.csv"
        frame = pd.DataFrame()
        list_ = []
        index_list =[]
        for file_ in  glob.glob(self.clr_og[0] + '/*/outs/metrics_summary.csv'):
            index_list.append(file_.split('/')[-3])
            df = pd.read_csv(file_,index_col=None, header=0)
            list_.append(df)
        frame = pd.concat(list_)
        frame.insert(0, 'batch', index_list)
        final_frame = frame.sort_values('batch')
        # frame.reset_index(drop = True, inplace = True)
        final_frame = final_frame.set_index('batch')
        final_frame = final_frame[['Estimated Number of Cells',\
                'Mean Reads per Cell', 'Median Genes per Cell',\
               'Number of Reads','Sequencing Saturation',\
                'Total Genes Detected', 'Median UMI Counts per Cell']]

        # output_prefix = self.clr_out.split('/')[-3]
        final_frame.to_excel(self.output_dir[0] + 'preprocessing_summary/metrics_summary_' + self.clr_out + '.xlsx')

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
            # if self.output_dir:
            #     sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'percent_mito', title=key)
            #     plt.savefig(self.output_dir[0] + "/"+key +"/preprocessing_figures/" +key+"_percent_mito_vs_n_counts")
            #     sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'n_genes', title=key)
            #     plt.savefig(self.output_dir[0]+"/"+key +"/preprocessing_figures/" +key+"_n_genes_vs_n_count")
            # else:
            #     sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'percent_mito', title=key)
            #     plt.savefig(key + "/preprocessing_figures/" +key+"_percent_mito_vs_n_counts")
            #     sc.pl.scatter(mito_stats_dict[key], x = 'n_counts', y = 'n_genes', title=key)
            #     plt.savefig(key + "/preprocessing_figures/" +key+"_n_genes_vs_n_count")
            print("\n")
        return mito_stats_dict

    def compute_upper_lower_gene_thresholds(self, title, x, y):
        n_counts_vs_n_genes_pd = pd.DataFrame({'x':x, 'y':y})
        n_counts_vs_n_genes_pd.plot(x='x', y='y', kind='hist')
        if self.output_dir:
            plt.savefig(self.output_dir[0]  + '/' + title  +'/preprocessing_figures/' +title+"_n_counts_vs_genes_hist.pdf")
        else:
            plt.savefig(title + "/preprocessing_figures/" + title + "_n_counts_vs_genes_hist.pdf")
        count, bins = np.histogram(y)

        print([(x,y) for x,y in zip(count,bins)])
        up_thrsh_genes = [(x,y) for x,y in zip(count,bins) if x>10][-1][1]
        print([(x,y) for x,y in zip(count,bins) if y<1000])
        low_thrsh_genes = [(x,y) for x,y in zip(count,bins)][1][1]
        return up_thrsh_genes, low_thrsh_genes

    def compute_mitochondria_threshold(self, title, x_, y_):
        n_counts_vs_mito_pct_pd = pd.DataFrame({'x': x_, 'y': y_})
        n_counts_vs_mito_pct_pd.plot(x="x", y="y", kind='hist')
        if self.output_dir:
            plt.savefig(self.output_dir[0] +"/"+title+ "/preprocessing_figures/" + title + '_n_counts_vs_mito_pct.pdf')
        else:
            plt.savefig(title+ "/preprocessing_figures/" + title + '_n_counts_vs_mito_pct.pdf')
        count, bins = np.histogram(y_)
        print("Mito filter!")
        print([(x,y) for x,y in zip(count, bins)])
        print([(x,y) for x,y in zip(count, bins) if x>=50])
        thrsh_mito = [(x,y) for x,y in zip(count, bins) if x>=50][-1][-1]
        print(thrsh_mito)
        return thrsh_mito

    def generate_preprocessing_mito_fig(self, mito_filtered_cell_dict, key):
        # mito threshold boundaries
        mito_col = np.where(mito_filtered_cell_dict[key].obs['percent_mito']<self.output_summary[key]['thrsh_mito'],'grey','r')
        # Scatter plot
        mito_scatter = plt.scatter(mito_filtered_cell_dict[key].obs['n_counts'], mito_filtered_cell_dict[key].obs['percent_mito'],\
        s=1, c=mito_col)
        # Red lines
        plt.hlines(y=self.output_summary[key]['thrsh_mito'], xmin=-1, xmax=50000,   color='red')
        # plt.text(50000,self.output_summary[key]['thrsh_mito']+1, s='Mito Threshold: ' + str(self.output_summary[key]['thrsh_mito']))
        legend_elements = [Line2D([0], [0], color='r', lw=1, label='Mito Threshold: {0:.2f}'.format(round(self.output_summary[key]['thrsh_mito'], 2)))]
        plt.legend(handles=legend_elements, loc='upper right')
        plt.xlabel('n_counts', fontsize=18)
        plt.ylabel('percent_mito', fontsize=16)
        if self.output_dir:
            plt.savefig(self.output_dir[0]+"/"+key +"/preprocessing_figures/" +key+"_percent_mito_vs_n_counts")
            plt.close()
        else:
            plt.savefig(key + "/preprocessing_figures/" +key+"_percent_mito_vs_n_counts")
            plt.close()
        plt.close()

    def generate_preprocessing_gene_threshold_fig(self,mito_filtered_cell_dict, key):
        # Cell threshold boundaries
        col = np.where(mito_filtered_cell_dict[key].obs['n_genes']>self.output_summary[key]['up_thrsh_genes'], 'r',\
                    np.where(mito_filtered_cell_dict[key].obs['n_genes']<self.output_summary[key]['low_thrsh_genes'],'blue','grey'))
        # Scatter plot
        cell_scatter = plt.scatter(mito_filtered_cell_dict[key].obs['n_counts'],\
        mito_filtered_cell_dict[key].obs['n_genes'],\
        s=1, c=col)
        # Red lines
        plt.hlines(y=self.output_summary[key]['up_thrsh_genes'], xmin=-1, xmax=50000,   color='red')
        plt.hlines(y=self.output_summary[key]['low_thrsh_genes'], xmin=-1, xmax=50000,   color='blue')
        # plt.text(50000,self.output_summary[key]['up_thrsh_genes']+1, s='Upper Threshold: ' + str(self.output_summary[key]['up_thrsh_genes']))
        # plt.text(50000,self.output_summary[key]['low_thrsh_genes']+1, s='Lower Threshold: ' + str(self.output_summary[key]['low_thrsh_genes']))
        legend_elements = [Line2D([0], [0], color='r', lw=1, label='Upper Gene Threshold: {0:.2f}'.format(round(self.output_summary[key]['up_thrsh_genes'],2))),\
        Line2D([0], [0], color='blue', lw=1, label='Lower Gene Threshold: {0:.2f}'.format(round(self.output_summary[key]['low_thrsh_genes'],2)))]
        plt.legend(handles=legend_elements, loc='upper right')
        plt.xlabel('n_counts', fontsize=18)
        plt.ylabel('n_genes', fontsize=16)
        if self.output_dir:
            plt.savefig(self.output_dir[0]+"/"+key +"/preprocessing_figures/" +key+"_n_genes_vs_n_count")
            plt.close()
        else:
            plt.savefig(key + "/preprocessing_figures/" +key+"_n_genes_vs_n_count")
            plt.close()
        plt.close()

    def mitochondria_filtering(self, mito_stats_dict):
        mito_filtered_cell_dict = copy.deepcopy(mito_stats_dict)

        print("5. Starting filter for upper gene threshold, lower gene threshold, and mitochondria threshold. ")
        for key in mito_filtered_cell_dict.keys():
            print("Batch: ", key)
            up_thrsh_genes, low_thrsh_genes = self.compute_upper_lower_gene_thresholds(key, np.array(mito_filtered_cell_dict[key].obs['n_counts']), np.array(mito_filtered_cell_dict[key].obs['n_genes']))
            thrsh_mito = self.compute_mitochondria_threshold(key, np.array(mito_filtered_cell_dict[key].obs['n_counts']), np.array(mito_filtered_cell_dict[key].obs['percent_mito']))

            # Upper threshold
            if self.up_thrsh_genes and self.low_thrsh_genes:
                self.output_summary[key]['up_thrsh_genes'] = self.up_thrsh_genes
                self.output_summary[key]['low_thrsh_genes'] = self.low_thrsh_genes
                print(" - Upper gene threshold for " + key + " was set to: " + str(self.up_thrsh_genes))
                print(" - Lower gene threshold for " + key + " was set to: " + str(self.low_thrsh_genes))
                self.generate_preprocessing_gene_threshold_fig(mito_filtered_cell_dict, key)
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']< self.up_thrsh_genes[0], :]
            else:
                self.output_summary[key]['up_thrsh_genes'] = up_thrsh_genes
                self.output_summary[key]['low_thrsh_genes'] = low_thrsh_genes
                print(" - Upper gene threshold for " + key + " was set to: " + str(up_thrsh_genes))
                print(" - Lower gene threshold for " + key + " was set to: " + str(low_thrsh_genes))
                self.generate_preprocessing_gene_threshold_fig(mito_filtered_cell_dict, key)
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']< up_thrsh_genes, :]

            # Lower threshold
            # if self.low_thrsh_genes:
            #     print(" - Lower gene threshold for " + key + " was set to: " + str(self.low_thrsh_genes))
            #     self.output_summary[key]['low_thrsh_genes'] = self.low_thrsh_genes
            #     mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']> self.low_thrsh_genes[0], :]
            # else:
            #     print(" - Lower gene threshold for " + key + " was set to: " + str(low_thrsh_genes))
            #     self.output_summary[key]['low_thrsh_genes'] = low_thrsh_genes
            #     mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['n_genes']> low_thrsh_genes, :]


            # Mitochondria threshold
            if self.thrsh_mito:
                print(" - Mitochondria gene cutoff percentage for " + key + " was set to: " + str(self.thrsh_mito))
                self.output_summary[key]['thrsh_mito'] = self.thrsh_mito
                self.generate_preprocessing_mito_fig(mito_filtered_cell_dict, key)
                mito_filtered_cell_dict[key] = mito_filtered_cell_dict[key][mito_filtered_cell_dict[key].obs['percent_mito']< self.thrsh_mito[0], :]
            else:
                print(" - Mitochondria gene cutoff percentage for " + key + " was set to: " + str(thrsh_mito))
                self.output_summary[key]['thrsh_mito'] = thrsh_mito
                self.generate_preprocessing_mito_fig(mito_filtered_cell_dict, key)
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
                self.output_summary[key]['counts_per_cell_after'] = self.cell_counts[0]
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
            self.output_summary[key]['number_of_varible_genes_found'] = sum(gene_dispersion_dict[key].gene_subset)
            print('Number of variable genes identified in ' + key + ': ', sum(gene_dispersion_dict[key].gene_subset))
            print(variable_gene_filtered_cell_dict[key])
            if self.output_dir:
                sc.pl.filter_genes_dispersion(gene_dispersion_dict[key])
                plt.savefig(self.output_dir[0] + "/" + key + "/preprocessing_figures/" + key +"_gene_dispersion_vs mean_expression")
            else:
                sc.pl.filter_genes_dispersion(gene_dispersion_dict[key])
                plt.savefig(key + "/preprocessing_figures/" + key +"_gene_dispersion_vs mean_expression")
            self.check_marker_gene(variable_gene_filtered_cell_dict[key])
            print("\n")
        return variable_gene_filtered_cell_dict, gene_dispersion_dict

    def log_transformation(self, variable_gene_filtered_cell_dict):
        log_filtered_cell_dict = copy.deepcopy(variable_gene_filtered_cell_dict)
        print("8. Log transform single cell data.")
        for key in log_filtered_cell_dict.keys():
            print("Batch: ", key)
            sc.pp.log1p(log_filtered_cell_dict[key], copy=True).write(self.output_dir[0] + "/" + key + "/gene_matrices/" + key+"_raw.h5ad")
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

    def generate_data_frame(self):
        columns = set()
        index = []
        for batch, data in self.output_summary.items():
            index.append(batch)
            for column_key in data.keys():
                columns.add(column_key)
        return pd.DataFrame(index=index, columns=list(columns))

    def populate_data_frame(self):
        report_df = self.generate_data_frame()
        for batch, data in self.output_summary.items():
            for column_key in data.keys():
                report_df[column_key][batch] = self.output_summary[batch][column_key]
        return report_df

    def export_excel_report(self):
        print("Exporting preprocessing summary report.")
        report_df = self.populate_data_frame()
        report_df.to_excel(self.output_dir[0] + 'preprocessing_summary/' + self.clr_out + "_preprocessing_summary_table.xlsx")

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
        self.export_excel_report()
        self.generate_metric_summary()
        print("\n")

        # Skip variable gene filtering
        print("Rerunning pipeline at step 7 but skipping variable gene filtering.")
        log_filtered_cell_dict_v2 = self.log_transformation(normalized_mito_filtered_cell_dict)
        regress_out_filtered_cell_dict_v2 = self.regress_out(log_filtered_cell_dict_v2)
        self.output_h5_file(regress_out_filtered_cell_dict_v2, '_log_regressed')
        # print("\n")

        # Skip variable gene filtering and regression
        # print("Rerunning pipeline at step 7 but skipping variable gene filtering and regression.")
        # log_filtered_cell_dict_v3 = self.log_transformation(normalized_mito_filtered_cell_dict)
        # self.output_h5_file(log_filtered_cell_dict_v3, '_log')


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
    path_prefix = "/*/outs/filtered_gene_bc_matrices"
    # path_prefix = "/*"
    if glob.glob(args.clr_out[0]+path_prefix):
        cell_path_dict = {}
        for matrix_path in glob.glob(args.clr_out[0]+path_prefix):
            cell_path_dict[matrix_path.split("/")[-3]] = {'filename_data': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1].split('.')[-1] == 'mtx'], \
            'filename_genes': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1] == 'genes.tsv'], \
            'filename_barcodes': [value for value in glob.glob(matrix_path+"/*/*") if value.split('/')[-1] == 'barcodes.tsv']}
        return cell_path_dict
    else:
        print(INVALID_DIRECTORY_MSG%(args.clr_out[0]))




# error messagesEF09W6D
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
    parser.add_argument("--clr_out", type = str, nargs = 1,
                        help = "Path to the cellRanger output directory/directories.")

    parser.add_argument("--gene_id_conversion_file", type = str, nargs = 1,
                        help = "Path to the Ensemble to Hugo gene ID conversion file.")

    parser.add_argument("--only_output_unprocessed_h5ad", type = str, nargs=1, help="Output a h5ad file of the unprocessed data.")

    parser.add_argument("--load_h5ad", type = str, nargs = 1,
                        help = "Path to a filtered gene matrix h5 file.")

    parser.add_argument("--output_unprocessed_h5ad", type = str, nargs = 1, help = "Only create a annData h5 file from the filtered matrices output directory.")

    parser.add_argument("--output_dir", type = str, nargs = 1,
                        help = "Output directory for processed single cell data and figures.")

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

    parser.add_argument('--marker_gene_file', type=argparse.FileType('r'),
                        help = 'Path to the txt file containing marker genes.')

    # parse the arguments from standard input
    args = parser.parse_args()
    if args.clr_out != None and args.gene_id_conversion_file != None:
        print('Running Scanpy version', sc.__version__)
        sc.logging.print_memory_usage()
        print("\n")

        execute = Single_Cell_Data_Wrangling(args.clr_out, generate_cell_batch_dictionary(args), generate_gene_mapping_dictionary(args), \
                                            args.min_cells, args.min_genes, args.counts_per_cell_after,
                                            args.thrsh_mito, args.up_thrsh_genes,\
                                            args.low_thrsh_genes, args.output_dir,\
                                             args.load_h5ad, args.output_unprocessed_h5ad,\
                                             args.only_output_unprocessed_h5ad, args.marker_gene_file)
        execute.handler()

    else:
        print("You are required to provide a path to the matrix file(s) and gene conversion id file.")


if __name__ == '__main__':
    main()
