import glob
import json
import argparse

class Check_CIRM_Download(object):

    def __init__(self, path_input, gene_converserion_path):
        self.path_input = path_input
        self.asterick = "/*"

        self.file_check = {'input_path': self.path_input}
        self.gene_converserion_path = gene_converserion_path

    def check_input(self):
        check =  glob.glob(self.path_input+self.asterick)
        if check:
            if len(check) == 1:
                print("Found parent directory: ", check.split("/")[-1])

            else:
                for p_dir in check:
                        print("Parent directory found named: ", p_dir.split("/")[-1])
                        self.file_check[p_dir.split("/")[-1]] = \
                        {'parent_dir_path': p_dir + "/",\
                        'cellRanger_path': None,\
                        'cell_batch_names': [],\

                        'cell_batch_outs_dir': None,\
                        'missing_cell_batch_outs_dir': [],\

                        'filtered_gene_bc_matrices_h5_file_path': None,\
                        'missing_filtered_gene_bc_matrices_h5_file': None,\

                        'filtered_gene_bc_matrices_dir_path': None,\
                        'missing_filtered_gene_bc_matrices_dir': None,\


                        'scanpy_command_line_args': {'recommended':[], \
                        'default':[]}}
                        self.cellRanger_check(p_dir.split("/")[-1])
                        if self.gene_converserion_path:
                            self.generate_scanpy_commands(p_dir.split("/")[-1])
                        print("\n")
                if self.gene_converserion_path:
                    self.export_bash(list(self.file_check.keys()))
        else:
            print("No parent directories were found!")

    def cellRanger_check(self, p_dir):
        if glob.glob(self.file_check[p_dir]['parent_dir_path'] + "/*"):
            print("\t-cellRanger directory found named: ", glob.glob(self.file_check[p_dir]['parent_dir_path'] + "/*")[0].split("/")[-1])
            self.file_check[p_dir]['cellRanger_path'] = glob.glob(self.file_check[p_dir]['parent_dir_path'] + "/*")[0]
            self.cell_batch_check(p_dir)
        else:
            print("\t-Warning " + p_dir + " does not have a output cellRanger directory.")
            pass

    def cell_batch_check(self, p_dir):
        print("\t\t-Batch directories found:")
        for cell_batches in glob.glob(self.file_check[p_dir]['cellRanger_path']+self.asterick):
            print("\t\t\t", cell_batches.split("/")[-1])
            self.file_check[p_dir]['cell_batch_names'].append(cell_batches.split("/")[-1])
        self.outs_check(p_dir)

    def outs_check(self, p_dir):
        print("\n")
        print("\t-Checking that batches have a cellRanger outs directory:")
        outs_dir = [cell_batches.split()[0].split('/')[-2] \
        for cell_batches\
         in glob.glob(self.file_check[p_dir]['cellRanger_path']+self.asterick+"/outs")]
        missing_out = set(outs_dir) ^ set(self.file_check[p_dir]['cell_batch_names'])

        self.file_check[p_dir]['cell_batch_outs_dir'] = [cell_batches for cell_batches\
         in glob.glob(self.file_check[p_dir]['cellRanger_path']+self.asterick+"/outs")]

        if len(missing_out) > 0:
            print("\t\t-Warning these batches are missing the cellRanger outs directory:")
            for items in list(missing_out):
                print("\t\t\t Skipping cellRanger check on batch: ", items)
                self.file_check[p_dir]['missing_cell_batch_outs_dir'].append(items)

        else:
            print("\t\t-All batches were found to have an outs directory")

        self.filtered_gene_bc_matrices_h5_file_check(p_dir)
        self.filtered_gene_bc_matrices_h5_dir_check(p_dir)

    def filtered_gene_bc_matrices_h5_file_check(self, p_dir):
        self.file_check[p_dir]['filtered_gene_bc_matrices_dir_path'] = \
        [item+'/filtered_gene_bc_matrices' for item in self.file_check[p_dir]['cell_batch_outs_dir'] \
        if glob.glob(item+'/filtered_gene_bc_matrices_h5.h5')]

        self.file_check[p_dir]['missing_filtered_gene_bc_matrices_h5_file'] = \
        [item.split("/")[-2] for item in self.file_check[p_dir]['cell_batch_outs_dir'] \
        if not glob.glob(item+'/filtered_gene_bc_matrices_h5.h5')]

        if not self.file_check[p_dir]['missing_filtered_gene_bc_matrices_h5_file']:
            print("\t\t-All batches that were found to have an outs directory have a filtered_gene_bc_matrices_h5 file")
        else:
            print("\t\t-Warning these batches have an outs directory but are missing a filtered_gene_bc_matrices_h5 file:")
            for items in self.file_check[p_dir]['missing_filtered_gene_bc_matrices_h5_file']:
                print(items.split("/")[0])

    def filtered_gene_bc_matrices_h5_dir_check(self, p_dir):
        filtered_gene_bc_matrices_files = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']

        self.file_check[p_dir]['filtered_gene_bc_matrices_h5_file_path'] = \
        [glob.glob(item+"/filtered_gene_bc_matrices/GRCh38/*") for item in self.file_check[p_dir]['cell_batch_outs_dir'] \
        if set([ch.split("/")[-1] for ch in glob.glob(item+"/filtered_gene_bc_matrices/GRCh38/*")])\
        == set(filtered_gene_bc_matrices_files)]

        self.file_check[p_dir]['missing_filtered_gene_bc_matrices_dir'] = \
        [item.split("/")[-2] for item in self.file_check[p_dir]['cell_batch_outs_dir'] \
        if set([ch.split("/")[-1] for ch in glob.glob(item+"/filtered_gene_bc_matrices/GRCh38/*")])\
        != set(filtered_gene_bc_matrices_files)]

        if not self.file_check[p_dir]['missing_filtered_gene_bc_matrices_dir']:
            print("\t\t-All batches that were found to have an outs directory have a filtered_gene_bc_matrices directory \n \t\t with barcodes.tsv, genes.tsv, and matrix.mtx files.")
        else:
            print("\t\t-Warning these batches have an outs directory but are missing a filtered_gene_bc_matrices directory:")
            for items in self.file_check[p_dir]['missing_filtered_gene_bc_matrices_h5_file']:
                print(items.split("/")[0])

    def generate_scanpy_commands(self, p_dir):

        print("\n")
        print("Generating command line arguments for SCANPY script\n")

        self.file_check[p_dir]['scanpy_command_line_args']['recommended'].append('python scanpy_preprocessing.py --clr_out' + self.file_check[p_dir]['cellRanger_path']\
        + '/ --gene_id_conversion_file ' + self.gene_converserion_path[0] + ' --output_unprocessed_h5ad true --counts_per_cell_after 1e4')

        print("Recommended:")
        print('\tpython scanpy_preprocessing.py --clr_out' + self.file_check[p_dir]['cellRanger_path']\
        + '/ --gene_id_conversion_file ' + self.gene_converserion_path[0] + ' --output_unprocessed_h5ad true --counts_per_cell_after 1e4')

        print("\n")
        print("Default:")
        self.file_check[p_dir]['scanpy_command_line_args']['default'].append('python scanpy_preprocessing.py --clr_out' + self.file_check[p_dir]['cellRanger_path']\
        + '/ --gene_id_conversion_file ' + self.gene_converserion_path[0])
        print('\tpython scanpy_preprocessing.py --matrix_file' + self.file_check[p_dir]['cellRanger_path']\
        + '/ --gene_id_conversion_file ' + self.gene_converserion_path[0])

    def export_json(self, p_dir):
        with open('cirm_summary.json', 'w') as outfile:
            json.dump(self.file_check, outfile)

    def export_bash(self, p_dir):
        bash_output = open('cirm_args.sh', 'w')
        bash_output.write('#!/bin/bash')
        bash_output.write('\n')
        for key, value in self.file_check.items():
            if key != 'input_path':
                bash_output.write("\n\n")
                bash_output.write(self.file_check[key]['scanpy_command_line_args']['recommended'][0])
                bash_output.write("\n\n")
                bash_output.write(self.file_check[key]['scanpy_command_line_args']['default'][0])
        bash_output.close()


def main():
    # create a parser object
    parser = argparse.ArgumentParser(description = "Check CIRM download.")
    # defining arguments for parser object
    parser.add_argument("-p", type = str, nargs = 1,
                        help = "Path to the cellRanger output directory/directories. Include the backslash after the raw directory.")

    parser.add_argument("-g", type = str, nargs = 1,
                        help = "Path to the gene coversion file.")

    args = parser.parse_args()
    if args.p:
        if glob.glob(args.p[0]):
            execute = Check_CIRM_Download(args.p[0], args.g)
            execute.check_input()
        else:
            print("The path provided is not valid!")
    else:
        print("Please add a path to the CIRM dataset download or gene conversion file!")

if __name__ == '__main__':
    main()
