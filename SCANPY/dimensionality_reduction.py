import matplotlib
import pylab as plt
matplotlib.pyplot.switch_backend('agg')
from errno import EEXIST
from os import makedirs,path
import glob
import argparse
import numpy as np
import pandas as pd
import scanpy.api as sc
import copy
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

class PCA_Analysis(object):

    def __init__(self, ah5_path):

        self.adata_dict = {}
        for processed_file in glob.glob(ah5_path+"/*.h5ad"):
            if len(processed_file.split('/')[-1].split('_')) == 1:

                self.adata_dict[processed_file.split('/')[-1].split('.')[0]] = sc.read(processed_file)

    def mkdir_p(self, mypath):
        '''Creates a directory. equivalent to using mkdir -p on the command line'''
        try:
            makedirs(mypath)
        except OSError as exc: # Python >2.5
            if exc.errno == EEXIST and path.isdir(mypath):
                pass
            else: raise
    def scale_data(self):

        scaled_adata_dict = copy.deepcopy(self.adata_dict)
        for batches, strc in scaled_adata_dict.items():
            scaled_adata_dict[batches] = StandardScaler().fit_transform(strc.X)
        return scaled_adata_dict

    def generate_model(self, scaled_adata_dict):

        for batches, strc in scaled_adata_dict.items():
            pca_model = PCA(n_components=0.99, svd_solver='full')
            X_pca = pca_model.fit_transform(strc)
            self.adata_dict[batches].obsm['X_pca'] = X_pca
            self.adata_dict[batches].varm['PCs'] = pca_model.components_.T
            self.adata_dict[batches].uns['pca'] = {}
            self.adata_dict[batches].uns['pca']['variance'] = pca_model.explained_variance_
            self.adata_dict[batches].uns['pca']['variance_ratio'] =  pca_model.explained_variance_ratio_


    def plot_explained_variance_ratio(self):
        self.mkdir_p('pca_Analysis')
        for batches, strc in self.adata_dict.items():
            plt.plot(np.cumsum(self.adata_dict[batches].uns['pca']['variance_ratio']))
            plt.title(batches)
            plt.xlabel('number of components')
            plt.ylabel('cumulative explained variance')
            plt.savefig('pca_Analysis/explained_variance_ratio_{}.png'.format(batches), format="PNG")
            plt.close()

    def handler(self):
        scaled_adata_dict = self.scale_data()
        self.generate_model(scaled_adata_dict)
        self.plot_explained_variance_ratio()


def main():
    # create a parser object
    parser = argparse.ArgumentParser(description = "Principle component analysis.")
    # defining arguments for parser object
    parser.add_argument("-p_h5ad", type = str, nargs = 1,
                        help = "Path to the parent directory of processed h5ad files.")
    args = parser.parse_args()
    if args.p_h5ad:
        if glob.glob(args.p_h5ad[0]):
            execute = PCA_Analysis(args.p_h5ad[0])
            execute.handler()
        else:
            print("The path provided is not valid!")
    else:
        print("Please add a path to the processed ah5 files!")

if __name__ == '__main__':
    main()
