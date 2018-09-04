import matplotlib
import pylab as plt
matplotlib.pyplot.switch_backend('agg')
from errno import EEXIST
from os import makedirs,path
import os
import glob
import json
import argparse
import numpy as np
import pandas as pd
import louvain
import scanpy.api as sc
import igraph as ig

class Single_Cell_Clustering(object):

    def __init__(self, pca_h5ad, n_pcs, n_neighbors):

        self.adata_dict = {}
        self.prefix = pca_h5ad

        for processed_file in glob.glob(pca_h5ad+"*/principle_component_matrices/*.h5ad"):
            self.adata_dict[processed_file.split('/')[-3]] = sc.read(processed_file)

        for output_dir in glob.glob(pca_h5ad+"*/"):
            if output_dir.split("/")[-1] == 'preprocessing_summary':
                pass
            else:
                try:
                    if not os.path.exists(output_dir + '/cluster_matrices'):
                        os.makedirs(output_dir + '/cluster_matrices')

                    if not os.path.exists(output_dir + "/cluster_analysis"):
                        os.makedirs(output_dir + "/cluster_analysis")

                    # if not os.path.exists(output_dir + "/cluster_analysis/tSNE"):
                    #     os.makedirs(output_dir + "/cluster_analysis/tSNE")
                    #
                    # if not os.path.exists(output_dir + "/cluster_analysis/umap"):
                    #     os.makedirs(output_dir + "/cluster_analysis/umap")
                    #
                    # if not os.path.exists(output_dir + "/cluster_analysis/louvain"):
                    #     os.makedirs(output_dir + "/cluster_analysis/louvain")

                except OSError:
                    print("Error creating directory")


        self.n_pcs = n_pcs
        self.n_neighbors = n_neighbors


    def compute_tsne(self):
        for batches, stc in self.adata_dict.items():
            print(batches)
            sc.tl.tsne(self.adata_dict[batches],  random_state = 2, n_pcs = self.n_pcs,\
            use_fast_tsne=True, n_jobs=1)
            # sc.pl.tsne(self.adata_dict[batches], color='batch_name', title=batches)
            # plt.savefig(self.prefix + batches + '/cluster_analysis/tSNE/' + batches + '_pca_tsne_.png', format="PNG")
            # plt.close()

    def compute_neighborhood_graph_and_umap(self):
        for batches,stc in self.adata_dict.items():
            sc.pp.neighbors(self.adata_dict[batches], n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)
            sc.tl.umap(self.adata_dict[batches])

            # sc.pl.umap(self.adata_dict[batches], color='batch_name', title=batches)
            # plt.savefig(self.prefix + batches + '/cluster_analysis/umap/' + batches + '_umap.png', format="PNG")
            # plt.close()

    def louvain_clustering(self):
        for batches, stc in self.adata_dict.items():
            sc.pp.neighbors(self.adata_dict[batches], n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)
            sc.tl.louvain(self.adata_dict[batches])
            #
            # sc.pl.tsne(self.adata_dict[batches], color='louvain', title=batches)
            # plt.savefig(self.prefix + batches + '/cluster_analysis/louvain/{}_tSNE_louvain.png'.format(batches), format="PNG")
            # plt.close()

            # sc.pl.umap(self.adata_dict[batches], color='louvain', title=batches)
            # plt.savefig(self.prefix + batches + '/cluster_analysis/louvain/{}_umap_louvain.png'.format(batches), format="PNG")
            # plt.close()

    def export_h5ad(self):
        for batches, stc in self.adata_dict.items():
            self.adata_dict[batches].write(self.prefix + "/" + batches + "/cluster_matrices/" + batches +"_clustered.h5ad")


    def handler(self):
        self.compute_tsne()
        self.compute_neighborhood_graph_and_umap()
        self.louvain_clustering()
        self.export_h5ad()

def main():
    # create a parser object
    parser = argparse.ArgumentParser(description = "Cluster single cell data.")
    # defining arguments for parser object
    parser.add_argument("-pca_h5ad", type = str, nargs = 1,
                        help = "Path to the parent directory containing the cell batches PCA h5ad files.")

    parser.add_argument("-n_pcs", type = float, nargs = 1, default = 20,
                        help = "Number of principal components for tSNE. Default is set at 20.")

    parser.add_argument("-n_neighbors", type = float, nargs = 1, default = 15,
                        help = "The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved.  Default is set at 15.")


    args = parser.parse_args()
    if args.pca_h5ad:
        if glob.glob(args.pca_h5ad[0]):
            execute = Single_Cell_Clustering(args.pca_h5ad[0], args.n_pcs, args.n_neighbors)
            execute.handler()
        else:
            print("The path provided is not valid!")
    else:
        print("Please add a path to the single cell pca h5ad file!")

if __name__ == '__main__':
    main()
