import matplotlib
import pylab as plt
matplotlib.pyplot.switch_backend('agg')
from errno import EEXIST
from os import makedirs,path
import glob
import json
import argparse
import numpy as np
import pandas as pd
import louvain
import scanpy.api as sc
import igraph as ig

class Single_Cell_Clustering(object):

    def __init__(self, pca_h5ad, n_pcs, n_neighbors, marker_gene_json):

        self.adata_dict = {}
        for processed_file in glob.glob(pca_h5ad+"/*.h5ad"):
            if processed_file.split('/')[-1].split('_')[-1] == 'pca.h5ad':
                self.adata_dict[processed_file.split('/')[-1].split('_')[0]] = sc.read(processed_file)
        self.n_pcs = n_pcs
        self.n_neighbors = n_neighbors


        if marker_gene_json:
            with open(marker_gene_json[0], 'r') as f:
                self.marker_gene_dict = json.load(f)

    def mkdir_p(self, mypath):
        '''Creates a directory. equivalent to using mkdir -p on the command line'''
        try:
            makedirs(mypath)
        except OSError as exc: # Python >2.5
            if exc.errno == EEXIST and path.isdir(mypath):
                pass
            else: raise

    def compute_tsne(self):
        self.mkdir_p('cluster_Analysis')
        for batches, stc in self.adata_dict.items():
            sc.tl.tsne(self.adata_dict[batches],  random_state = 2, n_pcs = self.n_pcs,\
            use_fast_tsne=True, n_jobs=1)

            if [self.marker_gene_dict]:
                for genes in self.marker_gene_dict[batches]:
                    sc.pl.tsne(self.adata_dict[batches], color=[genes])
                    plt.savefig('cluster_Analysis/'+ batches +'_tSNE_' +genes+'_.png', format="PNG")
                    plt.close()

    def compute_neighborhood_graph_and_umap(self):
        for batches,stc in self.adata_dict.items():
            sc.pp.neighbors(self.adata_dict[batches], n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)
            sc.tl.umap(self.adata_dict[batches])
            if [self.marker_gene_dict]:
                for genes in self.marker_gene_dict[batches]:
                    sc.pl.umap(self.adata_dict[batches], color=[genes])
                    plt.savefig('cluster_Analysis/' + batches + '_umap_'+genes+'_.png', format="PNG")
                    plt.close()

    def louvain_clustering(self):
        for batches, stc in self.adata_dict.items():
            sc.pp.neighbors(self.adata_dict[batches], n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)
            sc.tl.louvain(self.adata_dict[batches])

            sc.pl.tsne(self.adata_dict[batches], color='louvain')
            plt.savefig('cluster_Analysis/{}_tSNE_louvain.png'.format(batches), format="PNG")
            plt.close()

            sc.pl.umap(self.adata_dict[batches], color='louvain')
            plt.savefig('cluster_Analysis/{}_umap_louvain.png'.format(batches), format="PNG")
            plt.close()

    def diffusion_map()

    def handler(self):
        self.compute_tsne()
        self.compute_neighborhood_graph_and_umap()
        self.louvain_clustering()

def main():
    # create a parser object
    parser = argparse.ArgumentParser(description = "Cluster visualization of single cell dataset.")
    # defining arguments for parser object
    parser.add_argument("-pca_h5ad", type = str, nargs = 1,
                        help = "Path to the parent directory containing the cell batches PCA h5ad files.")

    parser.add_argument("-n_pcs", type = float, nargs = 1, default = 20,
                        help = "Number of principal components for tSNE. Default is set at 20.")

    parser.add_argument("-n_neighbors", type = float, nargs = 1, default = 15,
                        help = "The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved.  Default is set at 15.")

    parser.add_argument("-marker_gene_json", type = str, nargs = 1,
                        help = "A JSON dictionary containing the batch names and their marker genes.")

    args = parser.parse_args()
    if args.pca_h5ad:
        if glob.glob(args.pca_h5ad[0]):
            execute = Single_Cell_Clustering(args.pca_h5ad[0],\
             args.n_pcs, args.n_neighbors, args.marker_gene_json)
            execute.handler()
        else:
            print("The path provided is not valid!")
    else:
        print("Please add a path to the single cell pca h5ad file!")

if __name__ == '__main__':
    main()
