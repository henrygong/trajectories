
import matplotlib
import pylab as plt
matplotlib.pyplot.switch_backend('agg')

import glob
import argparse
import scanpy.api as sc
from errno import EEXIST
from os import makedirs,path
import os
import subprocess

class Single_Cell_Cluster_Figures(object):

    def __init__(self, cluster_h5ad, marker_gene_file, cell_type):

        self.prefix = cluster_h5ad
        self.cell_type = cell_type

        for clusters in glob.glob(self.prefix + "/*"):
            if clusters.split("/")[-1] != 'preprocessing_summary':
                try:
                    if not os.path.exists(clusters + "/cluster_analysis"):
                        os.makedirs(clusters + "/cluster_analysis")

                    # Create marker gene analysis directory
                    if not os.path.exists(clusters + "/cluster_analysis/marker_gene_analysis"):
                        os.makedirs(clusters + "/cluster_analysis/marker_gene_analysis")

                    if not os.path.exists(clusters + "/cluster_analysis/marker_gene_analysis/tSNE"):
                        os.makedirs(clusters + "/cluster_analysis/marker_gene_analysis/tSNE")

                    if not os.path.exists(clusters + "/cluster_analysis/marker_gene_analysis/louvain"):
                        os.makedirs(clusters + "/cluster_analysis/marker_gene_analysis/louvain")

                    if not os.path.exists(clusters + "/cluster_analysis/marker_gene_analysis/umap"):
                        os.makedirs(clusters + "/cluster_analysis/marker_gene_analysis/umap")

                    if not os.path.exists(clusters + "/cluster_analysis/marker_gene_analysis/heatmap"):
                        os.makedirs(clusters + "/cluster_analysis/marker_gene_analysis/heatmap")

                    if not os.path.exists(clusters + "/cluster_analysis/marker_gene_analysis/dotplot"):
                        os.makedirs(clusters + "/cluster_analysis/marker_gene_analysis/dotplot")

                    if not os.path.exists(clusters + "/cluster_analysis/marker_gene_analysis/cluster_gene_rankings"):
                        os.makedirs(clusters + "/cluster_analysis/marker_gene_analysis/cluster_gene_rankings")

                    # Create cell type analysis directory
                    if not os.path.exists(clusters + "/cluster_analysis/cell_type_analysis"):
                        os.makedirs(clusters + "/cluster_analysis/cell_type_analysis")

                    if not os.path.exists(clusters + "/cluster_analysis/cell_type_analysis/tSNE"):
                        os.makedirs(clusters + "/cluster_analysis/cell_type_analysis/tSNE")

                    if not os.path.exists(clusters + "/cluster_analysis/cell_type_analysis/louvain"):
                        os.makedirs(clusters + "/cluster_analysis/cell_type_analysis/louvain")

                    if not os.path.exists(clusters + "/cluster_analysis/cell_type_analysis/umap"):
                        os.makedirs(clusters + "/cluster_analysis/cell_type_analysis/umap")

                except OSError:
                    print("Error creating directory.")

        # Create dictionary containing cell cluster cluster cluster_matrices
        self.cluster_matrices_dict = {}
        for cluster_matrices in glob.glob(self.prefix + "/*/cluster_matrices/*"):
            self.cluster_matrices_dict[cluster_matrices.split('/')[-3]] = sc.read(cluster_matrices)

        # Create marker gene list with provided txt file
        self.marker_gene_file = marker_gene_file
        if self.marker_gene_file:
            self.marker_gene_list = list(map(lambda x: x.strip(), marker_gene_file.readlines()))



    def cluster_figures(self):
        for batches, stc in self.cluster_matrices_dict.items():
            if self.marker_gene_file:
                intersecting_genes = list(set(self.cluster_matrices_dict[batches].var.index.tolist()).intersection(set(self.marker_gene_list)))

                sc.pl.tsne(self.cluster_matrices_dict[batches], color=intersecting_genes)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/tSNE/' + batches + '_pca_tsne.png', format="PNG")
                plt.close()

                sc.pl.umap(self.cluster_matrices_dict[batches], color=intersecting_genes)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/umap/' + batches + '_umap.png', format="PNG")
                plt.close()

                sc.pl.tsne(self.cluster_matrices_dict[batches], color='louvain', title=batches)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/louvain/{}_tSNE_louvain.png'.format(batches), format="PNG")
                plt.close()

                sc.pl.umap(self.cluster_matrices_dict[batches], color='louvain', title=batches)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/louvain/{}_umap_louvain.png'.format(batches), format="PNG")
                plt.close()

                sc.pl.heatmap(self.cluster_matrices_dict[batches], var_names=intersecting_genes, groupby='louvain', log=True)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/heatmap/{}_louvain_heatmap.png'.format(batches), format="PNG")
                plt.close()

                sc.pl.heatmap(self.cluster_matrices_dict[batches], var_names=intersecting_genes, groupby='batch_name', log=True)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/heatmap/{}_heatmap.png'.format(batches), format="PNG")
                plt.close()

                sc.pl.dotplot(self.cluster_matrices_dict[batches], var_names=intersecting_genes, groupby='louvain', log=True)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/dotplot/{}_louvain_plot.png'.format(batches), format="PNG")
                plt.close()

                sc.pl.dotplot(self.cluster_matrices_dict[batches], var_names=intersecting_genes, groupby='batch_name', log=True)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/dotplot/{}_dotplot.png'.format(batches), format="PNG")
                plt.close()

                sc.tl.rank_genes_groups(self.cluster_matrices_dict[batches], 'louvain')
                sc.pl.rank_genes_groups(self.cluster_matrices_dict[batches], n_genes=30)
                plt.savefig(self.prefix + batches + '/cluster_analysis/marker_gene_analysis/cluster_gene_rankings/{}_gene_rankings.png'.format(batches), format="PNG")
                plt.close()

            if self.cell_type:
                sc.pl.tsne(self.cluster_matrices_dict[batches], color='batch_name', title=batches)
                plt.savefig(self.prefix + batches + '/cluster_analysis/cell_type_analysis/tSNE/' + batches + '_pca_tsne.png', format="PNG")
                plt.close()

                sc.pl.umap(self.cluster_matrices_dict[batches], color='batch_name', title=batches)
                plt.savefig(self.prefix + batches + '/cluster_analysis/cell_type_analysis/umap/' + batches + '_umap.png', format="PNG")
                plt.close()

                sc.pl.tsne(self.cluster_matrices_dict[batches], color='louvain', title=batches)
                plt.savefig(self.prefix + batches + '/cluster_analysis/cell_type_analysis/louvain/{}_tSNE_louvain.png'.format(batches), format="PNG")
                plt.close()

                sc.pl.umap(self.cluster_matrices_dict[batches], color='louvain', title=batches)
                plt.savefig(self.prefix + batches + '/cluster_analysis/cell_type_analysis/louvain/{}_umap_louvain.png'.format(batches), format="PNG")
                plt.close()


def main():
    # create a parser object
    parser = argparse.ArgumentParser(description = "Generate cluster figures.")
    # defining arguments for parser object
    parser.add_argument("-cluster_h5ad", type = str, nargs = 1,
                        help = "Path to the parent directory containing the cell cluster h5ad files.")

    parser.add_argument('-marker_gene_file', type=argparse.FileType('r'),
                        help = 'Path to the txt file containing marker genes.')

    parser.add_argument('-cell_type', action='store_false', default=True)

    args = parser.parse_args()

    if args.cluster_h5ad:
        if glob.glob(args.cluster_h5ad[0]):
            execute = Single_Cell_Cluster_Figures(args.cluster_h5ad[0],\
             args.marker_gene_file, args.cell_type)
            execute.cluster_figures()
        else:
            print("The path provided is not valid!")
    else:
        print("Please add a path to the single cell cluster h5ad file!")

if __name__ == '__main__':
    main()
