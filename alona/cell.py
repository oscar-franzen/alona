"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import sys
import os
import re
import logging
import pandas as pd
import numpy as np
from joblib import dump, load

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.pyplot import figure
from sklearn.covariance import MinCovDet

from .log import (log_info, log_debug, log_error)
from .constants import (OUTPUT, ORANGE, GENOME)
from .utils import (get_alona_dir, get_time)

from .alonabase import AlonaBase
#from .clustering import AlonaClustering

class AlonaCell(AlonaBase):
    """
    Cell data class.
    """

    def __init__(self):
        self.data = None
        self.data_norm = None
        self.data_ERCC = None
        self.low_quality_cells = None
        self.rRNA_genes = None
        self.pred = None
        super().__init__()

        # make matplotlib more quiet
        logging.getLogger('matplotlib').setLevel(logging.WARNING)

    def validate_counts(self):
        """ Basic data validation of gene expression values. """
        log_debug('Running validate_counts()')
        data_format = self.params['dataformat']

        if (np.any(self.data.dtypes != 'int64') and data_format == 'raw'):
            log_error(msg='Non-count values detected in data matrix while data format is \
set to raw read counts.')
        elif (np.any(self.data.dtypes == 'int64') and data_format == 'log2'):
            log_error(msg='Count values detected in data matrix while data format is \
set to log2.')
        elif self.params['dataformat'] == 'log2':
            if np.any(self.data > 1000):
                log_error(msg='Data do not appear to be log2 transformed.')
        else:
            log_debug('validate_counts() finished without a problem')

        log_debug('Finished validate_counts()')

    def remove_empty(self):
        """ Removes empty cells and genes """
        data_zero = self.data == 0

        cells = data_zero.sum(axis=0)
        genes = data_zero.sum(axis=1)

        total_genes = self.data.shape[0]
        total_cells = self.data.shape[1]

        if np.sum(cells == total_cells) > 0:
            log_info('%s empty cells will be removed' % (np.sum(cells == total_cells)))
            self.data = self.data[self.data.columns[np.logical_not(cells == total_cells)]]
        if np.sum(genes == total_genes) > 0:
            log_info('%s empty genes will be removed' % (np.sum(genes == total_genes)))
            self.data = self.data.loc[self.data.index[np.logical_not(genes == total_genes)]]

        log_info('Current dimensions: %s' % str(self.data.shape))

    def remove_mito(self):
        """ Remove mitochondrial genes. """
        if self.params['nomito']:
            log_debug('Removing mitochondrial genes')

            mt_count = self.data.index.str.contains('^mt-', regex=True, case=False)

            if np.sum(mt_count) > 0:
                log_info('detected and removed %s mitochondrial gene(s)' % np.sum(mt_count))
                self.data = self.data.loc[self.data.index[np.logical_not(mt_count)]]

    def read_counts_per_cell_barplot(self):
        """ Generates a bar plot of read counts per cell. """
        if self.params['dataformat'] == 'raw':
            min_reads = self.params['minreads']
            cell_counts = self.data.sum(axis=0)

            plt.clf()

            figure(num=None, figsize=(5, 5))

            no_cells_remove = np.sum(np.array(sorted(cell_counts, reverse=True)) < \
                 min_reads)
            colors = [ORANGE]*(len(cell_counts)-no_cells_remove) + \
                     ['#999999']*no_cells_remove

            plt.bar(np.arange(len(cell_counts)), sorted(cell_counts, reverse=True),
                    color=colors)
            plt.ylabel('raw read counts')
            plt.xlabel('cells (sorted on highest to lowest)')
            plt.axhline(linewidth=1, color='r', y=min_reads)
            plt.legend(handles=[
                mpatches.Patch(color='#E69F00', label='Passed'),
                mpatches.Patch(color='#999999', label='Removed')
            ])
            plt.title('sequencing reads')
            
            if self.params['timestamp']:
                plt.figtext(0.05, 0, get_time(), size=5)
            
            plt.savefig(self.get_wd() + \
                OUTPUT['FILENAME_BARPLOT_RAW_READ_COUNTS'], bbox_inches='tight')
            plt.close()

    def genes_expressed_per_cell_barplot(self):
        """ Generates a bar plot of number of expressed genes per cell. """
        genes_expressed = self.data.apply(lambda x: sum(x > 0), axis=0)

        plt.clf()
        figure(num=None, figsize=(5, 5))
        plt.bar(np.arange(len(genes_expressed)), sorted(genes_expressed, reverse=True),
                color=[ORANGE]*len(genes_expressed))
        plt.ylabel('number of genes')
        plt.xlabel('cells (sorted on highest to lowest)')
        plt.title('expressed genes')
        
        if self.params['timestamp']:
            plt.figtext(0.05, 0, get_time(), size=5)

        plt.savefig(self.get_wd() + OUTPUT['FILENAME_BARPLOT_GENES_EXPRESSED'],
                    bbox_inches='tight')
        plt.close()

    def simple_filters(self):
        """
        Removes cells with too few reads (set by `--minreads`).
        Removes "underexpressed" genes (set by `--minexpgenes`). If the value set by
        --minexpgenes is a float, then at least that fraction of cells must express the
        gene; if integer then minimum that number of cells must express the gene.
        """
        if self.params['dataformat'] == 'raw':
            min_reads = self.params['minreads']
            cell_counts = self.data.sum(axis=0)
            log_info('Keeping %s out of %s cells' % (
                np.sum(cell_counts > min_reads), len(cell_counts)))

            self.data = self.data[self.data.columns[cell_counts > min_reads]]

            if self.data.shape[1] < 100:
                log_error(self, msg='After removing cells with < %s reads, \
                   less than 100 reads remain. Please adjust --minreads' % min_reads)
        if self.params['minexpgenes'] > 0:
            log_debug('Filtering genes based on --minexpgenes')
            thres = self.params['minexpgenes']
            
            if thres.is_integer():
                genes_expressed = self.data.apply(lambda x: sum(x > 0), axis=1)
                target_genes = genes_expressed[genes_expressed>thres].index
                d = '{0:,g}'.format(np.sum(genes_expressed <= thres))
                log_info('Removing %s genes.' % d)
                self.data = self.data[self.data.index.isin(target_genes)]
            else:
                genes_expressed = self.data.apply(lambda x: sum(x > 0)/len(x), axis=1)
                d = '{0:,g}'.format(np.sum(genes_expressed <= thres))
                log_info('Removing %s genes.' % d)
                self.data = self.data[genes_expressed > thres]

    def print_dimensions(self):
        """ Prints the new dimensions after quality filtering. """
        log_info('Data dimensions after filtering: %s genes and %s cells' % \
        (self.data_norm.shape[0], self.data_norm.shape[1]))

    def rpkm(self, data):
        """ Normalize expression values as RPKM. """
        """ doi:10.1186/s13059-016-0881-8 """
        log_debug('Normalizing data to RPKM')
        if self.params['species'] == 'human':
            # TODO: Implement RPKM for human. This should be executed _before_ mapping
            # to mouse orthologs (16-May-2019).
            log_info('RPKM for "--species human" is not implemented at the moment.')
            raise NotImplementedError('RPKM for human is not implemented at the moment.')

        # the file contains _combined_ lengths of all exons of the particular gene
        exon_lengths = pd.read_csv(get_alona_dir()+GENOME['MOUSE_EXON_LENGTHS'], delimiter=' ',
                                   header=None)
        exon_lengths.columns = ['gene', 'length']

        # intersects and sorts
        temp = pd.merge(data, exon_lengths, how='inner', left_on=data.index,
                        right_on=exon_lengths['gene'].str.extract('^(.+)\.[0-9]+')[0])

        temp.index = temp['gene']
        exon_lengths = temp['length']
        temp = temp.drop(['gene', 'length'], axis=1)
        temp = temp[temp.columns[1:]]

        # gene length in kilobases
        kb = exon_lengths/1000

        def _foo(x):
            s = sum(x)/1000000
            rpm = x/s
            rpkm = rpm/kb

            return rpkm

        _q = temp.apply(_foo, axis=0) # 0 applying to each column

        data_norm = np.log2(_q+1)
        return data_norm

    @staticmethod
    def _dump(d, fn='foo.joblib'):
        """ Serialize data. """
        log_debug('Writing dump file %s' % fn)
        dump(d, fn)

    def normalize(self, data, fn_out='', input_type='raw', mrnafull=False):
        """ Normalizes gene expression values. """
        log_debug('Inside normalize()')
        
        remove_low_quality = self.params['qc_auto']

        if fn_out != '' and os.path.exists(fn_out):
            log_debug('Loading data matrix from file')
            return load(fn_out)

        data_cp = data.copy()
        
        if remove_low_quality:
            data_cp = data_cp.drop(self.low_quality_cells, axis=1)

        if not mrnafull and input_type == 'raw':
            col_sums = data_cp.apply(lambda x: sum(x), axis=0)
            data_norm = (data_cp / col_sums) * 10000
            data_norm = np.log2(data_norm+1)
        elif mrnafull and input_type == 'raw':
            data_norm = self.rpkm(data_cp)
        elif input_type == 'rpkm':
            log_debug('normalization() Running log2')
            data_norm = np.log2(data_cp+1)
        else:
            data_norm = data_cp
            log_debug('Normalization is not needed.')

        if fn_out != '':
            self._dump(data_norm, fn_out)

        log_debug('Finished normalize()')
        return data_norm

    def lift_ERCC(self):
        """ Moves ERCC (if present) to a separate matrix. """
        log_debug('Inside lift_ERCC()')
        self.data_ERCC = self.data[self.data.index.str.contains('^ERCC[_-]\S+$')]
        self.data = self.data[np.logical_not(self.data.index.str.contains('^ERCC[_-]\S+$'))]
        log_debug('Finishing lift_ERCC()')

    def load_rRNA_genes(self):
        """ Detect rRNA genes. """
        log_debug('Inside load_rRNA_genes()')
        # TODO: Add support for human rRNA genes (12 Jun 2019).
        rRNA_genes = {}
        with open(get_alona_dir() + GENOME['MOUSE_GENOME_ANNOTATIONS'], 'r') as fh:
            for line in fh:
                if re.search('"rRNA"', line):
                    gene = re.search('^gene_id "(.*?)\..*";', line.split('\t')[8]).group(1)
                    rRNA_genes[gene] = 1
        self.rRNA_genes = rRNA_genes
        log_debug('Finished load_rRNA_genes()')

    def find_low_quality_cells(self):
        """
        Finds low quality cells using five metrics:

            1. log-transformed number of molecules detected
            2. the number of genes detected
            3. the percentage of reads mapping to ribosomal
            4. mitochondrial genes
            5. ERCC recovery (if available)

        In detail, this function computes Mahalanobis distances from the quality metrics.
        A robust estimate of covariance is used in the Mahalanobis function.
        Cells with Mahalanobis distances of three standard deviations from the mean are
        considered outliers.

        Input should be raw read counts.
        """

        if self.params['qc_auto'] == 'no':
            return

        log_debug('Inside filter_cells_auto()')

        seed = self.params['seed']
        data = self.data
        rRNA_genes = self.rRNA_genes
        data_ERCC = self.data_ERCC

        reads_per_cell = data.sum(axis=0) # no. reads/cell
        no_genes_det = np.sum(data > 0, axis=0)
        data_rRNA = data.loc[data.index.intersection(rRNA_genes)]
        data_mt = data[data.index.str.contains('^mt-', regex=True, case=False)]

        perc_rRNA = data_rRNA.sum(axis=0)/reads_per_cell*100
        perc_mt = data_mt.sum(axis=0)/reads_per_cell*100
        perc_ERCC = data_ERCC.sum(axis=0)/reads_per_cell*100

        qc_mat = pd.DataFrame({'reads_per_cell' : np.log(reads_per_cell),
                               'no_genes_det' : no_genes_det,
                               'perc_rRNA' : perc_rRNA,
                               'perc_mt' : perc_mt,
                               'perc_ERCC' : perc_ERCC})

        robust_cov = MinCovDet(random_state=seed).fit(qc_mat)
        mahal_dists = robust_cov.mahalanobis(qc_mat)

        MD_mean = np.mean(mahal_dists)
        MD_sd = np.std(mahal_dists)

        thres_lower = MD_mean - MD_sd * 3
        thres_upper = MD_mean + MD_sd * 3

        res = (mahal_dists < thres_lower) | (mahal_dists > thres_upper)

        self.low_quality_cells = data.columns[res].values

        fn = self.get_wd() + OUTPUT['FILENAME_QC_SCORE']
        pd.DataFrame(mahal_dists, index=data.columns).to_csv(fn, header=None)

        log_info('%s low quality cells were removed' % len(self.low_quality_cells))
        log_debug('Finished filter_cells_auto()')
        
    def remove_genes_by_pattern(self):
        """ Remove genes matching regexp specified by exclude_gene. """
        log_debug('Entering remove_genes_by_pattern()')
        re_inp = self.params['exclude_gene']
        if re_inp:
            regexp = re.compile(re_inp)
            r = self.data.index.str.contains(regexp)
            self.data = self.data[np.logical_not(r)]
            log_info('Removed %s genes based on regexp.' % np.sum(r))
        log_debug('Exiting remove_genes_by_pattern()')

    def load_annotations(self):
        """ Load gene annotations from file and make the matrix same length as the data. """
        anno_path = self.params['annotations']
        if not anno_path:
            return
        if not os.path.exists(anno_path):
            return
        anno = pd.read_csv(anno_path,sep='\t', header=None, names=['gene','desc'])
        anno = anno[anno['gene'].isin(self.data_norm.index)]
        
        l = np.logical_not(self.data_norm.index.isin(anno['gene']))
        missing = self.data_norm.index[l]
        r = pd.DataFrame({ 'gene' : missing, 'desc' : ['no annotation']*len(missing) })
        anno = anno.append(r, ignore_index=True)
        anno.index=anno['gene']
        anno = anno.reindex(self.data_norm.index)
        self.anno = anno

    def load_data(self):
        """ Load expression matrix. """
        norm_mat_path = self.get_wd() + '/normdata.joblib'
        if os.path.exists(norm_mat_path):
            self.data_norm = load(norm_mat_path)
            return

        log_debug('loading expression matrix')
        self.data = pd.read_csv(self.get_matrix_file(),
                                delimiter=self._delimiter,
                                header=0 if self._has_header else None)

        if self._has_gene_id_column_id or not self._has_header:
            self.data.index = self.data[self.data.columns[0]]
            self.data = self.data.drop(self.data.columns[0], axis=1)

        # remove duplicate genes (another check)
        if np.any(self.data.index.duplicated(False)):
            log_info('%s duplicated genes still detected, removing them.' % np.sum(
                self.data.index.duplicated(False)))
            self.data = self.data.iloc[np.logical_not(self.data.index.duplicated(False))]

        # initialize data
        self.validate_counts()
        self.remove_empty()
        self.remove_mito()
        self.lift_ERCC()
        self.remove_genes_by_pattern()
        self.read_counts_per_cell_barplot()
        self.simple_filters()
        self.genes_expressed_per_cell_barplot()
        self.load_rRNA_genes()
        self.find_low_quality_cells()

        # normalize gene expression values
        dt = self.params['dataformat']
        mf = self.params['mrnafull']
        wd = self.get_wd()
        self.data_norm = self.normalize(self.data, wd + '/normdata.joblib',
                                            mrnafull = mf, input_type = dt)
        self.data_ERCC = self.normalize(self.data_ERCC, wd + '/normdata_ERCC.joblib',
                                            mrnafull = mf, input_type = dt)
        self.load_annotations()
        self.print_dimensions()

        log_debug('Done loading expression matrix')

    def analysis(self):
        """ Runs the analysis pipeline. """
        log_debug('Running analysis...')
        
        embedding_method = self.params['embedding']
        embedding_path = self.get_wd() + OUTPUT['FILENAME_EMBEDDING_PREFIX'] + \
                             embedding_method + '.csv'
        pca_path = self.get_wd() + OUTPUT['FILENAME_PCA']

        if os.path.exists(embedding_path) and os.path.exists(pca_path):
            log_debug('Loading embeddings from file')
            self.pca_components = pd.read_csv(pca_path, header=None,
                                                                index_col=0)
            self.embeddings = pd.read_csv(embedding_path, header=None,
                                                            index_col=0)
        else:
            self.find_variable_genes()
            self.PCA(pca_path)
            self.embedding(embedding_path)

        self.cluster()
        self.findMarkers()
        
        self.median_exp()
        self.load_markers()
        self.CTA_RANK_F(marker_plot=True)
        self.download_model()
        self.run_model_pred()
        
        self.cell_scatter_plot(title='Colored by cluster')
        self.cell_scatter_plot_w_gene_overlay()
        self.genes_exp_per_cluster(title='Colored by cluster')
        self.violin_top()
