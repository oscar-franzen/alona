"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import re
import os
import sys
import subprocess

# don't use `magic`, libmagic is not installed by default in MacOS
#import magic

import pandas as pd
import scipy.io

from .log import (log_info, log_debug, log_error, log_warning)
from .constants import (GENOME, OUTPUT)
from .utils import (get_alona_dir, random_str, is_binary)

class AlonaBase():
    """
    AlonaBase class
    """

    def __init__(self):
        self.params = None

    def set_params(self, params):
        """ Set initial parameters. """
        if params is None:
            params = {}

        self.params = {
        }

        self._is_binary = None
        self._delimiter = None
        self._has_header = None
        self._has_gene_id_column_id = None
        #self.mouse_symbols = {}
        #self.mouse_ensembls = {}
        self.unmappable = []
        self.state = {}
        self.mat_data_file = 'input.mat'
        self.anno = None

        self.params.update(params)

        # Set default options
        if self.get_wd() is None:
            self.params['output_directory'] = 'alona_out_%s' % random_str()
        if self.params['loglevel'] == 'debug':
            log_debug('*** parameters *********************')
            for par in self.params:
                log_debug('%s : %s' % (par, self.params[par]))
            log_debug('************************************')

        # validate some of the input parameters
        if self.params['minreads'] < 0:
            log_error('--minreads must be a positive integer.')
        if self.params['minexpgenes'] < 0:
            log_error('--minexpgenes must be a positive value.')
        if self.params['nn_k'] < 0:
            log_error('--nn_k cannot be negative.')
        if self.params['prune_snn'] < 0 or self.params['prune_snn'] > 1:
            log_error('--prune_snn must have a value within [0,1]')
        if self.params['perplexity'] < 0:
            log_error('Perplexity cannot be less than 0.')
        if self.params['perplexity'] > 100:
            log_warning('Recommended values of --perplexity is 5-50.')
        if self.params['hvg_n'] <= 0:
            log_error('--hvg must be a positive value')

    def get_wd(self):
        """ Retrieves the name of the output directory. """
        return self.params['output_directory']

    def get_matrix_file(self):
        """ Returns the full path to the input data matrix. """
        return '%s/%s' % (self.get_wd(), self.mat_data_file)

    def create_work_dir(self):
        """ Creates a working directory for temporary and output files. """
        try:
            log_debug('creating output directory: %s' %self.get_wd())
            os.mkdir(self.get_wd())
            os.mkdir(self.get_wd() + '/plots')
            os.mkdir(self.get_wd() + '/csvs')
            os.mkdir(self.get_wd() + '/csvs/CTA_RANK_F')
            #os.mkdir(self.get_wd() + '/csvs/SVM')
        except FileExistsError:
            log_info('Output directory already exists (%s), resuming.' %
                     self.get_wd())

    def is_file_empty(self):
        """ Checks if the file is an empty file. """
        if os.stat(self.params['input_filename']).st_size == 0:
            log_error('Input file is empty.')

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return self

    def matrix_market(self):
        """
        Unpacks data in matrix market format (used by NCBI GEO). A sparse format.
        https://math.nist.gov/MatrixMarket/formats.html
        """
        log_debug('entering matrix_market()')
        input_file = self.params['input_filename']
        out_dir = self.get_wd()
        mat_out = self.get_matrix_file()
        
        if re.search('.tar.gz$', input_file):
            out = subprocess.check_output("tar -ztf %s" % (input_file), shell=True)
        else:
            out = subprocess.check_output("tar -tf %s" % (input_file), shell=True)
        
        files_found = 0
        for fn in out.decode('ascii').split('\n'):
            if fn != '' and (fn == 'barcodes.tsv' or 'genes.tsv' or 'matrix.mtx'):
                files_found += 1

        if files_found == 3:
            # unpack using system utils
            if re.search('.tar.gz$', input_file):
                cmd = 'tar -C %s -zxf %s' % (out_dir, input_file)
            else:
                cmd = 'tar -C %s -xf %s' % (out_dir, input_file)
            subprocess.check_output(cmd, shell=True)
            m = scipy.io.mmread('%s/matrix.mtx' % out_dir)
            m = pd.DataFrame(m.todense())
            bc = pd.read_csv('%s/barcodes.tsv' % out_dir, header=None)
            g = pd.read_csv('%s/genes.tsv' % out_dir, header=None, sep='\t')
            if g.shape[1] > 1:
                g = g[1] + '_' + g[0]
            m.index = g
            m.columns = bc[0]
            m.to_csv(mat_out, sep='\t')
        log_debug('exiting matrix_market()')

    def unpack_data(self):
        """ Unpacks compressed data and if no compression is used, symlinks to data. """
        input_file = self.params['input_filename']
        abs_path = os.path.abspath(input_file)
        mat_out = self.get_matrix_file()
        # Is the file binary?
        self._is_binary = is_binary(input_file)
        if os.path.exists(mat_out):
            log_info('Data unpacking has already been done.')
            return
        # .tar file?
        if re.search(r'(.tar.gz|tar)$', input_file):
            self.matrix_market()
            return
        if self._is_binary:
            log_debug('Input file is binary.')
            out = subprocess.check_output("file %s" % (input_file), shell=True)
            out = out.decode('ascii')
            if re.search(' gzip compressed data,', out):
                log_debug('gzip data detected.')
                # confirm integrity
                try:
                    cmd = 'gzip -t %s' % (abs_path)
                    if self.params['timeout']:
                        cmd = 'timeout 60s %s' % cmd
                    out = subprocess.check_output(cmd, shell=True,
                                                  stderr=subprocess.STDOUT)
                    out = out.decode('ascii')
                except subprocess.CalledProcessError as exc:
                    if re.search('unexpected end of file', exc.output.decode('ascii')):
                        log_error('Error: input file is corrupt.')
                # uncompress
                cmd = 'gzip -d -c %s > %s' % (abs_path, mat_out)
                if self.params['timeout']:
                    cmd = 'timeout 60s %s' % cmd
                os.system(cmd)
            elif re.search(' Zip archive data,', out):
                log_debug('zip data detected.')
                # confirm integrity
                try:
                    # don't use Zip -T because it only accepts files ending with .zip
                    cmd = 'gzip -t %s' % abs_path
                    if self.params['timeout']:
                        cmd = 'timeout 60s %s' % cmd
                    out = subprocess.check_output(cmd, shell=True,
                                                  stderr=subprocess.STDOUT)
                    out = out.decode('ascii')
                except subprocess.CalledProcessError as exc:
                    if re.search(' unexpected end of file', exc.output.decode('ascii')):
                        log_error('Error: input file is corrupt.')
                # check that the archive only contains one file
                cmd = 'unzip -v %s | wc -l' % abs_path
                if self.params['timeout']:
                    cmd = 'timeout 60s %s' % cmd
                out = subprocess.check_output(cmd, shell=True)
                out = out.decode('ascii').replace('\n', '')
                no_files = int(out) - 5
                if no_files != 1:
                    log_error('More than one file in input archive.')
                # Files  created  by  zip  can be uncompressed by gzip only if they have
                # a single member compressed with the 'deflation' method.
                cmd = 'zcat %s > %s' % (abs_path, mat_out)
                if self.params['timeout']:
                    cmd = 'timeout 60s %s' % cmd
                os.system(cmd)
            elif re.search(' bzip2 compressed data,', out):
                log_debug('bzip2 data detected.')
                # confirm integrity
                try:
                    cmd = 'bzip2 -t %s' % abs_path
                    if self.params['timeout']:
                        cmd = 'timeout 60s %s' % cmd
                    out = subprocess.check_output(cmd, shell=True,
                                                  stderr=subprocess.STDOUT)
                    out = out.decode('ascii')
                except subprocess.CalledProcessError as exc:
                    if re.search('file ends unexpectedly', exc.output.decode('ascii')):
                        log_error('Input file is corrupt.')
                # uncompress
                cmd = 'bzip2 -d -c %s > %s' % (abs_path, mat_out)
                if self.params['timeout']:
                    cmd = 'timeout 60s %s' % cmd
                os.system(cmd)
            else:
                log_error('Invalid format of the input file.')
        else:
            log_debug('Input file is not binary.')
            # Create a symlink to the data
            cmd = 'ln -sfn %s %s' % (abs_path, mat_out)
            os.system(cmd)

        #mag = magic.Magic(mime=True)
        # don't use `from_file` (it doesn't follow symlinks)
        #f_type = mag.from_buffer(open(mat_out, 'r').read(1024))
        #if f_type != 'text/plain':
        #    log_error('Input file is not plain text (found type=%s).' % f_type)

    #def cleanup(self):
    #    """ Removes temporary files. """
    #    if self.params['cleanup']:
    #        # remove temporary files
    #        log_debug('Cleaning up (`--cleanup` is set)')
    #        
    #        for item in ('input.mat',)+tuple(OUTPUT.values()):
    #            log_debug('removing %s' % garbage)
    #
    #            try:
    #                os.remove('%s/%s' % (self.get_wd(), garbage))
    #            except FileNotFoundError:
    #                log_debug('Not found: %s' % garbage)

    def _guess_delimiter(self):
        """ Determines the data delimiter character. """
        dcount = {' ' : 0, '\t' : 0, ',' : 0}
        i = 0
        fh = open(self.get_matrix_file(), 'r')
        for line in fh:
            dcount[' '] += line.count(' ')
            dcount['\t'] += line.count('\t')
            dcount[','] += line.count(',')
            i += 1
            if i > 10:
                break
        fh.close()
        d_sorted = sorted(dcount, key=dcount.get, reverse=True)
        return d_sorted[0]

    def get_delimiter(self):
        """ Figures out the data delimiter of the input data. """
        used_delim = ''
        if self.params['delimiter'] == 'auto':
            used_delim = self._guess_delimiter()
        else:
            used_delim = self.params['delimiter'].upper()
            if used_delim == 'TAB':
                used_delim = '\t'
            elif used_delim == 'SPACE':
                used_delim = ' '
        log_debug('delimiter is: "%s" (ASCII code=%s)' % (used_delim, ord(used_delim)))
        self._delimiter = used_delim
        return used_delim

    def has_header(self):
        """ Determines if the input data has a header. """
        ret = None
        if self.params['header'] == 'auto':
            with open(self.get_matrix_file(), 'r') as fh:
                first_line = next(fh).replace('\n', '')
                total = 0
                count_digit = 0
                for item in first_line.split(self._delimiter):
                    if item.replace('.', '', 1).isdigit():
                        count_digit += 1
                    total += 1
                ret = not total == count_digit
        else:
            ret = (self.params['header'] == 'yes')
        # if all fields are non-numerical, it's likely a header
        self._has_header = ret
        log_debug('has header: %s' % self._has_header)
        return ret

    def sanity_check_columns(self):
        """
        Sanity check on data integrity. Terminates if column count is not consistent.
        """
        with open(self.get_matrix_file(), 'r') as fh:
            if self._has_header:
                next(fh)
            cols = {}
            for line in fh:
                no_columns = len(line.split(self._delimiter))
                cols[no_columns] = 1
            if len(cols.keys()) > 1:
                log_error('Rows in your data matrix have different \
number of columns (every row must have the same number of columns).')
            log_info('%s columns detected.' % '{:,}'.format(cols.popitem()[0]-1))

    def sanity_check_genes(self):
        """ Sanity check on gene count. Terminates if gene count is too low. """
        with open(self.get_matrix_file(), 'r') as fh:
            if self._has_header:
                next(fh)
            count = 0
            for _ in fh:
                count += 1
            if count < 1000:
                log_error('Number of genes (n=%s) in the input data is too low.' % count)
            log_info('%s genes detected' % '{:,}'.format(count))

    def sanity_check_gene_dups(self):
        """ Checks for gene duplicates. """
        with open(self.get_matrix_file(), 'r') as fh:
            if self._has_header:
                next(fh)
            genes = {}
            for line in fh:
                gene = line.split(self._delimiter)[0]
                if not gene in genes:
                    genes[gene] = 1
                else:
                    genes[gene] += 1
            _cancel = 0
            for gene in genes:
                if genes[gene] > 1:
                    _cancel = 1
                    log_info('%s has duplicates' % gene)
            if _cancel:
                log_error('Gene duplicates detected.')

    #def load_mouse_gene_symbols(self):
    #    """ Loads genome annotations. """
        # with open(get_alona_dir() + GENOME['MOUSE_GENOME_ANNOTATIONS'], 'r') as fh:
            # for line in fh:
                # if not re.search('gene_id "ERCC_', line):
                    # m = re.search(r'gene_id "(.+?)_(.+?)\.[0-9]+"', line)
                    # symbol, ensembl = m.group(1), m.group(2)
                    # self.mouse_symbols[symbol] = '%s_%s' % (symbol, ensembl)
                    # self.mouse_ensembls[ensembl] = '%s_%s' % (symbol, ensembl)
        # with open(get_alona_dir() + GENOME['ENTREZ_GENE_IDS'], 'r') as fh:
            # next(fh)
            # for line in fh:
                # gene_id_as_number, ens = line.rstrip('\n').split(' ')
                # if gene_id_as_number != 'null':
                    # self.mouse_entrez[gene_id_as_number] = ens

    def check_gene_name_column_id_present(self):
        """ Checks if the header line has a column attribute. """
        log_debug('running check_gene_name_column_id_present()')
        with open(self.get_matrix_file(), 'r') as fh:
            header = next(fh)
            line2 = next(fh)
            self._has_gene_id_column_id = len(header.split(self._delimiter)) == \
                                          len(line2.split(self._delimiter))
            if self._has_gene_id_column_id:
                log_info('A column ID for the gene symbols was identified.')

    def prepare(self):
        """ Prepares data for analysis. """
        self.create_work_dir()
        settings_file = self.get_wd() + OUTPUT['FILENAME_SETTINGS']
        if os.path.exists(settings_file):
            with open(settings_file, 'r') as f:
                for line in f:
                    key, val = line.replace('\n', '').split('\t')
                    if key == 'input_filename':
                        a = val.split('/')[-1]
                        basename = self.params['input_filename'].split('/')[-1]
                        if a != basename:
                            log_error('Input filename is different with the filename \
specified in settings.txt (%s). Try using a different output directory' % a)
        with open(settings_file, 'w') as f:
            for p in self.params:
                f.write('%s\t%s\n' % (p, self.params[p]))
        if os.path.exists(self.get_wd() + '/normdata.joblib'):
            log_debug('prepare(): normdata.joblib detected, skipping some steps')
            return
        self.is_file_empty()
        self.unpack_data()
        self.get_delimiter()
        self.has_header()
        self.sanity_check_columns()
        self.sanity_check_genes()
        self.sanity_check_gene_dups()
        #self.load_mouse_gene_symbols()
        self.check_gene_name_column_id_present()
