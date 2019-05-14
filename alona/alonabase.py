"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Details:
 https://alona.panglaodb.se/

 Contact:
 Oscar Franzen, <p.oscar.franzen@gmail.com>
"""

import re
import os
import uuid
import inspect
import logging
import subprocess
import magic

from .log import (log_info, log_error)
from .exceptions import *

class AlonaBase():
    """
    AlonaBase class
    """

    # pylint: disable=too-many-instance-attributes

    params = None

    def __init__(self, params=None):
        if params is None:
            params = {}

        self.params = {
        }

        self._is_binary = None
        self._delimiter = None
        self._has_header = None
        self._has_gene_id_column_id = None

        self.mouse_symbols = {}
        self.mouse_ensembls = {}
        self.mouse_entrez = {}

        self.unmappable = []

        self.params.update(params)

        # Set default options
        if self.get_working_dir() is None:
            self.params['output_directory'] = 'alona_out_%s' % self.random()
        if self.params['loglevel'] == 'debug':
            logging.debug('*** parameters *********************')
            for par in self.params:
                logging.debug('%s : %s', par, self.params[par])
            logging.debug('************************************')

    def get_working_dir(self):
        return self.params['output_directory']

    def get_matrix_file(self):
        return '%s/input.mat' % self.get_working_dir()

    def random(self):
        """ Get random 8 character string """
        return str(uuid.uuid4()).split('-')[0]

    def create_work_dir(self):
        """ Creates a working directory for temporary and output files. """
        try:
            logging.debug('creating output directory: %s', self.get_working_dir())
            os.mkdir(self.get_working_dir())
        except FileExistsError:
            log_error(self, 'Error: Output directory already exists (%s)' %
                      self.get_working_dir())
            raise

    def is_file_empty(self):
        if os.stat(self.params['input_filename']).st_size == 0:
            raise FileEmptyError()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return self

    def is_binary(self, filename):
        with open(filename, 'rb') as fh:
            for block in fh:
                if b'\0' in block:
                    return True
        return False

    def unpack_data(self):
        """ Unpacks compressed data and if no compression is used, symlinks to data. """
        # Is the file binary?
        self._is_binary = self.is_binary(self.params['input_filename'])
        abs_path = os.path.abspath(self.params['input_filename'])
        mat_out = self.get_matrix_file()

        if self._is_binary:
            logging.debug('Input file is binary.')

            out = subprocess.check_output("file %s" % (self.params['input_filename']),
                                          shell=True)

            out = out.decode('ascii')

            if re.search(' gzip compressed data,', out):
                logging.debug('gzip data detected.')

                # confirm integrity
                try:
                    out = subprocess.check_output('gzip -t %s' % (abs_path), shell=True,
                                                  stderr=subprocess.STDOUT)

                    out = out.decode('ascii')
                except subprocess.CalledProcessError as exc:
                    if re.search('unexpected end of file', exc.output.decode('ascii')):
                        raise FileCorruptError('Error: input file is corrupt.')

                # uncompress
                os.system('gzip -d -c %s > %s' % (abs_path, mat_out))
            elif re.search(' Zip archive data,', out):
                logging.debug('zip data detected.')

                # confirm integrity
                try:
                    # don't use Zip -T because it only accepts files ending with .zip
                    out = subprocess.check_output("gzip -t %s" % (abs_path), shell=True,
                                                  stderr=subprocess.STDOUT)

                    out = out.decode('ascii')
                except subprocess.CalledProcessError as exc:
                    if re.search(' unexpected end of file', exc.output.decode('ascii')):
                        raise FileCorruptError('Error: input file is corrupt.')

                # check that the archive only contains one file
                out = subprocess.check_output('unzip -v %s | wc -l' % (abs_path),
                                              shell=True)

                out = out.decode('ascii').replace('\n', '')
                no_files = int(out) - 5

                if no_files != 1:
                    raise Exception('More than one file in input archive.')

                # Files  created  by  zip  can be uncompressed by gzip only if they have
                # a single member compressed with the 'deflation' method.
                os.system('zcat %s > %s' % (abs_path, mat_out))
            elif re.search(' bzip2 compressed data,', out):
                logging.debug('bzip2 data detected.')

                # confirm integrity
                try:
                    out = subprocess.check_output('bzip2 -t %s' % (abs_path), shell=True,
                                                  stderr=subprocess.STDOUT)

                    out = out.decode('ascii')
                except subprocess.CalledProcessError as exc:
                    if re.search('file ends unexpectedly', exc.output.decode('ascii')):
                        raise FileCorruptError('Input file is corrupt.')

                # uncompress
                os.system('bzip2 -d -c %s > %s' % (abs_path, mat_out))
            else:
                raise InvalidFileFormatError('Invalid format of the input file.')
        else:
            logging.debug('Input file is not binary.')
            # Create a symlink to the data
            cmd = 'ln -sfn %s %s' % (abs_path, mat_out)
            os.system(cmd)

        mag = magic.Magic(mime=True)
        # don't use `from_file` (it doesn't follow symlinks)
        f_type = mag.from_buffer(open(mat_out, 'r').read(1024))

        if f_type != 'text/plain':
            raise InputNotPlainTextError('Input file is not plain text (found type=%s).' %
                                         f_type)

    def cleanup(self):
        """ Removes temporary files. """
        # remove temporary files
        for garbage in ('input.mat',):
            logging.debug('removing %s', garbage)

            try:
                os.remove('%s/%s' % (self.get_working_dir(), garbage))
            except FileNotFoundError:
                logging.debug('Not found: %s', garbage)

    def __guess_delimiter(self):
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
            used_delim = self.__guess_delimiter()
        else:
            used_delim = self.params['delimiter'].upper()
            if used_delim == 'TAB':
                used_delim = '\t'
            elif used_delim == 'SPACE':
                used_delim = ' '

        logging.debug('delimiter is: "%s" (ASCII code=%s)', used_delim, ord(used_delim))
        self._delimiter = used_delim
        return used_delim

    def has_header(self):
        """ Figures out if the input data uses a header or not. """
        ret = None
        if self.params['header'] == 'auto':
            fh = open(self.get_matrix_file(), 'r')
            first_line = next(fh).replace('\n', '')

            total = 0
            count_digit = 0

            for item in first_line.split(self._delimiter):
                if item.replace('.', '', 1).isdigit():
                    count_digit += 1
                total += 1

            fh.close()

            ret = not total == count_digit
        else:
            ret = (self.params['header'] == 'yes')

        # if all fields are non-numerical, it's likely a header
        self._has_header = ret

        logging.debug('has header: %s', self._has_header)

        return ret

    def sanity_check_columns(self):
        """ Sanity check on data integrity. Raises an exception if column count is not
            consistent. """
        fh = open(self.get_matrix_file(), 'r')

        if self._has_header:
            next(fh)

        cols = {}

        for line in fh:
            no_columns = len(line.split(self._delimiter))
            cols[no_columns] = 1

        fh.close()

        if len(cols.keys()) > 1:
            raise IrregularColumnCountError('Rows in your data matrix have different number \
of columns (every row must have the same number of columns).')
        log_info('%s cells detected.' % '{:,}'.format(cols.popitem()[0]))

    def sanity_check_genes(self):
        """ Sanity check on gene count. Raises an exception if gene count is too low. """
        fh = open(self.get_matrix_file(), 'r')
        if self._has_header:
            next(fh)

        count = 0
        for line in fh:
            count += 1

        fh.close()

        if count < 1000:
            raise IrregularGeneCountError('Number of genes in the input data is too low.')
        log_info('%s genes detected' % '{:,}'.format(count))

    def ortholog_mapper(self):
        """ Maps mouse genes to the corresponding human ortholog.
            Only one-to-one orthologs are considered. """
        # human gene symbols to ens
        f = open(os.path.dirname(inspect.getfile(AlonaBase)) +
                 '/genome/hgnc_complete_set.txt', 'r')
        hs_symb_to_hs_ens = {}

        for line in f:
            if re.search(r'^\S+\t\S+\t', line) and re.search('(ENSG[0-9]+)', line):
                hs_symbol = line.split('\t')[1]
                hs_ens = re.search('(ENSG[0-9]+)', line).group(1)
                hs_symb_to_hs_ens[hs_symbol] = hs_ens
        f.close()

        # ortholog mappings
        f = open(os.path.dirname(inspect.getfile(AlonaBase)) +
                 '/genome/human_to_mouse_1_to_1_orthologs.tsv', 'r')
        next(f)

        human_to_mouse = {}
        for line in f:
            if re.search('\tortholog_one2one\t', line):
                foo = line.split('\t')
                human_ens = foo[0]
                mouse_ens = foo[1]
                human_to_mouse[human_ens] = mouse_ens
        f.close()

        f = open(self.get_matrix_file(), 'r')
        ftemp = open(self.get_matrix_file() + '.mapped2mouse.mat', 'w')
        ftemp2 = open(self.get_matrix_file() + '.genes_missing_mouse_orthologs', 'w')

        if self._has_header:
            header = next(f)
            ftemp.write(header)

        orthologs_found = 0

        for line in f:
            # remove quotes
            line = re.sub('"', '', line)
            foo = line.split(self._delimiter)

            gene = foo[0]

            if re.search('.+_ENSG[0-9]+', gene):
                gene = re.search('^.+_(ENSG[0-9]+)', gene).group(1)
            if human_to_mouse.get(gene, '') != '':
                new_gene_name = human_to_mouse[gene]
                ftemp.write('%s%s%s' % (new_gene_name,
                                        self._delimiter,
                                        self._delimiter.join(foo[1:])))
                orthologs_found += 1
            elif hs_symb_to_hs_ens.get(gene, '') != '':
                hs_ens = hs_symb_to_hs_ens[gene]
                mm_ens = human_to_mouse.get(hs_ens, '')
                orthologs_found += 1
                if mm_ens != '':
                    ftemp.write('%s%s%s' % (mm_ens,
                                            self._delimiter,
                                            self._delimiter.join(foo[1:])))
                else:
                    ftemp2.write('%s\n' % (gene))

        f.close()
        ftemp.close()
        ftemp2.close()

        log_info('mapped %s genes to mouse orthologs' % ('{:,}'.format(orthologs_found)))

        return self.get_matrix_file() + '.mapped2mouse.mat'

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
                raise GeneDuplicatesError('Gene duplicates detected.')

    def load_mouse_gene_symbols(self):
        """ Loads genome annotations. """
        fh = open(os.path.dirname(inspect.getfile(AlonaBase)) +
                  '/genome/Mus_musculus.GRCm38.gencode.vM17.primary_assembly.\
annotation.gene_level.ERCC.gtf', 'r')

        for line in fh:
            if not re.search('gene_id "ERCC_', line):
                m = re.search(r'gene_id "(.+?)_(.+?)\.[0-9]+"', line)
                symbol, ensembl = m.group(1), m.group(2)

                self.mouse_symbols[symbol] = '%s_%s' % (symbol, ensembl)
                self.mouse_ensembls[ensembl] = '%s_%s' % (symbol, ensembl)

        fh.close()

        fh = open(os.path.dirname(inspect.getfile(AlonaBase)) +
                  '/genome/MGI_Gene_Model_Coord.txt.C', 'r')
        next(fh)

        for line in fh:
            gene_id_as_number, ens = line.rstrip('\n').split(' ')

            if gene_id_as_number != 'null':
                self.mouse_entrez[gene_id_as_number] = ens

        fh.close()

    def map_input_genes(self):
        """ Maps gene symbols to internal gene symbols. """
        data = []
        logging.debug('Mapping genes to reference.')

        ftemp = open(self.get_matrix_file() + '.C', 'w')
        with open(self.get_matrix_file(), 'r') as fh:
            if self._has_header:
                header = next(fh)
                ftemp.write(header)

            for line in fh:
                data.append(line.replace('"', ''))

        genes_found = {}
        switch = 0
        total = 0

        # are gene symbols "Entrez"? these gene symbols consists of numbers only.
        is_entrez_gene_id = 1

        for line in data:
            foo = line.split(self._delimiter)
            gene = foo[0]
            if not re.search('^[0-9]+$', gene):
                is_entrez_gene_id = 0

        if is_entrez_gene_id: logging.debug('Gene symbols appear to be Entrez.')

        for line in data:
            total += 1

            foo = line.split(self._delimiter)
            gene = foo[0]

            # some have gene symbols like this: Abc__chr1
            if gene.find('__') > 0:
                gene = gene.split('__')[0]

            if is_entrez_gene_id:
                if self.mouse_entrez.get(gene, '') != '':
                    ens = self.mouse_entrez[gene]

                    if self.mouse_ensembls.get(ens, '') != '':
                        new_gene_name = self.mouse_ensembls.get(ens, '')
                        genes_found[gene] = 1

                        ftemp.write('%s%s%s' % (new_gene_name, delimiter,
                                                delimiter.join(foo[1:])))

                        switch = 1
                    else:
                        self.unmappable.append(gene)
                else:
                    self.unmappable.append(gene)
            elif re.search('^.+_ENSMUSG[0-9]+', gene):
                ensembl_id = re.search('^.+_(ENSMUSG[0-9]+)', gene).group(1)

                if self.mouse_ensembls.get(ensembl_id, '') != '':
                    genes_found[ensembl_id] = 1

                    if re.search(r'^\S+_ENSMUSG[0-9]+\.[0-9]+$', gene):
                        new_gene_name = re.search(r'^(\S+_ENSMUSG[0-9]+)\.[0-9]+$',
                                                  gene).group(1)
                        ftemp.write('%s%s%s' % (new_gene_name, self._delimiter,
                                                self._delimiter.join(foo[1:])))
                    else:
                        ftemp.write(line)

                    switch = 1
                else:
                    self.unmappable.append(gene)
            elif re.search('^ENSMUSG[0-9]+', gene):
                ensembl_id = re.search('^(ENSMUSG[0-9]+)', gene).group(1)

                if self.mouse_ensembls.get(ensembl_id, '') != '':
                    genes_found[ensembl_id] = 1

                    new_gene_name = mouse_ensembls[ensembl_id]
                    ftemp.write('%s%s%s' % (new_gene_name, self._delimiter,
                                            self._delimiter.join(foo[1:])))

                    switch = 1
                else:
                    self.unmappable.append(gene)
            else:
                # some pipelines separate aliases with ;
                foobar = gene.split(';')
                found = 0

                for item in foobar:
                    # try mapping using gene symbols
                    if self.mouse_symbols.get(item, '') != '':
                        genes_found[item] = 1
                        new_gene_name = self.mouse_symbols[item]

                        ftemp.write('%s%s%s' % (new_gene_name, self._delimiter,
                                                self._delimiter.join(foo[1:])))
                        switch = 1
                        found = 1
                        break

                if not found:
                    self.unmappable.append(gene)

        ftemp.close()

        if self.params['species'] == 'human':
            del self.unmappable[:]

            f = open(self.get_working_dir() +
                     '/input_clean.mat.genes_missing_mouse_orthologs', 'r')

            for line in f:
                self.unmappable.append(line.rstrip('\n'))

            f.close()

        if switch:
            os.system('mv %s/input.mat.C %s/input.mat' % (self.get_working_dir(),
                                                          self.get_working_dir()))

        if len(self.unmappable) == 0:
            log_info('All genes were mapped to internal symbols.')
        else:
            with open(self.get_working_dir() + '/unmappable.txt', 'w') as fh:
                for item in self.unmappable:
                    fh.write("%s\n" % item)

            log_info('%s gene(s) were not mappable.' % '{:,}'.
                     format(len(self.unmappable)))
            log_info('These have been written to unmappable.txt')
            
        if (total-len(self.unmappable)) < 500:
            raise TooFewGenesMappableError('Input data error. \
Too few genes were mappable (<500).')

    def check_gene_name_column_id_present(self):
        logging.debug('running check_gene_name_column_id_present()')
        
        with open(self.get_matrix_file(), 'r') as fh:
            header = next(fh)
            line2 = next(fh)

            self._has_gene_id_column_id = len(header.split(self._delimiter)) == \
                                          len(line2.split(self._delimiter))

            if self._has_gene_id_column_id:
                log_info('A column ID for the gene symbols was identified.')
