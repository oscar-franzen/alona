import re
import os
import uuid
import magic
import logging
import subprocess

from .exceptions import *

class alona_base(object):
    """
    alona_base class
    """
    
    params = None

    def __init__(self, params=None):
        if params is None:
            params = {}
            
        self.params = {
        }
        
        self._is_binary = None
        self._delimiter = None
        self._has_header = None
            
        self.params.update(params)
        
        # Set default options
        if self.params['output_directory'] == None:
            self.params['output_directory'] = 'alona_out_%s' % self.random()
            
        if self.params['delimiter'] == None:
            self.params['auto']
            
        if self.params['loglevel'] == 'debug':
            logging.debug('*** parameters *********************')
            for par in self.params:
                logging.debug('%s : %s' % (par,self.params[par]))
            logging.debug('************************************')
        
    def random(self):
        """Get random 8 character string"""
        return str(uuid.uuid4()).split('-')[0]
    
    def create_work_dir(self):
        """ Creates a working directory for temporary and output files. """
        try:
            logging.debug('creating output directory: %s '% self.params['output_directory'])
            os.mkdir(self.params['output_directory'])
        except FileExistsError:
            logging.error('Error: Output directory already exists (%s)' %
                           self.params['output_directory'])
            raise
        
        return
        
    def is_file_empty(self):
        if os.stat(self.params['input_filename']).st_size == 0:
            raise file_empty_error()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return self
        
    def is_binary(self, filename):
        with open(filename, 'rb') as f:
            for block in f:
                if b'\0' in block:
                    return True
        
        return False

    def unpack_data(self):
        """ Unpacks compressed data and if no compression is used, symlinks to data. """
        # Is the file binary?
        self._is_binary = self.is_binary(self.params['input_filename'])
        abs_path = os.path.abspath(self.params['input_filename'])
        mat_out = '%s/input.mat' % self.params['output_directory']
        
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
                    if re.search('unexpected end of file',exc.output.decode('ascii')):
                        raise file_corrupt('Error: input file is corrupt.')
                
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
                  if re.search(' unexpected end of file',exc.output.decode('ascii')):
                      raise file_corrupt('Error: input file is corrupt.')
                  
                # check that the archive only contains one file
                out = subprocess.check_output('unzip -v %s | wc -l' % (abs_path),
                                               shell=True)
                
                out = out.decode('ascii').replace('\n','')
                no_files = int(out) - 5
                
                if no_files != 1:
                    raise Exception('More than one file in input archive.')
                
                # Files  created  by  zip  can be uncompressed by gzip only if they have
                # a single member compressed with the 'deflation' method.
                os.system('zcat %s > %s' % (abs_path,mat_out))
            elif re.search(' bzip2 compressed data,', out):
                logging.debug('bzip2 data detected.')
                
                # confirm integrity
                try:
                    out = subprocess.check_output('bzip2 -t %s' % (abs_path), shell=True,
                                                   stderr=subprocess.STDOUT)
                    
                    out = out.decode('ascii')
                except subprocess.CalledProcessError as exc:
                    if re.search('file ends unexpectedly',exc.output.decode('ascii')):
                        raise file_corrupt('Input file is corrupt.')
                        
                # uncompress
                os.system('bzip2 -d -c %s > %s' % (abs_path, mat_out))
            else:
                raise invalid_file_format('Invalid format of the input file.')
                
        else:
            logging.debug('Input file is not binary.')
            
            # Create a symlink to the data
            cmd = 'ln -sfn %s %s' % (abs_path,mat_out)
            os.system(cmd)
            
        f = magic.Magic(mime=True)
        
        # don't use `from_file` (it doesn't follow symlinks)
        f_type = f.from_buffer(open(mat_out,'r').read(1024))
        
        if f_type != 'text/plain':
            raise input_not_plain_text('Input file is not plain text (found type=%s).' %
                                        f_type)

        return
        
    def cleanup(self):
        """ Removes temporary files. """
        
        # remove temporary files
        for garbage in ('input.mat',):
            logging.debug('removing %s' % garbage)
            os.remove('%s/%s' % (self.params['output_directory'],garbage))

    def __guess_delimiter(self):
        d = {}
        d[' '] = 0
        d['\t'] = 0
        d[','] = 0

        i = 0
        f = open('%s/input.mat' % self.params['output_directory'],'r')

        for line in f:
            d[' '] += line.count(' ')
            d['\t'] += line.count('\t')
            d[','] += line.count(',')

            i += 1

            if i > 10:
                break

        f.close()

        d_sorted = sorted(d, key=d.get, reverse=True)

        return d_sorted[0]

    def get_delimiter(self):
        """ Figures out the data delimiter of the input data. """
        d = ''
        
        if self.params['delimiter'] == 'auto':
            d = self.__guess_delimiter()
        else:
            d = self.params['delimiter'].upper()
          
            if d == 'TAB':
                d = '\t'
            elif d == 'SPACE':
                d = ' '
                
        logging.debug('delimiter is: "%s" (ASCII code=%s)' % (d,ord(d)))
        
        _delimiter = d
        
        return d
        
    def has_header(self):
        """ Figures out if the input data uses a header or not. """
        
        ret = None
        
        if self.params['header'] == 'auto':
            f = open('%s/input.mat' % self.params['output_directory'],'r')
            
            first_line = next(f).replace('\n','')
            
            total = 0
            count_digit = 0
            
            for item in first_line.split(self._delimiter):
                if item.replace('.','',1).isdigit():
                    count_digit += 1
                
                total += 1
              
            f.close()
            
            ret = not (total == count_digit)
        else:
            ret = (self.params['header'] == 'yes')
            
        # if all fields are non-numerical, it's likely a header
        self._has_header = ret
        
        logging.debug('has header: %s' % self._has_header)
        
        return ret

    def column_sanity_check(self):
        f = open('%s/input.mat' % self.params['output_directory'],'r')

        if self._has_header:
            next(f)

        cols = {}
          
        for line in f:
            cols[len(line.split(self._delimiter))] = 1
          
        f.close()
        
        if len(cols.keys()) > 1:
            raise irregular_column_count('Rows in your data matrix have different number \
of columns (every row must have the same number of columns).')
