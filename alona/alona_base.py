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
            
        self.params.update(params)
        
        # Set default options
        if self.params['output_directory'] == None:
            self.params['output_directory'] = 'alona_out_%s' % self.random()
        
    def random(self):
        """Get random 8 character string"""
        return str(uuid.uuid4()).split('-')[0]
    
    def create_work_dir(self):        
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
        
        if f.from_file(mat_out) != 'text/plain':
            raise input_not_plain_text('Input file is not plain text.')

        return
        
    def cleanup(self):
        # remove temporary files
        
        for garbage in ('input.mat',):
            logging.debug('removing %s' % garbage)
            os.remove('%s/%s' % self.params['output_directory'],garbage)
