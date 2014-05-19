import glob
import shutil
import os
import subprocess

import utils
import info

from chimerascan_common import *



def install():

    utils.rmtree(source_directory)
    
    with utils.CurrentDirectory(source_directory):

        subprocess.check_call('wget --no-check-certificate https://chimerascan.googlecode.com/files/chimerascan-0.4.5a.tar.gz'.split(' '))

        subprocess.check_call('tar -xzvf chimerascan-0.4.5a.tar.gz'.split(' '))

        os.chdir('chimerascan-0.4.5')

        # Patch the code
        with open('chimerascan/pysam/samtools/ksort.h', 'r') as ksort_file:
            ksort_code = ksort_file.read()
        ksort_code = ksort_code.replace('inline', 'static inline')
        with open('chimerascan/pysam/samtools/ksort.h', 'w') as ksort_file:
            ksort_file.write(ksort_code)

        subprocess.check_call('python setup.py install --prefix'.split(' ') + [install_directory])

sentinal.run('install', install)


def get_data():

    with utils.CurrentDirectory(data_directory):

        subprocess.check_call('rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz .'.split(' '))

        subprocess.check_call('tar -zxvf chromFa.tar.gz'.split(' '))

        with open(reference_genome_fasta, 'w') as genome_file:
            for chromosome_filename in glob.glob('chr*.fa'):
                print 'adding ' + chromosome_filename + ' to reference'
                chromosome = chromosome_filename[3:-3]
                if len(chromosome) == 1 or len(chromosome) == 2:
                    with open(chromosome_filename, 'r') as chromosome_file:
                        shutil.copyfileobj(chromosome_file, genome_file)

        subprocess.check_call('wget --no-check-certificate https://chimerascan.googlecode.com/files/hg19.ucsc_genes.txt.gz'.split(' '))

        subprocess.check_call('gunzip hg19.ucsc_genes.txt.gz'.split(' '))

sentinal.run('get_data', get_data)


def chimerascan_index():

    utils.rmtree(index_directory)

    subprocess.check_call(['python', chimerascan_index_bin, reference_genome_fasta, gene_models_filename, index_directory])

sentinal.run('chimerascan_index', chimerascan_index)


