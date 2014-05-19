import glob
import shutil
import os
import subprocess

import utils
import info

from tophat_fusion_common import *



def install():

    utils.rmtree(tophat_fusion_directory)
    
    with utils.CurrentDirectory(tophat_fusion_directory):

        subprocess.check_call('wget --no-check-certificate http://tophat.cbcb.umd.edu/downloads/tophat-2.0.11.OSX_x86_64.tar.gz'.split(' '))

        subprocess.check_call('tar -xzvf tophat-2.0.11.OSX_x86_64.tar.gz'.split(' '))

    extract_dir = os.path.join(tophat_fusion_directory, 'tophat-2.0.11.OSX_x86_64')

    utils.makedirs(bin_directory)

    os.symlink(os.path.join(extract_dir, 'tophat2'), tophat2_bin)
    os.symlink(os.path.join(extract_dir, 'tophat-fusion-post'), tophat_fusion_post_bin)

sentinal.run('install', install)

def get_data():

    # utils.rmtree(index_directory)
    
    with utils.CurrentDirectory(index_directory):

        # subprocess.check_call('wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip'.split(' '))

        # subprocess.check_call('unzip hg19.zip'.split(' '))

        # genome_fasta_filename = index_prefix+'.fa'
        # with open(genome_fasta_filename, 'w') as genome_fasta_file:
        #     subprocess.check_call(['bowtie2-inspect', index_prefix], stdout=genome_fasta_file)

        # subprocess.check_call('wget http://tophat.cbcb.umd.edu/downloads/test_data.tar.gz'.split(' '))

        # subprocess.check_call('tar zxvf test_data.tar.gz'.split(' '))

        subprocess.check_call([tophat2_bin, '-r', '20', index_prefix, os.path.join('test_data', 'reads_1.fq'), os.path.join('test_data', 'reads_2.fq')])

    utils.rmtree(data_directory)
    
    for url, filename in ((ref_gene_url, ref_gene_filename), (ens_gene_url, ens_gene_filename)):

        subprocess.check_call(['wget', url, '-O', filename+'.gz')

        subprocess.check_call(['gunzip', filename])

sentinal.run('get_data', get_data)


