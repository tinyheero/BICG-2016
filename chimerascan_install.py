import glob
import shutil
import os
import subprocess

import utils
import info
import chimerascan_info


Sentinal = utils.Sentinal(os.path.join(chimerascan_info.install_directory, 'sentinal_'))


with Sentinal('install') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(chimerascan_info.source_directory)
        
        with utils.CurrentDirectory(chimerascan_info.source_directory):

            subprocess.check_call('wget --no-check-certificate https://chimerascan.googlecode.com/files/chimerascan-0.4.5a.tar.gz'.split(' '))

            subprocess.check_call('tar -xzvf chimerascan-0.4.5a.tar.gz'.split(' '))

            os.chdir('chimerascan-0.4.5')

            # Patch the code
            with open('chimerascan/pysam/samtools/ksort.h', 'r') as ksort_file:
                ksort_code = ksort_file.read()
            ksort_code = ksort_code.replace('inline', 'static inline')
            with open('chimerascan/pysam/samtools/ksort.h', 'w') as ksort_file:
                ksort_file.write(ksort_code)

            subprocess.check_call('python setup.py install --prefix'.split(' ') + [chimerascan_info.install_directory])


with Sentinal('get_data') as sentinal:

    if sentinal.unfinished:

        with utils.CurrentDirectory(chimerascan_info.data_directory):

            os.symlink(info.hg19_filename, chimerascan_info.reference_genome_fasta)

            subprocess.check_call('wget --no-check-certificate https://chimerascan.googlecode.com/files/hg19.ucsc_genes.txt.gz'.split(' '))

            subprocess.check_call('gunzip hg19.ucsc_genes.txt.gz'.split(' '))


with Sentinal('chimerascan_index') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(index_directory)

        subprocess.check_call([sys.executable,
                               chimerascan_info.chimerascan_index_bin,
                               chimerascan_info.reference_genome_fasta,
                               chimerascan_info.gene_models_filename,
                               chimerascan_info.index_directory])

