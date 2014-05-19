import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'chimerascan')
source_directory = os.path.join(install_directory, 'src')
bin_directory = os.path.join(install_directory, 'bin')
data_directory = os.path.join(install_directory, 'data')
index_directory = os.path.join(install_directory, 'index')


if 'PYTHONPATH' not in os.environ:
    os.environ['PYTHONPATH'] = ''
os.environ['PYTHONPATH'] += os.path.join(install_directory, 'lib', 'python2.7', 'site-packages')


chimerascan_index_bin = os.path.join(bin_directory, 'chimerascan_index.py')
chimerascan_run_bin = os.path.join(bin_directory, 'chimerascan_run.py')


reference_genome_fasta = os.path.join(data_directory, 'hg19.fa')
gene_models_filename = os.path.join(data_directory, 'hg19.ucsc_genes.txt')


def results_directory(sample_id):
    return os.path.join(info.results_directory, 'chimerascan', sample_id)

