import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'tophat_fusion')

packages_directory = os.path.join(install_directory, 'packages')
bin_directory = os.path.join(install_directory, 'bin')
data_directory = os.path.join(install_directory, 'data')
blast_human_directory = os.path.join(data_directory, 'blast_human')


ref_gene_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
ref_gene_filename = os.path.join(data_directory, 'refGene.txt')

ens_gene_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz'
ens_gene_filename = os.path.join(data_directory, 'ensGene.txt')


results_directory = os.path.join(info.results_directory, 'tophat_fusion')

def sample_results_directory(sample_id):
    return os.path.join(results_directory, sample_id)

def sample_tophat_output_directory(sample_id):
    return os.path.join(sample_results_directory(sample_id), 'tophat_'+sample_id)

