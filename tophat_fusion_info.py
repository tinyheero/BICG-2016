import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'tophat_fusion')
bin_directory = os.path.join(install_directory, 'bin')
data_directory = os.path.join(install_directory, 'data')


tophat2_bin = os.path.join(bin_directory, 'tophat2')
tophat_fusion_post_bin = os.path.join(bin_directory, 'tophat-fusion-post')


ref_gene_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
ref_gene_filename = os.path.join(data_directory, 'refGene.txt')

ens_gene_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz'
ens_gene_filename = os.path.join(data_directory, 'ensGene.txt')


def results_directory(sample_id):
    return os.path.join(info.results_directory, 'tophat_fusion', sample_id)

def tophat_output_directory(sample_id):
    return os.path.join(results_directory(sample_id), 'tophat_'+sample_id)

