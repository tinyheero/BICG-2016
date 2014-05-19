import os

import info
import utils


tophat_fusion_directory = os.path.join(info.install_directory, 'tophat_fusion')
bin_directory = os.path.join(tophat_fusion_directory, 'bin')
index_directory = os.path.join(tophat_fusion_directory, 'index')


sentinal = utils.AutoSentinal(os.path.join(tophat_fusion_directory, 'sentinal_'))


tophat2_bin = os.path.join(bin_directory, 'tophat2')
tophat_fusion_post_bin = os.path.join(bin_directory, 'tophat-fusion-post')

index_prefix = '/Users/amcphers/Downloads/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'#os.path.join(index_directory, 'hg19')


def results_directory(sample_id):
    return os.path.join(info.results_directory, 'tophat_fusion', sample_id)

