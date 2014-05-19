import glob
import shutil
import os
import subprocess

import utils
import info

from tophat_fusion_common import *


def run_tophat(sample_id):

    utils.makedirs(tophat_output_directory(sample_id))

    subprocess.check_call([tophat2_bin,
                           '--mate-inner-dist', '100',
                           '--mate-std-dev', '30',
                           '--fusion-search',
                           '-o', tophat_output_directory(sample_id),
                           index_prefix,
                           info.fastq_filename(sample_id, '1'),
                           info.fastq_filename(sample_id, '2')])


def run_tophat_fusion_post(sample_id):

    with CurrentDirectory(results_directory(sample_id)):

        os.symlink(ref_gene_filename, )



for sample_id in info.rnaseq_samples:
    sentinal.run('run_tophat_'+sample_id, run_tophat, sample_id)


