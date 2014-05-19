import glob
import shutil
import os
import subprocess

import utils
import info

from tophat_fusion_common import *


def run(sample_id):

    subprocess.check_call([tophat2_bin,
                           index_prefix,
                           info.fastq_filename(sample_id, '1'),
                           info.fastq_filename(sample_id, '2'),
                           '-o', results_directory(sample_id)])


for sample_id in info.rnaseq_samples:
    sentinal.run('run_'+sample_id, run, sample_id)


