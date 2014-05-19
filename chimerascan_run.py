import glob
import shutil
import os
import subprocess

import utils
import info

from chimerascan_common import *


def run(sample_id):
    subprocess.check_call([chimerascan_run_bin, index_directory, info.fastq_filename(sample_id, '1'), info.fastq_filename(sample_id, '2'), results_directory(sample_id)])


for sample_id in info.rnaseq_samples:
    sentinal.run('run_'+sample_id, run, sample_id)


