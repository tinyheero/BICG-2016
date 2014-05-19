import glob
import shutil
import os
import subprocess

import utils
import info
import chimerascan_info


Sentinal = utils.Sentinal(os.path.join(chimerascan_info.results_directory, 'sentinal_'))


for sample_id in info.rnaseq_samples:

    with Sentinal('run_'+sample_id) as sentinal:

        if sentinal.unfinished:

            subprocess.check_call([chimerascan_info.chimerascan_run_bin, 
                                   chimerascan_info.index_directory,
                                   info.fastq_filename(sample_id, '1'),
                                   info.fastq_filename(sample_id, '2'),
                                   chimerascan_info.sample_results_directory(sample_id)])

