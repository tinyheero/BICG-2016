import glob
import shutil
import os
import subprocess

import utils
import info
import defuse_info


Sentinal = utils.Sentinal(os.path.join(defuse_directory, 'sentinal_'))


def results_directory(sample_id):
    return os.path.join(defuse_directory, sample_id)


for sample_id in info.rnaseq_samples:

    with Sentinal('run_'+sample_id) as sentinal:

        if sentinal.unfinished:

            utils.makedirs(results_directory(sample_id))
            subprocess.check_call([defuse_info.defuse_script,
                                   '-c', defuse_info.config_filename,
                                   '-1', info.fastq_filename(sample_id, '1'),
                                   '-2', info.fastq_filename(sample_id, '2'),
                                   '-o', results_directory(sample_id)])

