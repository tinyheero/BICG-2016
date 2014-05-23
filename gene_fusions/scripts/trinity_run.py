import glob
import shutil
import os
import subprocess

import utils
import info
import trinity_info


Sentinal = utils.Sentinal(os.path.join(trinity_info.results_directory, 'sentinal_'))


for sample_id in info.sim_samples:

    with Sentinal('run_trinity_'+sample_id) as sentinal:

        if sentinal.unfinished:

            subprocess.check_call([trinity_info.trinity_bin,
                                   '--seqType', 'fq',
                                   '--JM', '32G',
                                   '--CPU', '4',
                                   '--left', info.fastq_filename(sample_id, '1'),
                                   '--right', info.fastq_filename(sample_id, '2'),
                                   '--output', trinity_info.sample_results_directory(sample_id)])

