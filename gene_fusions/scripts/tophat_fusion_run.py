import glob
import shutil
import os
import subprocess

import utils
import info
import tophat_fusion_info


os.environ['PATH'] += ':' + tophat_fusion_info.bin_directory


Sentinal = utils.Sentinal(os.path.join(tophat_fusion_info.results_directory, 'sentinal_'))


for sample_id in info.rnaseq_samples:

    with Sentinal('run_tophat_'+sample_id) as sentinal:

        if sentinal.unfinished:

            utils.makedirs(tophat_fusion_info.sample_tophat_output_directory(sample_id))

            subprocess.check_call(['tophat2',
                                   '--mate-inner-dist', '100',
                                   '--mate-std-dev', '30',
                                   '--fusion-search',
                                   '-o', tophat_fusion_info.sample_tophat_output_directory(sample_id),
                                   info.hg19_index_base,
                                   info.fastq_filename(sample_id, '1'),
                                   info.fastq_filename(sample_id, '2')])

    with Sentinal('run_tophat_fusion_'+sample_id) as sentinal:

        if sentinal.unfinished:

            with utils.CurrentDirectory(tophat_fusion_info.sample_results_directory(sample_id)):

                utils.symlink(tophat_fusion_info.blast_human_directory)
                utils.symlink(tophat_fusion_info.ref_gene_filename)
                utils.symlink(tophat_fusion_info.ens_gene_filename)

                subprocess.check_call(['tophat-fusion-post',
                                       '--num-fusion-reads', '1',
                                       '--num-fusion-pairs', '2',
                                       '--num-fusion-both', '5',
                                       info.hg19_index_base])

