import os
import subprocess
import sys

import utils
import info


sentinal = utils.AutoSentinal(os.path.join(info.data_directory, 'sentinal_'))


sentinal.run('download_ensembl_gtf', utils.wget_file, info.ensembl_gtf_url, info.ensembl_gtf_filename)
sentinal.run('download_ensembl_chr20', utils.wget_file, info.ensembl_chr20_url, info.ensembl_chr20_filename)


for sra_id in info.sra_samples:

    sentinal.run('download_'+sra_id, utils.download_single_sra_dataset, sra_id, os.path.join(info.data_directory, sra_id))

    sentinal.run('extract_'+sra_id, utils.extract_single_sra_dataset, sra_id, os.path.join(info.data_directory, sra_id))


def create_simulated_rnaseq(sample_id, sim_seed):

    num_fusions = 10
    num_wildtype = 100

    subprocess.check_call([sys.executable, 'rnaseq_fusion_simulator.py',
        info.ensembl_gtf_filename,
        info.ensembl_chr20_filename,
        str(num_fusions),
        str(num_wildtype),
        os.path.join(info.sample_dir(sample_id), 'simulation'),
        info.fastq_filename(sample_id, '1'),
        info.fastq_filename(sample_id, '2'),
        '--chromosomes', '20',
        '--seed', str(sim_seed)])


for sample_id, sim_seed in zip(info.sim_samples, info.sim_seeds):

    sentinal.run('simulate_'+sample_id, create_simulated_rnaseq, sample_id, sim_seed)

