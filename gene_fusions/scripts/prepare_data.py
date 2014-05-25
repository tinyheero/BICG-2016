import os
import subprocess
import sys
import tarfile
import shutil

import utils
import info


utils.makedirs(info.data_directory)


Sentinal = utils.Sentinal(os.path.join(info.data_directory, 'sentinal_'))


for sra_id in info.sra_samples:

    with Sentinal('download_'+sra_id) as sentinal:

        if sentinal.unfinished:

            utils.download_single_sra_dataset(sra_id, os.path.join(info.data_directory, sra_id))

    with Sentinal('extract_'+sra_id) as sentinal:

        if sentinal.unfinished:

            utils.extract_single_sra_dataset(sra_id, os.path.join(info.data_directory, sra_id))


with Sentinal('download_ensembl_gtf') as sentinal:

    if sentinal.unfinished:

        utils.wget_file(info.ensembl_gtf_url, info.ensembl_gtf_filename+'.gz')
        subprocess.check_call(['gunzip', info.ensembl_gtf_filename+'.gz'])


with Sentinal('download_ensembl_chr20') as sentinal:

    if sentinal.unfinished:

        utils.wget_file(info.ensembl_chr20_url, info.ensembl_chr20_filename+'.gz')
        subprocess.check_call(['gunzip', info.ensembl_chr20_filename+'.gz'])


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

    with Sentinal('simulate_'+sample_id) as sentinal:

        if sentinal.unfinished:

            create_simulated_rnaseq(sample_id, sim_seed)

