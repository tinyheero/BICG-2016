import os
import subprocess
import sys
import tarfile
import shutil

import utils
import info


sentinal = utils.AutoSentinal(os.path.join(info.data_directory, 'sentinal_'))


sentinal.run('download_ensembl_gtf', utils.wget_file, info.ensembl_gtf_url, info.ensembl_gtf_filename)
sentinal.run('download_ensembl_chr20', utils.wget_file, info.ensembl_chr20_url, info.ensembl_chr20_filename)
sentinal.run('download_hg19', utils.wget_file, info.hg19_url, info.hg19_tar_filename)


def create_genome():

    with open(info.hg19_filename, 'w') as hg19_file, tarfile.open(info.hg19_tar_filename, 'r:gz') as tar:

        for tarinfo in tar:

            chromosome = tarinfo.name[3:-3]

            if chromosome in info.chromosomes:
                shutil.copyfileobj(tar.extractfile(tarinfo), hg19_file)

sentinal.run('extract_hg19', create_genome)


def index_genome():

    subprocess.check_call(['bowtie2-build', info.hg19_filename, info.hg19_index_base])

sentinal.run('index_hg19', index_genome)


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

