import os

config = {}
execfile('setup.config', config) 

tutorial_directory = config['tutorial_directory']
chromosomes = config['chromosomes']

data_directory = os.path.join(tutorial_directory, 'data')
install_directory = os.path.join(tutorial_directory, 'install')
results_directory = os.path.join(tutorial_directory, 'results')

sra_samples = ['SRR064439', 'SRR201779']

sim_samples = ['SIM001']
sim_seeds = [2014]

rnaseq_samples = sim_samples


def sample_dir(sample_id):
    return os.path.join(data_directory, sample_id)


def fastq_filename(sample_id, read_end):
    return os.path.join(sample_dir(sample_id), sample_id+'_'+read_end+'.fastq')

