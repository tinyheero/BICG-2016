import os

config = {}
execfile('paths.conf', config) 

tutorial_directory = config['tutorial_directory']

data_directory = os.path.join(tutorial_directory, 'data')
genome_directory = os.path.join(tutorial_directory, 'genome')
install_directory = os.path.join(tutorial_directory, 'install')
results_directory = os.path.join(tutorial_directory, 'results')

ensembl_gtf_url = 'ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz'
ensembl_chr20_url = 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.20.fa.gz'

ensembl_gtf_filename = os.path.join(genome_directory, 'Homo_sapiens.GRCh37.75.gtf.gz')
ensembl_chr20_filename = os.path.join(genome_directory, 'Homo_sapiens.GRCh37.75.dna.chromosome.20.fa.gz')

hg19_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'
hg19_tar_filename = os.path.join(genome_directory, 'chromFa.tar.gz')
hg19_filename = os.path.join(genome_directory, 'hg19.fa')
hg19_index_base = os.path.join(genome_directory, 'hg19')

chromosomes = [str(a) for a in xrange(1, 23)] + ['X', 'Y', 'MT']

sra_samples = ['SRR064439', 'SRR201779']

sim_samples = ['SIM001']
sim_seeds = [2014]

rnaseq_samples = sra_samples + sim_samples


def sample_dir(sample_id):
    return os.path.join(data_directory, sample_id)


def fastq_filename(sample_id, read_end):
    return os.path.join(sample_dir(sample_id), sample_id+'_'+read_end+'.fastq')
