import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'trinity')

packages_directory = os.path.join(install_directory, 'packages')
data_directory = os.path.join(install_directory, 'data')
bin_directory = os.path.join(install_directory, 'bin')
trinity_bin = os.path.join(install_directory, 'bin', 'trinity', 'Trinity')

gmap_index_directory = os.path.join(data_directory, 'gmap')

ensembl_chromosome_url = 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.{0}.fa.gz'
ensembl_chromosome_fasta = os.path.join(data_directory, 'Homo_sapiens.GRCh37.75.dna.chromosome.{0}.fa')
ensembl_genome_fasta = os.path.join(data_directory, 'Homo_sapiens.GRCh37.75.dna.chromosomes.fa')

ensembl_gtf_url = 'ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz'
ensembl_gtf_filename = os.path.join(data_directory, 'Homo_sapiens.GRCh37.75.gtf')
