import glob
import shutil
import os
import sys
import subprocess
import tarfile

import utils
import info
import tophat_fusion_info


Sentinal = utils.Sentinal(os.path.join(tophat_fusion_info.install_directory, 'sentinal_'))


with Sentinal('install') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(tophat_fusion_info.packages_directory)
        utils.makedirs(tophat_fusion_info.packages_directory)

        utils.rmtree(tophat_fusion_info.bin_directory)
        utils.makedirs(tophat_fusion_info.bin_directory)

        with utils.CurrentDirectory(tophat_fusion_info.packages_directory):

            if 'darwin' in sys.platform:
                package_name = 'tophat-2.0.11.OSX_x86_64'
            elif 'linux' in sys.platform:
                package_name = 'tophat-2.0.11.Linux_x86_64'

            subprocess.check_call('wget --no-check-certificate http://tophat.cbcb.umd.edu/downloads/{0}.tar.gz'.format(package_name).split(' '))

            subprocess.check_call('tar -xzvf {0}.tar.gz'.format(package_name).split(' '))

            extract_dir = os.path.join(tophat_fusion_info.packages_directory, package_name)

            with utils.CurrentDirectory(tophat_fusion_info.bin_directory):

                utils.symlink(os.path.join(extract_dir, 'tophat2'))
                utils.symlink(os.path.join(extract_dir, 'tophat-fusion-post'))


with Sentinal('download_genome_data') as sentinal:

    if sentinal.unfinished:

        utils.wget_file(tophat_fusion_info.ensembl_gtf_url, tophat_fusion_info.ensembl_gtf_filename+'.gz')
        subprocess.check_call(['gunzip', tophat_fusion_info.ensembl_gtf_filename+'.gz'])

        utils.wget_file(tophat_fusion_info.ucsc_genome_url, tophat_fusion_info.ucsc_genome_tar_filename)


with Sentinal('prepare_genome_data') as sentinal:

    if sentinal.unfinished:

        with open(tophat_fusion_info.ucsc_genome_fasta, 'w') as fasta_file, tarfile.open(tophat_fusion_info.ucsc_genome_tar_filename, 'r:gz') as tar:

            for tarinfo in tar:

                chromosome = tarinfo.name[3:-3]

                if chromosome in info.chromosomes:
                    shutil.copyfileobj(tar.extractfile(tarinfo), fasta_file)

        with open(tophat_fusion_info.ensembl_gtf_filename, 'r') as ensembl_gtf_file, open(tophat_fusion_info.ucsc_gtf_filename, 'w') as ucsc_gtf_file:

            for line in ensembl_gtf_file:

                chromosome = line.split()[0]

                if chromosome in info.chromosomes:
                    ucsc_gtf_file.write('chr'+line)


with Sentinal('bowtie_index') as sentinal:

    if sentinal.unfinished:

        subprocess.check_call(['bowtie-build', tophat_fusion_info.ucsc_genome_fasta, tophat_fusion_info.ucsc_genome_base])
        subprocess.check_call(['bowtie2-build', tophat_fusion_info.ucsc_genome_fasta, tophat_fusion_info.ucsc_genome_base])


with Sentinal('get_data') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(tophat_fusion_info.data_directory)
        utils.makedirs(tophat_fusion_info.data_directory)
        
        for url, filename in ((tophat_fusion_info.ref_gene_url, tophat_fusion_info.ref_gene_filename), (tophat_fusion_info.ens_gene_url, tophat_fusion_info.ens_gene_filename)):

            subprocess.check_call(['wget', url, '-O', filename+'.gz'])

            subprocess.check_call(['gunzip', filename])

        utils.rmtree(tophat_fusion_info.blast_directory)
        utils.makedirs(tophat_fusion_info.blast_directory)

        with utils.CurrentDirectory(tophat_fusion_info.blast_directory):

            subprocess.check_call(['wget', 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genomic.*.tar.gz'])
            subprocess.check_call(['wget', 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz'])

            for filename in glob.glob('*.tar.gz'):
                subprocess.check_call(['tar', '-xzvf', filename])


