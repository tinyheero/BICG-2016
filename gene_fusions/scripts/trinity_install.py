import glob
import shutil
import os
import sys
import subprocess

import utils
import info
import trinity_info


Sentinal = utils.Sentinal(os.path.join(trinity_info.install_directory, 'sentinal_'))


with Sentinal('install') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(trinity_info.packages_directory)
        utils.makedirs(trinity_info.packages_directory)

        utils.rmtree(trinity_info.bin_directory)
        utils.makedirs(trinity_info.bin_directory)

        with utils.CurrentDirectory(trinity_info.packages_directory):

            subprocess.check_call('wget --no-check-certificate http://sourceforge.net/projects/trinityrnaseq/files/trinityrnaseq_r20140413p1.tar.gz'.split(' '))

            subprocess.check_call('tar -xzvf trinityrnaseq_r20140413p1.tar.gz'.split(' '))

            extract_dir = os.path.join(trinity_info.packages_directory, 'trinityrnaseq_r20140413p1')

            with utils.CurrentDirectory(extract_dir):

                subprocess.check_call(['make'])

            try:
                os.remove(os.path.join(trinity_info.bin_directory, 'trinity'))
            except:
                pass
            os.symlink(extract_dir, os.path.join(trinity_info.bin_directory, 'trinity'))


with Sentinal('download_genome_data') as sentinal:

    if sentinal.unfinished:

        for chromosome in info.chromosomes:

            utils.wget_file(trinity_info.ensembl_chromosome_url.format(chromosome), trinity_info.ensembl_chromosome_fasta.format(chromosome)+'.gz')
            subprocess.check_call(['gunzip', trinity_info.ensembl_chromosome_fasta.format(chromosome)+'.gz'])


with Sentinal('prepare_genome_data') as sentinal:

    if sentinal.unfinished:

        with open(trinity_info.ensembl_genome_fasta, 'w') as genome_file:

            for chromosome in info.chromosomes:

                with open(trinity_info.ensembl_chromosome_fasta.format(chromosome), 'r') as chromosome_file:

                    shutil.copyfileobj(chromosome_file, genome_file)


with Sentinal('gmap_build') as sentinal:

    if sentinal.unfinished:

        utils.makedirs(trinity_info.gmap_index_directory)
        subprocess.check_call(['gmap_build', '-D', trinity_info.gmap_index_directory, '-d', 'hg19', trinity_info.ensembl_genome_fasta])


with Sentinal('download_gtf_data') as sentinal:

    if sentinal.unfinished:

        temp_gtf_filename = trinity_info.ensembl_gtf_filename+'.tmp'

        utils.wget_file(trinity_info.ensembl_gtf_url, temp_gtf_filename+'.gz')
        subprocess.check_call(['gunzip', temp_gtf_filename+'.gz'])

        with open(temp_gtf_filename, 'r') as temp_gtf_file, open(trinity_info.ensembl_gtf_filename, 'w') as gtf_file:

            for line in temp_gtf_file:

                chromosome = line.split()[0]

                if chromosome in info.chromosomes:
                    gtf_file.write(line)

        os.remove(temp_gtf_filename)

