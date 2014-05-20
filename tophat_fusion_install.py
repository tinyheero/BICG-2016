import glob
import shutil
import os
import sys
import subprocess

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

            extract_dir = os.path.join(tophat_fusion_info.install_directory, package_name)

            os.symlink(os.path.join(extract_dir, 'tophat2'), tophat_fusion_info.tophat2_bin)
            os.symlink(os.path.join(extract_dir, 'tophat-fusion-post'), tophat_fusion_info.tophat_fusion_post_bin)

            if 'darwin' in sys.platform:
                package_name = 'cufflinks-2.2.1.OSX_x86_64'
            elif 'linux' in sys.platform:
                package_name = 'cufflinks-2.2.1.Linux_x86_64'

            subprocess.check_call('wget --no-check-certificate http://cufflinks.cbcb.umd.edu/downloads/{0}.tar.gz'.format(package_name).split(' '))

            subprocess.check_call('tar -xzvf {0}.tar.gz'.format(package_name).split(' '))

            extract_dir = os.path.join(tophat_fusion_info.install_directory, package_name)

            os.symlink(os.path.join(extract_dir, 'cufflinks'), tophat_fusion_info.cufflinks_bin)
            os.symlink(os.path.join(extract_dir, 'cuffnorm'), tophat_fusion_info.cuffnorm_bin)
            os.symlink(os.path.join(extract_dir, 'gffread'), tophat_fusion_info.gffread_bin)


with Sentinal('get_data') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(tophat_fusion_info.data_directory)
        utils.makedirs(tophat_fusion_info.data_directory)
        
        for url, filename in ((tophat_fusion_info.ref_gene_url, tophat_fusion_info.ref_gene_filename), (tophat_fusion_info.ens_gene_url, tophat_fusion_info.ens_gene_filename)):

            subprocess.check_call(['wget', url, '-O', filename+'.gz'])

            subprocess.check_call(['gunzip', filename])

        utils.rmtree(tophat_fusion_info.blast_human_directory)
        utils.makedirs(tophat_fusion_info.blast_human_directory)

        with utils.CurrentDirectory(tophat_fusion_info.blast_human_directory):

            subprocess.check_call(['wget', 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genomic.*.tar.gz'])
            subprocess.check_call(['wget', 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz'])

            for filename in glob.glob('*.tar.gz'):
                subprocess.check_call(['tar', '-xzvf', filename])


