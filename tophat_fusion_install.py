import glob
import shutil
import os
import subprocess

import utils
import info
import tophat_fusion_info


Sentinal = utils.Sentinal(os.path.join(tophat_fusion_info.install_directory, 'sentinal_'))


with Sentinal('install') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(tophat_fusion_info.install_directory)
        
        utils.makedirs(tophat_fusion_info.install_directory)
        
        with utils.CurrentDirectory(tophat_fusion_info.install_directory):

            subprocess.check_call('wget --no-check-certificate http://tophat.cbcb.umd.edu/downloads/tophat-2.0.11.OSX_x86_64.tar.gz'.split(' '))

            subprocess.check_call('tar -xzvf tophat-2.0.11.OSX_x86_64.tar.gz'.split(' '))

        extract_dir = os.path.join(tophat_fusion_info.install_directory, 'tophat-2.0.11.OSX_x86_64')

        utils.makedirs(tophat_fusion_info.bin_directory)

        os.symlink(os.path.join(extract_dir, 'tophat2'), tophat2_bin)
        os.symlink(os.path.join(extract_dir, 'tophat-fusion-post'), tophat_fusion_post_bin)


with Sentinal('get_data') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(tophat_fusion_info.data_directory)

        utils.makedirs(tophat_fusion_info.data_directory)
        
        for url, filename in ((tophat_fusion_info.ref_gene_url, tophat_fusion_info.ref_gene_filename), (tophat_fusion_info.ens_gene_url, tophat_fusion_info.ens_gene_filename)):

            subprocess.check_call(['wget', url, '-O', filename+'.gz'])

            subprocess.check_call(['gunzip', filename])

