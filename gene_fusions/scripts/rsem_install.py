import glob
import shutil
import os
import sys
import subprocess

import utils
import info
import rsem_info


Sentinal = utils.Sentinal(os.path.join(rsem_info.install_directory, 'sentinal_'))


with Sentinal('install') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(rsem_info.packages_directory)
        utils.makedirs(rsem_info.packages_directory)

        utils.rmtree(rsem_info.bin_directory)
        utils.makedirs(rsem_info.bin_directory)

        with utils.CurrentDirectory(rsem_info.packages_directory):

            subprocess.check_call('wget --no-check-certificate http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.12.tar.gz'.split(' '))

            subprocess.check_call('tar -xzvf rsem-1.2.12.tar.gz'.split(' '))

            extract_dir = os.path.join(rsem_info.packages_directory, 'rsem-1.2.12')

            with utils.CurrentDirectory(extract_dir):

                subprocess.check_call(['make'])

            with utils.CurrentDirectory(rsem_info.bin_directory):

                utils.symlink(os.path.join(extract_dir, 'rsem-prepare-reference'))
                utils.symlink(os.path.join(extract_dir, 'rsem-calculate-expression'))


