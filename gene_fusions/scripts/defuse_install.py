import glob
import shutil
import os
import sys
import subprocess

import utils
import info
import defuse_info


Sentinal = utils.Sentinal(os.path.join(defuse_info.install_directory, 'sentinal_'))


with Sentinal('install') as sentinal:

    if sentinal.unfinished:

        utils.rmtree(defuse_info.packages_directory)
        utils.makedirs(defuse_info.packages_directory)

        utils.rmtree(defuse_info.bin_directory)
        utils.makedirs(defuse_info.bin_directory)

        with utils.CurrentDirectory(defuse_info.packages_directory):

            subprocess.check_call('wget --no-check-certificate http://downloads.sourceforge.net/project/defuse/defuse/0.6/defuse-0.6.1.tar.gz'.split(' '))

            subprocess.check_call('tar -xzvf defuse-0.6.1.tar.gz'.split(' '))

            extract_dir = os.path.join(defuse_info.packages_directory, 'defuse-0.6.1')

            with utils.CurrentDirectory(os.path.join(extract_dir, 'tools')):

                subprocess.check_call(['make'])

