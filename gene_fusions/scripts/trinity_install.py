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




