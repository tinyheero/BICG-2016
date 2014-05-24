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

        utils.makedirs(defuse_info.data_directory)

        with utils.CurrentDirectory(defuse_info.packages_directory):

            subprocess.check_call('wget --no-check-certificate http://downloads.sourceforge.net/project/defuse/defuse/0.6/defuse-0.6.1.tar.gz'.split(' '))

            subprocess.check_call('tar -xzvf defuse-0.6.1.tar.gz'.split(' '))

            extract_dir = os.path.join(defuse_info.packages_directory, 'defuse-0.6.1')

            with utils.CurrentDirectory(os.path.join(extract_dir, 'tools')):

                subprocess.check_call(['make'])

            config_filename = os.path.join(defuse_info.data_directory, 'config.txt')
            template_config_filename = os.path.join(extract_dir, 'scripts', 'config.txt')

            with open(config_filename, 'w') as config_file, open(template_config_filename, 'r') as template_config_file:

                for line in template_config_file:
                    if not line.startswith('#') and '=' in line:
                        key, value = line.split('=')
                        key = key.strip()
                        value = value.strip()

                        if key == 'source_directory':
                            line = 'source_director = {0}\n'.format(extract_dir)
                        elif key == 'dataset_directory':
                            line = 'dataset_directory = {0}\n'.format(defuse_info.data_directory)
                        elif '[path of your' in value:
                            line = '{0} = {1}\n'.format(key, value.split(' ')[3])
                    
                    config_file.write(line)
                    
