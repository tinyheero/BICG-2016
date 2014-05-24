import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'defuse')

packages_directory = os.path.join(install_directory, 'packages')
bin_directory = os.path.join(install_directory, 'bin')
data_directory = os.path.join(install_directory, 'data')

createref_script = os.path.join(bin_directory, 'scripts', 'create_reference_dataset.pl')
defuse_script = os.path.join(bin_directory, 'scripts', 'defuse.pl')
config_filename = os.path.join(data_directory, 'config.txt')

