import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'defuse')

packages_directory = os.path.join(install_directory, 'packages')
data_directory = os.path.join(install_directory, 'data')

createref_script = os.path.join(packages_directory, 'defuse', 'scripts', 'create_reference_dataset.pl')
defuse_script = os.path.join(packages_directory, 'defuse', 'scripts', 'defuse.pl')
config_filename = os.path.join(data_directory, 'config.txt')

