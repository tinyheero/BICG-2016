import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'trinity')

packages_directory = os.path.join(install_directory, 'packages')
bin_directory = os.path.join(install_directory, 'bin')
trinity_bin = os.path.join(install_directory, 'bin', 'trinity', 'Trinity')


results_directory = os.path.join(info.results_directory, 'trinity')

def sample_results_directory(sample_id):
    return os.path.join(results_directory, sample_id)

