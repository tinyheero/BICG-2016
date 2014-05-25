import glob
import shutil
import os
import subprocess

import utils
import info
import trinity_info


trinity_results_directory = os.path.join(info.results_directory, 'trinity')


Sentinal = utils.Sentinal(os.path.join(trinity_results_directory, 'sentinal_'))


for sample_id in info.rnaseq_samples:

    results_directory = os.path.join(trinity_results_directory, sample_id)

    with Sentinal('run_trinity_'+sample_id) as sentinal:

        if sentinal.unfinished:

            subprocess.check_call([trinity_info.trinity_bin,
                                   '--seqType', 'fq',
                                   '--JM', '32G',
                                   '--CPU', '4',
                                   '--left', info.fastq_filename(sample_id, '1'),
                                   '--right', info.fastq_filename(sample_id, '2'),
                                   '--output', results_directory])

    contig_fasta = os.path.join(results_directory, 'Trinity.fasta')
    contig_mapping_filename = os.path.join(results_directory, 'Trinity.gff3')

    with Sentinal('run_gmap_'+sample_id) as sentinal:

        if sentinal.unfinished:

            with open(contig_mapping_filename, 'w') as contig_mapping_file:

                subprocess.check_call(['gmap',
                                       '-f', 'gff3_gene',
                                       '-D', trinity_info.gmap_index_directory,
                                       '-d', 'hg19',
                                       contig_fasta],
                                      stdout=contig_mapping_file)
