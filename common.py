


sentinal = utils.AutoSentinal(os.path.join(chimerascan_directory, 'sentinal_'))

def get_data():

    with utils.CurrentDirectory(data_directory):

        subprocess.check_call('rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz .'.split(' '))

        subprocess.check_call('tar -zxvf chromFa.tar.gz'.split(' '))

        with open(reference_genome_fasta, 'w') as genome_file:
            for chromosome_filename in glob.glob('chr*.fa'):
                print 'adding ' + chromosome_filename + ' to reference'
                chromosome = chromosome_filename[3:-3]
                if len(chromosome) == 1 or len(chromosome) == 2:
                    with open(chromosome_filename, 'r') as chromosome_file:
                        shutil.copyfileobj(chromosome_file, genome_file)

        subprocess.check_call('wget --no-check-certificate https://chimerascan.googlecode.com/files/hg19.ucsc_genes.txt.gz'.split(' '))

        subprocess.check_call('gunzip hg19.ucsc_genes.txt.gz'.split(' '))

sentinal.run('get_data', get_data)

