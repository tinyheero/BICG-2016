import os
import argparse
import gzip
import random
import string
import subprocess
import shutil
import pandas as pd


min_coverage = 3
max_coverage = 100
min_length = 2000

def get_opener(filename):

    if filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    return opener


def read_gene_regions(gene_models_filename, chromosome_subset=None):
    
    gene_regions = list()
    
    with get_opener(gene_models_filename)(gene_models_filename, 'r') as gtf:
        
        for line in gtf:
            
            if line.startswith('#'):
                continue
                
            fields = line.rstrip().split('\t')
            
            if fields[2] != 'exon':
                continue
                
            keyvals = fields[-1].replace('; ', ';').rstrip(';').split(';')
            keyvals = [a.split(' ', 1) for a in keyvals]
            keyvals = [(a, b.strip('"').rstrip('"')) for a, b in keyvals]
            keyvals = dict(keyvals)
            
            keyvals['chromosome'] = fields[0]
            keyvals['start'] = int(fields[3])
            keyvals['end'] = int(fields[4])
            keyvals['strand'] = fields[6]
            
            if chromosome_subset is not None and keyvals['chromosome'] not in chromosome_subset:
                continue

            gene_regions.append(keyvals)
            
    gene_regions = pd.DataFrame(gene_regions)

    gene_regions['exon_number'] = gene_regions['exon_number'].astype(int)
    
    return gene_regions


def read_sequences(fasta):

    id = None
    sequences = []

    for line in fasta:

        line = line.rstrip()

        if len(line) == 0:
            continue

        if line[0] == '>':
            if id is not None:
                yield (id, ''.join(sequences))
            id = line[1:]
            sequences = []
        else:
            sequences.append(line)

    if id is not None:
        yield (id, ''.join(sequences))


def read_genome(genome_filename, chromosome_subset=None):

    genome = dict()

    with get_opener(genome_filename)(genome_filename, 'r') as genome_file:

        for chromosome, sequence in read_sequences(genome_file):

            chromosome = chromosome.split()[0]

            if chromosome_subset is None or chromosome in chromosome_subset:
                genome[chromosome] = sequence

    return genome


def reverse_complement(sequence):
    return sequence[::-1].translate(string.maketrans('ACTGactg','TGACtgac'))


def create_gene_sequence(genome, exons):

    sequence = ''

    for idx, row in exons.iterrows():

        exon_sequence = genome[row['chromosome']][row['start']-1:row['end']-1]

        if row['strand'] == '-':
            exon_sequence = reverse_complement(exon_sequence)

        sequence += exon_sequence

    return sequence


def simulate_reads(seq_id, sequence, coverage, output_dir, fragment_mean, fragment_stddev, read_length):

    sim_dir = os.path.join(output_dir, seq_id)
    sim_prefix = os.path.join(sim_dir, 'dgwsim')

    try:
        os.makedirs(sim_dir)
    except OSError as e:
        if e.errno != 17:
            raise

    fasta_filename = os.path.join(sim_dir, 'sequence.fa')

    with open(fasta_filename, 'w') as fasta_file:
        fasta_file.write('>{0}\n{1}'.format(seq_id, sequence))

    read_count = int(float(coverage) * float(len(sequence)) / float(fragment_mean))

    dwgsim_stdout_filename = os.path.join(sim_dir, 'dwgsim.out')
    dwgsim_stderr_filename = os.path.join(sim_dir, 'dwgsim.err')

    with open(dwgsim_stdout_filename, 'w') as dwgsim_stdout_file, open(dwgsim_stderr_filename, 'w') as dwgsim_stderr_file:
        subprocess.check_call(['dwgsim', 
                               '-H', 
                               '-y', '0',
                               '-z', str(random.randint(0, 1000000)),
                               '-N', str(int(read_count)),
                               '-d', str(int(fragment_mean)),
                               '-s', str(int(fragment_stddev)),
                               '-1', str(int(read_length)),
                               '-2', str(int(read_length)),
                               fasta_filename,
                               sim_prefix],
                              stdout=dwgsim_stdout_file,
                              stderr=dwgsim_stderr_file)

    return sim_prefix + '.bwa.read1.fastq', sim_prefix + '.bwa.read2.fastq'


def cat_files(out_filename, in_filenames):

    with open(out_filename, 'w') as out_file:

        for in_filename in in_filenames:

            with open(in_filename, 'r') as in_file:
                shutil.copyfileobj(in_file, out_file)


def simulate(args):

    try:
        os.makedirs(args.outdir)
    except OSError as e:
        if e.errno != 17:
            raise

    random.seed(args.seed)

    genome = read_genome(args.genome, args.chromosomes)

    genes = read_gene_regions(args.gtf)

    if args.chromosomes is not None:
        genes = genes[genes['chromosome'].isin(args.chromosomes)]

    multi_exon_ids = genes.loc[genes['exon_number'] != 1, 'transcript_id'].unique()
    genes = genes[genes['transcript_id'].isin(multi_exon_ids)]

    assert len(genes.index) > 0

    fusions_table = list()

    reads1_filenames = list()
    reads2_filenames = list()

    while len(fusions_table) < args.num_fusions:

        transcript_5p = random.choice(genes['transcript_id'].unique())
        transcript_3p = random.choice(genes['transcript_id'].unique())

        gene_5p = genes.loc[genes['transcript_id'] == transcript_5p, 'gene_id'].iloc[0]
        gene_3p = genes.loc[genes['transcript_id'] == transcript_3p, 'gene_id'].iloc[0]

        if gene_5p == gene_3p:
            continue

        exon_5p = random.choice(genes.loc[genes['transcript_id'] == transcript_5p, 'exon_number'].iloc[:-1].values)
        exon_3p = random.choice(genes.loc[genes['transcript_id'] == transcript_3p, 'exon_number'].iloc[1:].values)

        exons_5p = genes[(genes['transcript_id'] == transcript_5p) & (genes['exon_number'] <= exon_5p)]
        exons_3p = genes[(genes['transcript_id'] == transcript_3p) & (genes['exon_number'] >= exon_3p)]

        exons_5p = exons_5p.sort('exon_number')
        exons_3p = exons_3p.sort('exon_number')

        fusion_exons = pd.concat([exons_5p, exons_3p], ignore_index=True)

        fusion_sequence = create_gene_sequence(genome, fusion_exons)

        if len(fusion_sequence) < min_length:
            continue

        seq_id = 'fusion_{0}'.format(len(fusions_table))
        coverage = random.randint(min_coverage, max_coverage)

        fusion_exons['seq_id'] = seq_id
        fusion_exons['coverage'] = coverage

        reads1, reads2 = simulate_reads(seq_id, fusion_sequence, coverage, args.outdir, args.fragment_mean, args.fragment_stddev, args.read_length)

        fusions_table.append(fusion_exons)

        reads1_filenames.append(reads1)
        reads2_filenames.append(reads2)

    fusions_table = pd.concat(fusions_table, ignore_index=True)

    fusions_table_filename = os.path.join(args.outdir, 'fusions.tsv')

    fusions_table.to_csv(fusions_table_filename, sep='\t', index=False)

    wildtype_table = list()

    while len(wildtype_table) < args.num_wildtype:

        transcript = random.choice(genes['transcript_id'].unique())

        exons = genes[(genes['transcript_id'] == transcript)]

        exons = exons.sort('exon_number')

        wildtype_sequence = create_gene_sequence(genome, exons)

        if len(wildtype_sequence) < min_length:
            continue

        seq_id = 'wildtype_{0}'.format(len(wildtype_table))
        coverage = random.randint(min_coverage, max_coverage)

        reads1, reads2 = simulate_reads(seq_id, wildtype_sequence, coverage, args.outdir, args.fragment_mean, args.fragment_stddev, args.read_length)

        wildtype_table.append({'seq_id':seq_id, 'coverage':coverage})

        reads1_filenames.append(reads1)
        reads2_filenames.append(reads2)

    wildtype_table = pd.DataFrame(wildtype_table)

    wildtype_table_filename = os.path.join(args.outdir, 'wildtype.tsv')

    wildtype_table.to_csv(wildtype_table_filename, sep='\t', index=False)

    cat_files(args.reads1, reads1_filenames)
    cat_files(args.reads2, reads2_filenames)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='Ensembl gene models in gtf format')
    parser.add_argument('genome', help='Reference genome fasta')
    parser.add_argument('num_fusions', type=int, help='Number of fusions to simulate')
    parser.add_argument('num_wildtype', type=int, help='Number of wild type transcripts to simulate')
    parser.add_argument('outdir', help='Output directory')
    parser.add_argument('reads1', help='Fastq reads 1')
    parser.add_argument('reads2', help='Fastq reads 2')
    parser.add_argument('--seed', type=int, default=2014, help='RNG seed')
    parser.add_argument('--chromosomes', type=str, default=None, nargs='+', help='Restrict to specific chromosome')
    parser.add_argument('--fragment_mean', type=int, default=250., nargs='+', help='Mean length of fragments')
    parser.add_argument('--fragment_stddev', type=int, default=30., nargs='+', help='Standard deviation of fragment lengths')
    parser.add_argument('--read_length', type=int, default=75, nargs='+', help='Read lengths for simulated reads')
    args = parser.parse_args()

    simulate(args)


