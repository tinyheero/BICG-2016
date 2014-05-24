import os
import argparse
import gzip
import random
import string
import collections



def read_gtf(line):

    fields = line.rstrip().split('\t')
        
    keyvals = fields[-1].replace('; ', ';').rstrip(';').split(';')
    keyvals = [a.split('=', 1) for a in keyvals]
    keyvals = [(a, b.strip('"').rstrip('"')) for a, b in keyvals]
    keyvals = dict(keyvals)
    
    keyvals['chromosome'] = fields[0]
    keyvals['start'] = int(fields[3])
    keyvals['end'] = int(fields[4])
    keyvals['strand'] = fields[6]
    keyvals['type'] = fields[2]

    return keyvals


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input_gtf', help='Transcript alignments in gtf format')
    parser.add_argument('input_fasta', help='Transcript sequences in fasta format')
    parser.add_argument('output_fasta', help='Output fusion transcript sequences in fasta format')
    args = parser.parse_args()

    gene_alignments = collections.defaultdict(set)

    with open(args.input_gtf, 'r') as input_gtf:
        
        for line in input_gtf:

            if line.startswith('#'):
                continue

            gtf_info = read_gtf(line)

            if gtf_info['type'] != 'gene':
                continue

            gene_alignments[gtf_info['Name']].add(gtf_info['ID'])

    sequences = collections.defaultdict(list)

    with open(args.input_fasta, 'r') as input_fasta_file:

        name = None
        
        for line in input_fasta_file:

            if line.startswith('>'):
                name = line[1:].rstrip()
            elif name is not None:
                sequences[name].append(line)

    with open(args.output_fasta, 'w') as output_fasta_file:

        for name, sequence_lines in sequences.iteritems():

            short_name = name.split()[0]

            if len(gene_alignments[short_name]) > 1:

                output_fasta_file.write('>{0}\n'.format(name))
                for line in sequence_lines:
                    output_fasta_file.write(line)

