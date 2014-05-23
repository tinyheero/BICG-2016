# /usr/bin/env python
'''
This script will script will call somatic mutations from a VCF generate from a tumour/normal pair of BAMs 

Dependencies:
    PyVCF >= 0.6.3

@author Andrew Roth
@license GPLv3
'''
import csv
import vcf

def main(args):
    normal_column = args.normal_column
    
    if normal_column == 0:
        tumour_column = 1
    else:
        tumour_column = 0
     
    reader = vcf.Reader(open(args.vcf_file))

    for record in reader:
        normal_call = record.samples[normal_column]
        
        tumour_call = record.samples[tumour_column]
        
        # Skip germline positions
        if normal_call.is_variant:
            continue
        
        # If tumour is variant call somatic
        if not tumour_call.is_variant:
            continue
        
        # Skip low quality calls in tumour
        if tumour_call.data.GQ < args.min_genotype_quality:
            continue 
        
        print record.CHROM, record.POS

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()

    parser.add_argument('vcf_file',
                        help='''Path to VCF file containing calls from the normal genome. This will be used to identify
                        SNPs.''')
    
    parser.add_argument('--normal_column', type=int, default=0,
                        help='''0 based index of the column containing the normal sample. Choices are 0 or 1. The tumour
                        column will be assumed to be the other index.''')
    
    parser.add_argument('--min_genotype_quality', type=int, default=99)
    
    args = parser.parse_args()
    
    main(args)
