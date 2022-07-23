#!/usr/bin/env python


__author__ = 'duceppemo'
__version__ = 'v0.1'


import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from time import time
from spoligotyper_methods import SpoligoMethods
import pkg_resources


class Spoligo(object):
    def __init__(self, args):
        # I/O
        if args.r2:
            self.fastq_list = [args.r1, args.r2]
        else:
            self.fastq_list = [args.r1]
        self.output = args.output
        self.min_count = args.min_count

        # Performance
        self.cpu = args.threads

        # Data - Need to create a package for that
        self.spoligo_db = pkg_resources.resource_filename('dep', 'spoligotype_db.txt')
        self.spoligo_spacers = pkg_resources.resource_filename('dep', 'spoligo_spacers.fasta')

        # Run
        self.run()

    def run(self):
        start_time = time()

        # Get sample name
        sample = '.'.join(os.path.basename(self.fastq_list[0]).split('.')[:-1])  # basename minus what is after last "."
        if self.fastq_list[0].endswith('gz'):  # Need to drop what after the last 2 "."
            sample = '.'.join(sample.split('.')[:-1])
        if '_R1' in sample:
            sample = '_'.join(sample.split("_")[:-1])

        print('Spoligotyping {}'.format(sample))

        # Count spacers
        spacer_count_dict = SpoligoMethods.count_spacers_seal(self.spoligo_spacers, self.output, self.fastq_list,
                                                              sample, self.cpu)

        # Convert to binary
        binary_code = SpoligoMethods.count_binary(spacer_count_dict, self.min_count)

        # Convert binary to octal
        octal_code = SpoligoMethods.binary_to_octal(binary_code)

        # Convert binary to hexadecimal
        hex_code = SpoligoMethods.binary_to_hexa(binary_code)

        # Convert binary to spoligotype
        sb_code = SpoligoMethods.binary_to_sbcode(binary_code, self.spoligo_db)

        # Print report on scree and to file
        SpoligoMethods.print_report(self.output, sample, spacer_count_dict, binary_code, octal_code, hex_code, sb_code)

        # Print elapsed time
        print('Elapsed time: {}'.format(SpoligoMethods.elapsed_time(time() - start_time)))


if __name__ == "__main__":
    max_cpu = cpu_count()

    parser = ArgumentParser(description='In silico spoligotyping from WGS data for Mycobacterium bovis.')
    parser.add_argument('-r1', metavar='/path/to/sample_R1.[fastq|fasta]',
                        required=True, type=str,
                        help='R1 fastq file from paired-end or single end sequencing data. Can be a fasta file. '
                             'Gzipped or not. Mandatory.')
    parser.add_argument('-r2', metavar='/path/to/R2/fastq',
                        required=False, type=str,
                        help='R2 fastq file from paired-end. Optional.')
    parser.add_argument('-m', '--min-count', metavar='5',
                        required=False, default=4, type=int,
                        help='Minimum count to consider a spacer to be present. Set to 1 if "-r1" is a fasta file. '
                             'Default 5. Optional.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder/',
                        required=True, type=str,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional'.format(max_cpu))
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.abspath(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Spoligo(arguments)
