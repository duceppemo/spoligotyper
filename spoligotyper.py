#!/usr/bin/env python


__author__ = 'duceppemo'
__version__ = 'v0.1'


import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from time import time
from spoligotyper_methods import SpoligoMethods


class Spoligo(object):
    def __init__(self, args):
        # I/O
        if args.r2:
            self.fastq_list = [args.r1, args.r2]
        else:
            self.fastq_list = [args.r1]
        self.output = args.output

        # Performance
        self.cpu = args.threads

        # Data
        self.spoligo_db = './dependencies/spoligotype_db.txt'
        self.spacer_folder = './dependencies/spacers'
        self.spoligo_spacers = './dependencies/spoligo_spacers.fasta'

        # Run
        self.run()

    def run(self):
        start_time = time()

        sample = os.path.basename(self.fastq_list[0]).split('_')[0]
        print('Spoligotyping {}'.format(sample))

        # Count spacers
        spacer_count_dict = SpoligoMethods.count_spacers_seal(self.spoligo_spacers, self.output, self.fastq_list,
                                                              sample, self.cpu)

        # Convert to binary
        binary_code = SpoligoMethods.spoligo(spacer_count_dict)

        # Convert binary to octal
        octal_code = SpoligoMethods.binary_to_octal(binary_code)

        # Convert binary to spoligotype
        sb_code = SpoligoMethods.binary_to_sbcode(binary_code, self.spoligo_db)

        # Print report on scree and to file
        SpoligoMethods.print_report(self.output, sample, spacer_count_dict, binary_code, octal_code, sb_code)

        # Print elapsed time
        print('Elapsed time: {}'.format(SpoligoMethods.elapsed_time(time() - start_time)))


if __name__ == "__main__":
    max_cpu = cpu_count()

    parser = ArgumentParser(description='In silico spoligotyping from WGS data for Mycobacterium bovis.')
    parser.add_argument('-r1', metavar='/path/to/R1/fastq/or/fasta',
                        required=True,
                        help='R1 fastq file from paired-end or single end sequencing data. Can be a fasta file too. '
                             'Gzipped or not. Mandatory.')
    parser.add_argument('-r2', metavar='/path/to/R2/fastq',
                        required=False,
                        help='R2 fastq file from paired-end. Optional.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder/',
                        required=True,
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
