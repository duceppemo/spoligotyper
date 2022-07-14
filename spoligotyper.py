#!/usr/bin/env python


__author__ = 'duceppemo'
__version__ = 'v0.1'


import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from collections import OrderedDict
import subprocess


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
        self.spoligo_fasta = './dependencies/spoligo_spacers.fasta'

        # Run
        self.run()

    def run(self):
        print('')
        # Count spacers
        count_summary_dict = OrderedDict(SpoligoMethods.count_spacer_occurence(self.fastq_list, self.spoligo_fasta))

        # Convert to binary
        binary_list = SpoligoMethods.spoligo(count_summary_dict)

        # Convert binary to octal
        octal_list = SpoligoMethods.binary_to_octal(binary_list)

        # Convert return serotype
        SpoligoMethods.spoligo(count_summary_dict)


class SpoligoMethods(object):
    @staticmethod
    def count_spacer_occurence(fastq_list, spoligo_fasta):
        if len(fastq_list) == 1:
            cmd = ['bbduk.sh',
                   'in={}'.format(fastq_list[0]),
                   'ref={}'.format(spoligo_fasta),
                   'k=21', 'rcomp=t',
                   'edist=1',  # up to 1 missmatch
                   'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
                   'stats=stats.txt',
                   'ow=t']
        else:
            cmd = ['bbduk.sh',
                   'in={}'.format(fastq_list[0]),
                   'in2={}'.format(fastq_list[1]),
                   'ref={}'.format(spoligo_fasta),
                   'k=21', 'rcomp=t',
                   'edist=1',  # up to 1 missmatch
                   'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
                   'stats=stats.txt',
                   'ow=t']

        subprocess.run(cmd)

        # Parse stats file
        """
        #File /home/bioinfo/analyses/vsnp3_test_spoligo/SRR16058435_R1.fastq.gz
        #Total  822714
        #Matched    799 0.09712%
        #Name   Reads   ReadsPct
        spacer25    62  0.00754%
        spacer02    48  0.00583%
        """

        count_summary = dict()
        with open('stats.txt', 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue

                if line.startswith('#'):
                    continue
                else:
                    field_list = line.split('\t')
                    count_summary[field_list[0]] = int(field_list[1])

        # Fill any spacers with zero counts
        # Parse spacer fasta file
        spacer_list = list()
        with open(spoligo_fasta, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    spacer_list.append(line[1:])  # Drop the leading ">"

        for spacer_id in spacer_list:
            if spacer_id not in count_summary:
                count_summary[spacer_id] = 0

        return count_summary

    @staticmethod
    def binary_to_octal(binary):
        # binary_len = len(binary)
        i = 0
        ie = 1
        octal = ""
        while ie < 43:
            ie = i + 3
            region = binary[i:ie]
            region_len = len(region)
            i += 3
            if int(region[0]) == 1:
                if region_len < 2:  # for the lone spacer 43.  When present needs to be 1 not 4.
                    oct = 1
                else:
                    oct = 4
            else:
                oct = 0
            try:
                if int(region[1]) == 1:
                    oct += 2
                if int(region[2]) == 1:
                    oct += 1
            except IndexError:
                pass
            octal = octal + str(oct)
        return (octal)

    @staticmethod
    def spoligo(count_summary_dict):
        spoligo_binary_dictionary = {}
        call_cut_off = 4
        for k, v in count_summary_dict.items():
            if v > call_cut_off:
                spoligo_binary_dictionary.update({k: 1})
            else:
                spoligo_binary_dictionary.update({k: 0})
        spoligo_binary_dictionary = OrderedDict(sorted(spoligo_binary_dictionary.items()))
        spoligo_binary_list = []
        for v in spoligo_binary_dictionary.values():
            spoligo_binary_list.append(v)
        binary_list = ''.join(str(e) for e in spoligo_binary_list)  # sample_binary correct

        return binary_list

    @staticmethod
    def binary_to_sbcode(binary_list, count_summary_dict, spoligo_db):
        with open(spoligo_db) as f:  # put into dictionary or list
            for line in f:
                line = line.rstrip()
                sbcode = line.split()[1]
                db_binarycode = line.split()[2]
                if binary_list == db_binarycode:
                    found = True
                    sbcode = sbcode
        if not found:
            if binary_list == '0000000000000000000000000000000000000000000':
                sbcode = "spoligo not found, binary all zeros, see spoligo file"
            else:
                sbcode = "Not Found"
        count_summary_list = []
        for spacer, count in count_summary_dict.items():
            count_summary_list.append(count)


if __name__ == "__main__":
    max_cpu = cpu_count()

    parser = ArgumentParser(description='In silico spoligotyping from WGS data for Mycobacterium bovis.')
    parser.add_argument('-r1', metavar='/path/to/reference_organelle/genome.fasta',
                        required=True,
                        help='R1 fastq file from paired-end or single end sequencing data. Mandatory.')
    parser.add_argument('-r2', '--input', metavar='/path/to/input/folder/',
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