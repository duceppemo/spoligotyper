#!/usr/bin/env python


__author__ = 'duceppemo'
__version__ = 'v0.1'


import os
import mmap
from argparse import ArgumentParser
from multiprocessing import cpu_count
from collections import OrderedDict
import subprocess
from glob import glob
from time import time


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
        self.parallel = args.parallel

        # Data
        self.spoligo_db = './dependencies/spoligotype_db.txt'
        self.spacer_folder = './dependencies/spacers'

        # Run
        self.run()

    def run(self):
        start_time = time()

        sample = os.path.basename(self.fastq_list[0]).split('_')[0]
        print('Spoligotyping {}'.format(sample))

        # Count spacers
        # spacer_count_dict = SpoligoMethods.run_bbduk_parallel(self.fastq_list, self.spacer_folder,
        #                                                       self.cpu, self.parallel)

        spacer_count_dict = SpoligoMethods.count_spacers_bbduk(self.output, self.fastq_list, self.spacer_folder, sample)

        # Original count method
        # count_summary = SpoligoMethods.count_spacers_regex(self.fastq_list)

        # Convert to binary
        binary_code = SpoligoMethods.spoligo(spacer_count_dict)

        # Convert binary to octal
        octal_code = SpoligoMethods.binary_to_octal(binary_code)

        # Convert binary to spoligotype
        sb_code = SpoligoMethods.binary_to_sbcode(binary_code, self.spoligo_db)

        # TODO: Get fastq file info to add to report (number of reads, average read size, GC%, etc.)

        # Print report
        report_file = self.output + '/' + sample + '_spoligotyping.txt'

        spacer_counts = ':'.join([str(x) for x in list(spacer_count_dict.values())])
        header = 'Sample\tSpacerCount\tBinary\tOctal\tSpoligotype\n'
        results = '{}\t{}\t{}\t{}\t{}\n'.format(sample, spacer_counts, binary_code, octal_code, sb_code)

        with open(report_file, 'w') as f:
            f.write(header)
            f.write(results)

        print('\nElapsed time: {}'.format(SpoligoMethods.elapsed_time(time() - start_time)))


class SpoligoMethods(object):
    @staticmethod
    def run_bbduk(mmap_obj, spacer_file, cpu):
        spacer = os.path.basename(spacer_file).split('.')[0]
        print('Checking for {}'.format(spacer))

        # stdout = ''
        # if len(fastq_list) == 1:
        #     cmd = ['bbduk.sh',
        #            # 'in={}'.format(fastq_list[0]),
        #            'in=stdin.fq.gz',
        #            'int=f',
        #            'ref={}'.format(spacer_file),
        #            'k=25', 'rcomp=t',
        #            'hdist=1',  # up to 1 missmatch
        #            'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
        #            'ow=t',
        #            'threads={}'.format(cpu)]
        #     p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #     sdtout = p.communicate(input=mmap_obj)[0]  # communicate returns a tuple (stdout_data, stderr_data)
        #     test =1
        # else:
        #     cmd = ['bbduk.sh',
        #            'in={}'.format(fastq_list[0]),
        #            'in2={}'.format(fastq_list[1]),
        #            'ref={}'.format(spacer_file),
        #            'k=25', 'rcomp=t',
        #            'hdist=1',  # up to 1 missmatch
        #            'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
        #            'ow=t',
        #            'threads={}'.format(cpu)]
        #     # https://fredrikaverpil.github.io/2013/10/11/catching-string-from-stdout-with-python/
        #     p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #     # p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #     sdtout = p.communicate()[0]

        cmd = ['bbduk.sh',
               'in=stdin.fq.gz',
               'int=f',
               'ref={}'.format(spacer_file),
               'k=25', 'rcomp=t',
               'hdist=1',  # up to 1 missmatch
               'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
               'ow=t',
               'threads={}'.format(cpu)]
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate(input=mmap_obj)  # communicate returns a tuple (stdout_data, stderr_data)

        spacer_dict = dict()
        for line in stdout.decode('ascii').splitlines():
            if line.startswith('Contaminants'):
                count = int(line.split()[1])
                spacer_dict[spacer] = count
                return spacer_dict

    @staticmethod
    def run_bbduk_parallel(fastq_list, spacer_folder, cpu, parallel):
        spacer_file_list = glob(spacer_folder + '/*.fasta')
        spacer_dict = dict()

        # with futures.ThreadPoolExecutor(max_workers=2) as executor:
        #     with open(fastq_list[0], 'rb') as file_obj:
        #         with mmap.mmap(file_obj.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
        #             args = ((mmap_obj, spacer_file, "8") for spacer_file in spacer_file_list)
        #             for results in executor.map(lambda x: SpoligoMethods.run_bbduk(*x), args):
        #                 spacer_dict.update(results)

        for spacer_file in spacer_file_list:
            with open(fastq_list[0], 'rb') as file_obj:
                with mmap.mmap(file_obj.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
                    spacer_dict.update(SpoligoMethods.run_bbduk(mmap_obj, spacer_file, int(cpu / parallel)))

        spacer_dict = OrderedDict(sorted(spacer_dict.items()))
        return spacer_dict

    @staticmethod
    def count_spacers_bbduk(output_folder, fastq_list, spacer_folder, sample):
        # List all spacer files
        spacer_file_list = glob(spacer_folder + '/*.fasta')

        # Stats file
        stats_file = output_folder + '/' + sample + '_stats.tsv'

        # To store the spacer counts
        spacer_count_dict = dict()

        # Loop spacer files one after the other
        for spacer_file in spacer_file_list:
            spacer = os.path.basename(spacer_file).split('.')[0]
            print('Checking for {}'.format(spacer))

            if len(fastq_list) == 1:
                cmd = ['bbduk.sh',
                       'in={}'.format(fastq_list[0]),
                       'ref={}'.format(spacer_file),
                       'k=25', 'rcomp=t',
                       'hdist=1',  # up to 1 missmatch
                       'maskmiddle=f']  # Do not treat the middle base of a kmer as a wildcard
                       # 'ow=t',
                       # 'stats={}'.format(stats_file), 'nzo=t']
            else:
                cmd = ['bbduk.sh',
                       'in={}'.format(fastq_list[0]),
                       'in2={}'.format(fastq_list[1]),
                       'ref={}'.format(spacer_file),
                       'k=25', 'rcomp=t',
                       'hdist=1',  # up to 1 missmatch
                       'maskmiddle=f']  # Do not treat the middle base of a kmer as a wildcard
                       # 'ow=t',
                       # 'stats={}'.format(stats_file), 'nzo=t']

            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout = p.communicate()[0]

            for line in stdout.decode('ascii').splitlines():
                if line.startswith('Contaminants'):
                    count = int(line.split()[1])
                    spacer_count_dict[spacer] = count

        # Parse stats file
        """
        #File /home/bioinfo/analyses/vsnp3_test_spoligo/SRR16058435_R1.fastq.gz
        #Total  822714
        #Matched    799 0.09712%
        #Name   Reads   ReadsPct
        spacer25    62  0.00754%
        spacer02    48  0.00583%
        """

        # spacer_count_dict = dict()
        # with open(stats_file, 'r') as f:
        #     for line in f:
        #         line = line.rstrip()
        #         if not line:
        #             continue
        #
        #         if line.startswith('#'):
        #             continue
        #         else:
        #             field_list = line.split('\t')
        #             spacer_count_dict[field_list[0]] = int(field_list[1])

        # Fill any spacers with zero counts
        # Parse spacer fasta file
        # spacer_list = list()
        # with open(spoligo_fasta, 'r') as f:
        #     for line in f:
        #         line = line.rstrip()
        #         if not line:
        #             continue
        #         if line.startswith('>'):
        #             spacer_list.append(line[1:])  # Drop the leading ">"
        #
        # for spacer_id in spacer_list:
        #     if spacer_id not in spacer_count_dict:
        #         spacer_count_dict[spacer_id] = 0

        # Order dictionary by keys
        spacer_count_dict = OrderedDict(sorted(spacer_count_dict.items()))

        return spacer_count_dict

    @staticmethod
    def binary_to_octal(binary_rep):
        i = 0
        ie = 1
        octal_rep = ""
        while ie < 43:
            ie = i + 3
            region = binary_rep[i:ie]
            region_len = len(region)
            i += 3
            if int(region[0]) == 1:
                if region_len < 2:  # for the lone spacer 43.  When present needs to be 1 not 4.
                    octal = 1
                else:
                    octal = 4
            else:
                octal = 0
            try:
                if int(region[1]) == 1:
                    octal += 2
                if int(region[2]) == 1:
                    octal += 1
            except IndexError:
                pass
            octal_rep = octal_rep + str(octal)
        return octal_rep

    @staticmethod
    def spoligo(count_summary_dict):
        binary_rep = ''
        for k, v in count_summary_dict.items():
            if v > 4:  # call cut off
                binary_rep += '1'
            else:
                binary_rep += '0'
        return binary_rep

    @staticmethod
    def binary_to_sbcode(binary_code, spoligo_db):
        # Parse spolifo_db file into dictionary
        binary_to_spoligo_dict = dict()
        with open(spoligo_db) as f:  # put into dictionary or list
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                spoligo, binary_rep = line.split(" ")[1:]
                binary_to_spoligo_dict[binary_rep] = spoligo

        if binary_code in binary_to_spoligo_dict:
            return binary_to_spoligo_dict[binary_code]
        else:
            return "Spoligo not found"

    @staticmethod
    def elapsed_time(seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formatted time string
        """
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(round(value)), name) for name, value in periods if value)
        return time_string


if __name__ == "__main__":
    max_cpu = cpu_count()

    parser = ArgumentParser(description='In silico spoligotyping from WGS data for Mycobacterium bovis.')
    parser.add_argument('-r1', metavar='/path/to/reference_organelle/genome.fasta',
                        required=True,
                        help='R1 fastq file from paired-end or single end sequencing data. Mandatory.')
    parser.add_argument('-r2', metavar='/path/to/input/folder/',
                        required=False,
                        help='R2 fastq file from paired-end. Optional.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder/',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel. Default is 2. Optional.')
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.abspath(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Spoligo(arguments)
