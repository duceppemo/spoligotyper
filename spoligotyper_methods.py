import os
from collections import OrderedDict
import subprocess
import pathlib


class SpoligoMethods(object):
    @staticmethod
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def count_spacers_seal(ref, output_folder, fastq_list, sample, cpu):
        # Stats file
        stats_file = output_folder + '/' + sample + '_stats.tsv'

        # Create output folder
        SpoligoMethods.make_folder(output_folder)

        if len(fastq_list) == 1:
            cmd = ['seal.sh',
                   'in={}'.format(fastq_list[0]),
                   'ref={}'.format(ref),
                   'k=25',
                   'rcomp=t',
                   'hdist=1',  # up to 1 missmatch
                   'clearzone=999999',
                   'ambiguous=all',
                   'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
                   'ow=t',
                   'stats={}'.format(stats_file),
                   'nzo=f',
                   'threads={}'.format(cpu)]
        else:
            cmd = ['seal.sh',
                   'in={}'.format(fastq_list[0]),
                   'in2={}'.format(fastq_list[1]),
                   'ref={}'.format(ref),
                   'k=25',
                   'rcomp=t',
                   'hdist=1',  # up to 1 missmatch
                   'clearzone=999999',
                   'ambiguous=all',
                   'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
                   'ow=t',
                   'stats={}'.format(stats_file),
                   'nzo=f',
                   'threads={}'.format(cpu)]

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Parse stats file
        """
        #File /home/bioinfo/analyses/vsnp3_test_spoligo/SRR16058435_R1.fastq.gz
        #Total  822714
        #Matched    799 0.09712%
        #Name   Reads   ReadsPct
        spacer25    62  0.00754%
        spacer02    48  0.00583%
        """

        spacer_count_dict = dict()
        with open(stats_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue

                if line.startswith('#'):
                    continue
                else:
                    field_list = line.split('\t')
                    spacer_count_dict[field_list[0]] = int(field_list[1])

        # Remove temporary stat file
        os.remove(stats_file)

        return OrderedDict(sorted(spacer_count_dict.items()))  # Order dictionary by keys

    @staticmethod
    def count_binary(count_summary_dict, min_count):
        binary_rep = ''
        for k, v in count_summary_dict.items():
            if v >= min_count:  # call cut off
                binary_rep += '1'
            else:
                binary_rep += '0'
        return binary_rep

    @staticmethod
    def binary_to_octal(binary_rep):
        oct_rep = ''
        for octal_bloc in [binary_rep[i:i + 3] for i in range(0, len(binary_rep), 3)]:
            oct_rep += str(oct(int(octal_bloc, 2)))[-1]
        return oct_rep

    @staticmethod
    def binary_to_hexa(binary_rep):
        hex_bloc_list = [(0, 7), (7, 14), (14, 21), (21, 28), (28, 36), (36, 43)]
        hex_rep = ''
        for hex_bloc in [binary_rep[hex_bloc_list[i][0]:hex_bloc_list[i][1]] for i in range(0, len(hex_bloc_list))]:
            hex_rep += str(hex(int(hex_bloc, 2)))[-2:].upper().replace('X', '0') + '-'
        return hex_rep[:-1]

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
    def print_report(output_folder, sample, spacer_count_dict, binary_code, octal_code, hex_code, sb_code):
        report_file = output_folder + '/' + sample + '_spoligotyping.txt'

        spacer_counts = ':'.join([str(x) for x in list(spacer_count_dict.values())])
        header = 'Sample\tSpacerCount\tBinary\tOctal\tHexadecimal\tSpoligotype\n'
        results = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(sample, spacer_counts, binary_code, octal_code, hex_code, sb_code)

        with open(report_file, 'w') as f:
            f.write(header)
            f.write(results)

        print('\n{}{}\nResults saved in {}\n'.format(header, results, report_file))

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