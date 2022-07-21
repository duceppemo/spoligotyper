import os
from collections import OrderedDict
import subprocess


class SpoligoMethods(object):
    @staticmethod
    def count_spacers_seal(ref, output_folder, fastq_list, sample, cpu):
        # Stats file
        stats_file = output_folder + '/' + sample + '_stats.tsv'

        if len(fastq_list) == 1:
            cmd = ['seal.sh',
                   'in={}'.format(fastq_list[0]),
                   'ref={}'.format(ref),
                   'k=21',
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
                   'k=21',
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
    def print_report(output_folder, sample, spacer_count_dict, binary_code, octal_code, sb_code):
        report_file = output_folder + '/' + sample + '_spoligotyping.txt'

        spacer_counts = ':'.join([str(x) for x in list(spacer_count_dict.values())])
        header = 'Sample\tSpacerCount\tBinary\tOctal\tSpoligotype\n'
        results = '{}\t{}\t{}\t{}\t{}\n'.format(sample, spacer_counts, binary_code, octal_code, sb_code)

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