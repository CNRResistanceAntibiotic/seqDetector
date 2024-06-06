#!/usr/bin/python3
import os
import re
import argparse
from collections import OrderedDict
from subprocess import Popen, PIPE, STDOUT

import numpy as np
import pandas as pd


def write_site_file(position, output_dir):
    pattern = re.compile('([0-9a-zA-Z_-]+):*([0-9]*)-*([0-9]*)')
    match = pattern.match(position)
    if match:
        txt = ''
        site_file = os.path.join(os.path.dirname(output_dir), 'site_file.txt')
        for item in match.groups():
            if item != '':
                txt = txt + f'\t{item}'
        with open(site_file, 'w') as out_f:
            out_f.write(txt[1:])
        return site_file
    else:
        print('\nPositions not identified\n')
        print('Position formats:\n'
              '\t<contig or chromosome name>\n'
              '\t<contig or chromosome name>:<position in chromome or contig>\n'
              '\t<contig or chromosome name>:<start position in chromome or contig>:<end position in chromome or contig>')
        exit(1)


def bam_count(bam_file, fasta_ref, output_dir, q=0, b=0, feature_name='', site_file='', force=False):
    sample = os.path.basename(bam_file).split(".")[0]
    if feature_name == '' and site_file == '':
        out_file = os.path.join(output_dir, f'{sample}_whole_genome_raw.csv')
        cmd = f'$(which bam-readcount) -w 0 -q {q} -b {b} -i -f {fasta_ref} {bam_file} > {out_file}'
    elif feature_name != '' and site_file == '':
        out_file = os.path.join(output_dir, f'{sample}_{feature_name}_raw.csv')
        cmd = f'$(which bam-readcount) -w 0 -q {q} -b {b} -i -f {fasta_ref} {bam_file} > {out_file}'
    else:
        out_file = os.path.join(output_dir, f'{sample}_{feature_name}_raw.csv')
        cmd = f'$(which bam-readcount) -w 0 -q {q} -b {b} -i -l {site_file} -f {fasta_ref} {bam_file} > {out_file}'
    if not os.path.exists(out_file) or force:
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT).stdout.read()
        log_info = f"Command line executed: {cmd}\n\n\n{process.decode('utf-8')}"
        print(log_info)
    else:
        print(f'\nRead count file {out_file} already done.\n')
    return out_file


def calculate_stats(data):
    result = OrderedDict()
    result['perc>=30'] = round(100 * len([x for x in data if x >= 30]) / float(len(data)), 2)
    result['perc>=20'] = round(100 * len([x for x in data if x >= 20]) / float(len(data)), 2)
    result['perc>=10'] = round(100 * len([x for x in data if x >= 10]) / float(len(data)), 2)
    result['mean'] = round(sum(data) / float(len(data)), 2)
    result['50_perc'] = round(np.percentile(np.array(data), 50), 2)
    result['25_perc'] = round(np.percentile(np.array(data), 25), 2)
    result['75_perc'] = round(np.percentile(np.array(data), 75), 2)
    result['max'] = max(data)
    result['min'] = min(data)
    return result


def bam_count_stats(bam_count_file, feature_name, header, output_dir, bam_file):

    sample = os.path.basename(bam_file).split(".")[0]

    if feature_name == '':
        out_file = os.path.join(output_dir, f'{sample}_whole_genome_stats')
    else:
        out_file = os.path.join(output_dir, f'{sample}_{feature_name}_stats')
    with open(bam_count_file) as count_f:
        ctgs = []
        result_stat = []
        result_data = []
        start = 0
        end = 0
        base_nb = 0

        for line in count_f:
            line = line.strip().split('\t')
            ctg = line[0]

            if ctg not in ctgs:
                if len(ctgs) > 0:
                    d = OrderedDict([('ID', ctgs[-1]), ('start', start), ('end', end), ('size', base_nb)])
                    result_data.append(d)
                    s = []
                    # for key, value in calculate_stats(ctg_depth).iteritems():
                    #    s.append(('All_depth_{0}'.format(key), value))
                    for key, value in calculate_stats(ctg_ref_depth).items():
                        s.append((f'Ref_depth_{key}', value))
                    for key, value in calculate_stats(ctg_ref_qual).items():
                        s.append((f'Ref_quali_{key}', value))
                    s = OrderedDict(s)
                    result_stat.append(s)
                ctg_ref_depth, ctg_ref_qual = [], []
                # ctg_ref_depth = []
                ctgs.append(ctg)
                start = -1
                end = 1
                base_nb = 0

            pos = int(line[1])
            if start == -1:
                start = pos
            if pos >= end:
                end = pos
            base_nb += 1

            ref = line[2]
            # ctg_depth.append(int(line[3]))

            for item in line[4:]:
                data = dict(zip(header.split(':'), item.split(':')))
                base = data['base']
                if base == ref:
                    ctg_ref_depth.append(int(data['count']))
                    ctg_ref_qual.append(round(float(data['avg_base_quality']), 2))

        # print(start)
        # print(end)
        # print(base_nb)
        # print(ctgs)
        if ctgs:
            d = OrderedDict([('ID', ctgs[-1]), ('start', start), ('end', end), ('size', base_nb)])
            result_data.append(d)
            s = []
            # for key, value in calculate_stats(ctg_depth).iteritems():
            #    s.append(('All_depth_{0}'.format(key), value))
            for key, value in calculate_stats(ctg_ref_depth).items():
                s.append((f'Ref_depth_{key}', value))
            for key, value in calculate_stats(ctg_ref_qual).items():
                s.append((f'Ref_quali_{key}', value))
            s = OrderedDict(s)
            result_stat.append(s)

            df_data = pd.DataFrame.from_dict(result_data)
            df_stat = pd.DataFrame.from_dict(result_stat)
            weighted_mean = ((df_stat.multiply(df_data['size'], axis=0).sum()).div(df_data['size'].sum(), axis=0)).round(2)
            df_stat = df_stat.append(weighted_mean, ignore_index=True)
            df_stat = df_stat.round(2)
            df_data = df_data.append([{'ID': 'Overall', 'size': df_data['size'].sum(), 'end': '-', 'start': '-'}],
                                     ignore_index=True)
            df = pd.concat([df_data, df_stat], axis=1)
            df.sort_values('size', inplace=True)
            df.to_csv(f'{out_file}.csv', sep='\t', header=True, index=True)
            df.to_html(f'{out_file}.html', index=True)
            return f'{out_file}.csv'
        else:
            return ''


def bam_count_extract(bam_count_file, feature_name, header, output_dir, bam_file):

    sample = os.path.basename(bam_file).split(".")[0]
    if feature_name == '':
        out_file = os.path.join(output_dir, sample + f'_{"whole_genome"}_count')
    else:
        out_file = os.path.join(output_dir, sample + f'_{feature_name}_count')

    result_data = []
    with open(bam_count_file) as count_f:
        for line in count_f:
            line = line.strip().split('\t')
            ctg = line[0]
            pos = int(line[1])
            ref = line[2]
            depth = int(line[3])
            data = [('ID', ctg), ('position', pos), ('reference', ref), ('total_depth', depth)]
            for item in line[4:]:
                d = dict(zip(header.split(':'), item.split(':')))
                base = d['base']
                if base != '=':
                    data.append((f'{base}_depth', int(d['count'])))
                    data.append((f'{base}_quality', round(float(d['avg_base_quality']), 2)))
            data = OrderedDict(data)
            result_data.append(data)

        df = pd.DataFrame.from_dict(result_data)
        # print(df)
        df.to_csv(f'{out_file}.csv', sep='\t', header=True, index=True)
        df.to_html(f'{out_file}.html', index=True)
    if result_data:
        return f'{out_file}.csv'
    else:
        return ''


def pre_main(args):
    bam_file = args.bam_file
    fasta_ref = args.fas_file
    force = args.Force
    mapping_qual = int(args.mapping_qual)
    base_qual = int(args.base_qual)
    stats = args.stat
    data = args.data
    position = args.position
    feature_name = args.name
    output_dir = args.output
    main(bam_file, fasta_ref, position, output_dir, feature_name, force, mapping_qual, base_qual, stats, data)


def main(bam_file, fasta_ref, position, output_dir, feature_name, force=False, mapping_qual=0, base_qual=0, stats=False,
         data=False):

    mapping_qual = int(mapping_qual)
    base_qual = int(base_qual)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    while feature_name == '' and position != '':
        # feature_name = raw_input('Enter a feature name for the position and validate by enter:\n')
        print("YOU NEED TO SELECT A FEATURE")
        exit(1)
    if feature_name == '':
        position = ''
    if position != '':
        site_file = write_site_file(position, output_dir)
    bam_count_file = bam_count(bam_file, fasta_ref, output_dir, mapping_qual, base_qual, feature_name, site_file, force)
    if os.path.exists(site_file):
        os.remove(site_file)

    header = 'base:count:avg_mapping_quality:avg_base_quality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:' \
             'avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:' \
             'avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end'

    out_file_list = [bam_count_file]
    if stats:
        stat_file = bam_count_stats(bam_count_file, feature_name, header, output_dir, bam_file)
        if stat_file:
            out_file_list.append(stat_file)
        else:
            print("\n Error ! No stats file produced \n")
    if data:
        data_file = bam_count_extract(bam_count_file, feature_name, header, output_dir, bam_file)
        if data_file:
            out_file_list.append(data_file)
        else:
            print("\n Error ! No data file produced \n")
    return out_file_list


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='tsv2fasta - Version ' + version())
    parser.add_argument('-b', '--bamFile', dest="bam_file", default='data/CNR2470.bam', help='Bam input file')
    parser.add_argument('-f', '--fasFile', dest="fas_file", default='data/CNR2470.fasta', help='Fasta reference file')
    parser.add_argument('-n', '--nameFeature', dest="feature_name", default='data/CNR2470.bam', help='Name of data')
    parser.add_argument('-p', '--position', dest="position", default='ctg_7:30-35',
                        help='Selection option: <ctg name>:<start base>-<end base>')
    parser.add_argument('-mq', '--mapQual', dest="mapping_qual", default='0', help='Minimum mapping quality [0]')
    parser.add_argument('-bq', '--basQual', dest="base_qual", default='0', help="Minimum base quality [0]")
    parser.add_argument('-st', '--stats', dest="stats", action="store_true", default=False,
                        help="Produce depth and quality stats")
    parser.add_argument('-dt', '--data', dest="data", action="store_true", default=False,
                        help="Extract depth and quality data")
    parser.add_argument('-F', '--Force', dest="Force", action="store_true", default=False,
                        help="Overwrite the detected files")
    parser.add_argument('-o', '--out', dest="output", default="",
                        help="Output directory")

    parser.add_argument('-V', '--version', action='version', version='bam2data-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)
