#!/usr/bin/python3
import argparse
import glob
import itertools
import operator
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from seqdetector.blast_functions import make_blastn_database, run_blastn, load_blastn_result, \
    dna_extract_quality_and_depth, cds_extract_quality_and_depth
from seqdetector.misc_functions import load_fasta, make_dmnd_database, dna_global_alignemnt, cds_global_alignment,\
    load_sample_name


def run_diam(dmnd_db, query_file, pass_pid=70, pass_pcv=70, threads=8, force=True, out_file=""):
    if not force and os.path.exists(out_file):
        print(f'\nResult file {out_file} already exists')
    else:
        fmt = '6 qseqid qframe sseqid slen qstart qend sstart send length pident nident ppos positive mismatch ' \
              'gapopen gaps qseq sseq full_sseq'
        cmd = f"$(which diamond) blastx --more-sensitive -k 0 -p {threads} -d {dmnd_db} -q {query_file}" \
              f" --id {pass_pid} --subject-cover {pass_pcv} -f {fmt} -o {out_file} --masking no"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
        print(f"\n{cmd}\n{process.decode('utf-8')}")
    return out_file


def load_dmnd_result(result_file, target_file):
    print("\n***** Start loading Diamond Result *****")
    target_dic = load_fasta(target_file)
    dmnd_results = []
    header = ['qid', 'strand', 'tid', 'tlen', 'qstart', 'qend', 'tstart', 'tend', 'alen', 'pid', 'nid', 'ppos', 'npos',
              'mismatch', 'gapopen', 'gap', 'qseq', 'tseq', 'fulltseq']

    with open(result_file) as inf:
        for line in inf:
            data = dict(zip(header, line.strip().split('\t')))
            # qid = data['qid']
            tid = data['tid']
            tdes = {}
            for item in target_dic[tid].description.split(';'):
                key, value = item.split(':')
                if key == 'prot_snp':
                    key = 'known_prot_snp'
                tdes[key] = value
            item_dic = {'seqtype': 'cds', 'func': 'divers', 'mech': 'divers', 'known_prot_snp': ''}
            for item in item_dic.keys():
                if item in tdes.keys():
                    data[item] = tdes[item]
                else:
                    data[item] = item_dic[item]
            if int(data['strand']) < 0:
                data['strand'] = -1
            else:
                data['strand'] = 1
            for item in ['tlen', 'qstart', 'qend', 'tstart', 'tend', 'alen', 'nid', 'npos', 'mismatch', 'gapopen',
                         'gap']:
                data[item] = int(data[item])
            if data['strand'] == -1:
                data['qstart'], data['qend'] = data['qend'], data['qstart']
            for item in ['pid', 'ppos']:
                data[item] = round(float(data[item]), 2)
            data['pcv'] = round(100 * len(data['tseq'].replace('-', '')) / float(len(data['fulltseq'])), 2)
            data['tseq'] = "-"
            dmnd_results.append(data)

    print("***** End loading Diamond Result *****\n")
    return dmnd_results


def overlap_filter(results, taxonomy, pass_overlap=50):
    filtered_results = []
    ctgs = list(set([d['qid'] for d in results]))
    ctgs.sort()
    for ctg in ctgs:
        subset_results = []
        for d in results:
            if d['qid'] == ctg:
                subset_results.append(d)
        print(ctg, len(subset_results), 'features -> kept:', end=' ')
        print("")
        comparisons = list(itertools.combinations(subset_results, 2))
        # print ctg, len(subset_result.keys()), len(comparisons)
        del_list = []
        for data1, data2 in comparisons:
            if subset_results.index(data1) not in del_list and subset_results.index(data2) not in del_list:

                qid1 = data1['qid']
                tid1 = data1['tid']
                pos1 = range(data1['qstart'], data1['qend'] + 1)
                qid2 = data2['qid']
                tid2 = data2['tid']
                pos2 = range(data2['qstart'], data2['qend'] + 1)
                intersection = len([x for x in pos1 if x in pos2])
                if intersection >= pass_overlap:
                    pattern = re.compile('.+_\[(.+)\]::.+')
                    match1 = pattern.match(tid1)
                    match2 = pattern.match(tid2)

                    # case to save mutation in ARM gene (mutation give supplementary resistance)
                    if "known_prot_snp" in data1 and "known_prot_snp" in data2:
                        if len(data1["known_prot_snp"]) >= 1 and len(data2["known_prot_snp"]) == 0:
                            continue
                        elif len(data2["known_prot_snp"]) >= 1 and len(data1["known_prot_snp"]) == 0:
                            continue
                    elif "known_dna_snp" in data1 and "known_dna_snp" in data2:
                        if len(data1["known_dna_snp"]) >= 1 and len(data2["known_dna_snp"]) == 0:
                            continue
                        elif len(data2["known_dna_snp"]) >= 1 and len(data1["known_dna_snp"]) == 0:
                            continue

                    if match1 and match2 and taxonomy != '':
                        if taxonomy.lower() in match1.groups()[0].lower() and \
                                taxonomy.lower() in match2.groups()[0].lower():
                            score1 = (data1['nid'] + data1['npos'] - data1['gap']) / float(data1['tlen'])
                            score2 = (data2['nid'] + data2['npos'] - data2['gap']) / float(data2['tlen'])

                            if score1 >= score2:
                                del_list.append(subset_results.index(data2))
                            else:
                                del_list.append(subset_results.index(data1))
                        elif taxonomy.lower() not in match2.groups()[0].lower():
                            del_list.append(subset_results.index(data2))
                        elif taxonomy.lower() not in match1.groups()[0].lower():
                            del_list.append(subset_results.index(data1))
                    elif match1 and taxonomy != '':
                        if taxonomy.lower() not in match1.groups()[0].lower():
                            del_list.append(subset_results.index(data2))
                    elif match2 and taxonomy != '':
                        if taxonomy.lower() not in match2.groups()[0].lower():
                            del_list.append(subset_results.index(data1))
                    else:
                        score1 = (data1['nid'] + data1['npos'] - data1['gap']) / float(data1['tlen'])
                        score2 = (data2['nid'] + data2['npos'] - data2['gap']) / float(data2['tlen'])
                        if score1 >= score2:
                            del_list.append(subset_results.index(data2))
                        else:
                            del_list.append(subset_results.index(data1))

        del_list = list(set(del_list))
        del_list.sort()
        del_list.reverse()
        for item in del_list:
            del subset_results[item]
        filtered_results = filtered_results + subset_results
    return filtered_results


def taxonomy_filter(results, taxonomy):
    del_list = []
    for index, data in enumerate(results):
        tid = data['tid']
        pattern = re.compile('.+_\[(.+)\]::.+')
        match = pattern.match(tid)
        if match:
            if taxonomy.lower() not in match.groups()[0].lower():
                del_list.append(index)

    del_list = list(set(del_list))
    del_list.sort()
    del_list.reverse()
    for item in del_list:
        del results[item]
    return results


def show_cds_result(dmnd_results):
    for data in dmnd_results:
        print(f'\nTarget: {data["tid"]}\ttarget length: {data["tlen"]}\tstart: {data["tstart"]}\tend: {data["tend"]}')
        print( f'Query: {data["qid"]}\tquery start: {data["qstart"]}\tquery end: {data["qend"]}\tstrand: {data["strand"]}')
        print(f'Perc_pos: {data["ppos"]}\tPerc_id: {data["pid"]}\ttarget_cov: {data["pcv"]}')
        print('Blastx:')
        if data['strand'] > 0:
            try:
                print(f'Detected: {str(Seq(data["qseq"]).translate(table="Bacterial", cds=True))}')
            except Exception as e:
                print(e)
                print(f'Detected: {str(Seq(data["qseq"]).translate(table="Bacterial", cds=False))}')
        else:
            try:
                print(f'Detected: {str(Seq(data["qseq"]).reverse_complement().translate(table="Bacterial", cds=True))}')
            except Exception as e:
                print(e)
                print(
                    f'Detected: {str(Seq(data["qseq"]).reverse_complement().translate(table="Bacterial", cds=False))}')
        if 'qprot' and 'tprot' in data:
            print('Final alignment:')
            print(f'Detected: {data["qprot"]}')
            print(f'Database: {data["tprot"]}')
        if 'mean_qual' in data:
            print(
                f'Max base quality: {data["max_qual"]}\tMean base quality: {data["mean_qual"]}\tMin base quality:'
                f' {data["min_qual"]}')
            print(
                f'Max base depth: {data["max_depth"]}\tMean base depth: {data["mean_depth"]}\tMin base depth:'
                f' {data["min_depth"]}')
        if 'prot_snp' in data:
            snps = []
            for d in data['prot_snp']:
                if 'dna_data' in d:
                    snps.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]} [{d["dna_data"]}]')
                else:
                    snps.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]}')
            print(f'{len(data["prot_snp"])} prot snp(s): {", ".join(snps)}')
        if 'prot_sub' in data:
            subs = []
            for d in data['prot_sub']:
                if 'dna_data' in d:
                    subs.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]} [{d["dna_data"]}]')
                else:
                    subs.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]}')
            print(f'{len(data["prot_sub"])} prot sub(s): {", ".join(subs)}')
        if 'known_prot_snp' in data:
            print(f'DBse snp: {", ".join(data["known_prot_snp"].split("|"))}')
        if 'warning' in data:
            print(f'Sequence warning: {data["warning"]}')
        print('')


def show_dna_result(blastn_results):
    for data in blastn_results:
        print(f'\nTarget: {data["tid"]}\ttarget length: {data["tlen"]}\tstart: {data["tstart"]}\tend: {data["tend"]}')
        print(
            f'Query: {data["qid"]}\tquery start: {data["qstart"]}\tquery end: {data["qend"]}\tstrand: {data["strand"]}')
        print(f'Perc_pos: {data["ppos"]}\tPerc_id: {data["pid"]}\ttarget_cov: {data["pcv"]}')
        print('Blastn:')
        # if data['strand'] > 0:
        print(f'Detected: {data["qseq"]}')
        print(f'DBase   : {data["tseq"]}')
        if 'warning' in data:
            print(f'Sequence warning: {data["warning"]}')
        # else:
        #    print 'Detected: {0}'.format(str(Seq(data['qseq']).reverse_complement()))
        # if 'qdna' and 'tdna' in data:
        #    print 'Final alignment:'
        #    print 'Detected: {0}'.format(data['qdna'])
        #    print 'Database: {0}'.format(data['tdna'])
        if 'mean_qual' in data:
            print(
                f'Max base quality: {data["max_qual"]}\tMean base quality: {data["mean_qual"]}\tMin base quality:'
                f' {data["min_qual"]}')
            print(
                f'Max base depth: {data["max_depth"]}\tMean base depth: {data["mean_depth"]}\tMin base depth:'
                f' {data["min_depth"]}')
        if 'dna_snp' in data:
            snps = []
            for d in data['dna_snp']:
                # print d['dna_data']
                if 'dna_data' in d:
                    snps.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]} [{d["dna_data"]}]')
                else:
                    snps.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]}')
            print(f'{len(data["dna_snp"])} dna snp(s): {", ".join(snps)}')
        if 'dna_sub' in data:
            subs = []
            for d in data['dna_sub']:
                # print d['dna_data']
                if 'dna_data' in d:
                    subs.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]} [{d["dna_data"]}]')
                else:
                    subs.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]}')
            print(f'{len(data["dna_sub"])} dna sub(s): {", ".join(subs)}')
        if 'known_dna_snp' in data:
            print(f'DBse snp: {", ".join(data["known_dna_snp"].split("|"))}')
        print('')


def dna_data_to_dict(id, d, dna_data, item):
    pos = dna_data['p']
    strand = dna_data['f']
    base = dna_data['b']
    depth = dna_data['d']
    qual = dna_data['q']
    type_mut = dna_data['t']
    mut = ""

    if item in ['prot_snp', 'prot_sub']:
        mut = OrderedDict([('Feature', id), ('Mutation', f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]}'),
                           ('Mutation type', type_mut), ('DNA position', pos), ('DNA strand', strand), ('Codon', base),
                           ('Sequencing depth', depth), ('Base quality', qual)])
    elif item in ['dna_snp', 'dna_sub']:
        mut = OrderedDict([('Feature', id), ('Mutation', f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]}'),
                           ('Mutation type', type_mut), ('DNA position', pos), ('DNA strand', strand),
                           ('Codon/base', base),
                           ('Sequencing depth', depth), ('Base quality', qual)])
    return mut


def write_csv_html(merged_results, mut_prefix, id_prefix, pass_alarm_qual=20, pass_alarm_depth=30, pass_alarm_frac=0.9):

    header = ['Function', 'DBase name', '% ident', '% cov', 'Sequence Warning', 'Min depth', 'Mean depth', 'Max depth',
              'Min qual', 'Mean qual', 'Max qual', 'SNP', 'SUB', 'Known prot SNP', 'Known DNA SNP', 'DBase start',
              'DBase end', 'DBase length', 'Query name', 'Query start', 'Query end', 'Query strand', 'Query length',
              'Query DNA seq', 'Query prot seq', 'DBase dna seq', 'DBase prot seq']
    records = []
    for data in merged_results:
        if 'warning' in data:
            warning = data['warning']
        else:
            warning = 'No information'
        values = [data['func'], data['tid'], data['pid'], data['pcv'], warning]
        if 'min_depth' in data:
            values = values + [data['min_depth'], data['mean_depth'], data['max_depth']]
        else:
            values = values + ['-', '-', '-']
        if 'min_qual' in data:
            values = values + [data['min_qual'], data['mean_qual'], data['max_qual']]
        else:
            values = values + ['-', '-', '-']

        mutations = []
        snps = []

        for item in ['prot_snp', 'dna_snp']:
            if item in data:
                for d in data[item]:
                    alarm = []
                    depth_alarm = ''
                    qual_alarm = ''
                    frac_alarm = ''
                    if 'dna_data' in d:
                        if d['dna_data'] != '':
                            for dna_data in d['dna_data'].split(' | '):
                                dna_data = dict(x.split(':') for x in dna_data.strip().split('; '))
                                dna_data['t'] = 'ARM'
                                mut = dna_data_to_dict(data['tid'], d, dna_data, item)
                                mutations.append(mut)

                                if depth_alarm == '':
                                    depths = mut['Sequencing depth']
                                    for depth in depths.split(','):
                                        depth_ref, depth_total = depth.split('/')
                                        if int(depth_ref) < pass_alarm_depth:
                                            depth_alarm = f'depth<{pass_alarm_depth}'
                                            break

                                if frac_alarm == '':
                                    depths = mut['Sequencing depth']
                                    for depth in depths.split(','):
                                        depth_ref, depth_total = depth.split('/')
                                        if int(depth_total) != 0:
                                            frac = int(depth_ref) / float(depth_total)
                                            if frac < pass_alarm_frac:
                                                frac_alarm = f'frac<{pass_alarm_frac}'
                                                break
                                        else:
                                            frac_alarm = 'depth_total:0'
                                            break

                                if qual_alarm == '':
                                    quals = mut['Base quality']
                                    for qual in quals.split(','):
                                        if float(qual) < pass_alarm_qual:
                                            qual_alarm = f'qual<{pass_alarm_qual}'
                                            break

                            for a in [frac_alarm, depth_alarm, qual_alarm]:
                                if a != '':
                                    alarm.append(a)
                            if alarm:
                                alarm = ','.join(alarm)
                            else:
                                alarm = ''
                        else:
                            alarm = ''
                    else:
                        alarm = ''

                    if item == 'prot_snp':
                        if alarm == '':
                            snps.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]}')
                        else:
                            snps.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]} #{alarm}#')
                    elif item == 'dna_snp':
                        if alarm == '':
                            snps.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]}')
                        else:
                            snps.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]} #{alarm}#')

        snps = ', '.join(snps)
        subs = []
        for item in ['prot_sub', 'dna_sub']:
            if item in data:
                for d in data[item]:
                    alarm = []
                    depth_alarm = ''
                    qual_alarm = ''
                    frac_alarm = ''
                    if 'dna_data' in d:
                        if d['dna_data'] != '':
                            for dna_data in d['dna_data'].split(' | '):
                                dna_data = dict(x.split(':') for x in dna_data.strip().split('; '))
                                dna_data['t'] = 'Unknown'
                                mut = dna_data_to_dict(data['tid'], d, dna_data, item)
                                mutations.append(mut)

                                if depth_alarm == '':
                                    depths = mut['Sequencing depth']
                                    for depth in depths.split(','):
                                        depth_ref, depth_total = depth.split('/')
                                        try:
                                            if int(depth_ref) < pass_alarm_depth:
                                                depth_alarm = f'depth<{pass_alarm_depth}'
                                                break
                                        except Exception as e:
                                            print(e)
                                            depth_alarm = 'depth_trbl'
                                            break

                                if frac_alarm == '':
                                    depths = mut['Sequencing depth']
                                    for depth in depths.split(','):
                                        depth_ref, depth_total = depth.split('/')
                                        try:
                                            if int(depth_total) != 0:
                                                frac = int(depth_ref) / float(depth_total)
                                                if frac < pass_alarm_frac:
                                                    frac_alarm = f'frac<{pass_alarm_frac}'
                                                    break
                                            else:
                                                frac_alarm = 'depth_total:0'
                                                break
                                        except Exception as e:
                                            print(e)
                                            frac_alarm = 'frac_trbl'
                                            break

                                if qual_alarm == '':
                                    quals = mut['Base quality']
                                    for qual in quals.split(','):
                                        try:
                                            if float(qual) < pass_alarm_qual:
                                                qual_alarm = f'qual<{pass_alarm_qual}'
                                                break
                                        except Exception as e:
                                            print(e)
                                            qual_alarm = 'qual_trbl'
                                            break

                            for a in [frac_alarm, depth_alarm, qual_alarm]:
                                if a != '':
                                    alarm.append(a)
                            if alarm:
                                alarm = ','.join(alarm)
                            else:
                                alarm = ''
                        else:
                            alarm = ''
                    else:
                        alarm = ''

                    if item == 'prot_sub':
                        if alarm == '':
                            subs.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]}')
                        else:
                            subs.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]} #{alarm}#')
                    elif item == 'dna_sub':
                        if alarm == '':
                            subs.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]}')
                        else:
                            subs.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]} #{alarm}#')

                    if mutations:
                        df = pd.DataFrame.from_dict(mutations)
                        tid = data['tid'].replace(':', '_').replace('(', '').replace(')', '') .replace('\'', 'pr')
                        df.to_html(f"{mut_prefix}_{tid}_{data['qid']}_{data['qstart']}-{data['qend']}_mut.html",
                                   index=False)
                        df.to_csv(f'{mut_prefix}_{tid}_{data["qid"]}_{data["qstart"]}-{data["qend"]}_mut.csv', sep='\t',
                                  index=False)
        subs = ', '.join(subs)
        if 'known_prot_snp' in data:
            known_prot_snp = data['known_prot_snp'].replace('|', ', ')
        else:
            known_prot_snp = ''

        if 'known_dna_snp' in data:
            known_dna_snp = data['known_dna_snp'].replace('|', ', ')
        else:
            known_dna_snp = ''

        values = values + [snps, subs, known_prot_snp, known_dna_snp]
        values = values + [data['tstart'], data['tend'], data['tlen']]
        values = values + [data['qid'], data['qstart'], data['qend'], data['strand'], data['qend'] - data['qstart'] + 1]

        if 'tprot' in data:
            t_prot_seq = str(data['tprot'])
        else:
            t_prot_seq = ''

        if 'qprot' in data:
            q_prot_seq = str(data['qprot'])
        else:
            q_prot_seq = ''

        if data['strand'] > 0:
            q_dna_seq = str(data['qseq'])
        else:
            q_dna_seq = str(Seq(data['qseq']).reverse_complement())

        if 'tseq' in data:
            t_dna_seq = str(data['tseq'])
        else:
            t_dna_seq = ''

        values = values + [q_dna_seq, q_prot_seq, t_dna_seq, t_prot_seq]
        d = dict(zip(header, values))
        records.append(d)
    df = pd.DataFrame.from_dict(records)
    df.sort_values(['Function', 'DBase name'], ascending=[True, True], inplace=True)
    df[header].to_csv(f'{id_prefix}_results.csv', sep='\t', index=False)
    df[header].to_html(f'{id_prefix}_results.html', index=False)


def description(data):
    snps = []
    if 'prot_snp' in data:
        for d in data['prot_snp']:
            snps.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]}')
        snps = '|'.join(snps)
    elif 'dna_snp' in data:
        for d in data['dna_snp']:
            snps.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]}')
        snps = '|'.join(snps)

    subs = []
    if 'prot_sub' in data:
        for d in data['prot_sub']:
            subs.append(f'{d["t_aa"]}{d["t_prot_pos"]}{d["q_aa"]}')
        subs = '|'.join(subs)
    elif 'dna_sub' in data:
        for d in data['dna_sub']:
            subs.append(f'{d["t_base"]}{d["t_dna_pos"]}{d["q_base"]}')
        subs = '|'.join(subs)

    try:
        des = f'function:{data["func"]}, mechanism:{data["mech"]}, reference_sequence: {data["tid"].split("::")[-1]},' \
              f' perc_identity:{data["pid"]}, perc_coverage:{data["pcv"]}, min_depth:{data["min_depth"]},' \
              f' mean_depth:{data["mean_depth"]}, max_depth:{data["max_depth"]}, min_quality:{data["min_qual"]},' \
              f' mean_quality:{data["mean_qual"]}, max_quality:{data["max_qual"]}, known_sub:{snps}'
    except KeyError:
        des = f'function:{data["func"]}, mechanism:{data["mech"]}, reference_sequence: {data["tid"].split("::")[-1]},' \
              f' perc_identity:{data["pid"]}, perc_coverage:{data["pcv"]}, known_sub:{snps}'
    return des


def test_cds(data):
    test = 1
    if data['strand'] > 0:
        q_dna_seq = Seq(str(data['qseq']).replace('-', ''))
    else:
        q_dna_seq = Seq(str(data['qseq']).replace('-', '')).reverse_complement()
    try:
        x = q_dna_seq.translate(table='Bacterial', cds=True)
    except Exception as e:
        print(e)

        test = 0
    return test


def write_fasta(results, out_dir, out_prefix):
    aa_outfile = os.path.join(out_dir, f'{out_prefix}.faa')
    dna_outfile = os.path.join(out_dir, f'{out_prefix}.fna')
    aa_records_list = []
    dna_records_list = []
    for data in results:
        ID = f'{data["qid"]}__{data["qstart"]}__{data["qend"]}__f{data["strand"]}__{data["tid"]}'
        des = description(data)

        if 'qprot' in data:
            q_prot_seq = str(data['qprot']).replace('-', '')
            aa_rec = SeqRecord(Seq(q_prot_seq), id=ID, description=des)
            aa_records_list.append(aa_rec)

            if data['strand'] > 0:
                q_dna_seq = Seq(data['qseq'].replace('-', ''))
            else:
                q_dna_seq = Seq(data['qseq'].replace('-', '')).reverse_complement()
            dna_rec = SeqRecord(q_dna_seq, id=ID, description=des)
            dna_records_list.append(dna_rec)
        else:
            dna_rec = SeqRecord(Seq(str(data['qseq']).replace('-', '')), id=ID, description=des)
            dna_records_list.append(dna_rec)

    with open(aa_outfile, 'w') as aa_out_f, open(dna_outfile, 'w') as dna_out_f:
        SeqIO.write(aa_records_list, aa_out_f, 'fasta')
        SeqIO.write(dna_records_list, dna_out_f, 'fasta')


def write_gbk(results, query_dic, out_dir, out_prefix):
    rec_dic = {}
    results.sort(key=operator.itemgetter('qid'))
    for data in results:
        key = data['qid']
        if key not in rec_dic:
            rec_dic[key] = [data]
        else:
            rec_dic[key].append(data)

    # list keys of the dic
    keys_list = list(rec_dic.keys())

    # sort list of keys
    keys = sorted(keys_list)

    n = 0
    for key in keys:
        records = rec_dic[key]
        rec = SeqRecord(Seq(str(query_dic[key].seq)), id=key, name=key, description='',
                        annotations={"molecule_type": "DNA"})
        for data in records:
            if test_cds(data) == 1:
                feature = SeqFeature(FeatureLocation(data['qstart'] - 1, data['qend'], strand=data['strand']),
                                     type='CDS', qualifiers={})
            else:
                feature = SeqFeature(FeatureLocation(data['qstart'] - 1, data['qend'], strand=data['strand']),
                                     type='misc_feature', qualifiers={})

            if 'qprot' in data:
                # feature.qualifiers = {'locus_tag':'{0}_{1}'.format(out_prefix, n), 'product':data['tid'].
                # split('::')[0],
                #                  'note':description(data), 'translation':data['qprot']}
                feature.qualifiers = OrderedDict([('product', data['tid'].split('::')[0]),
                                                  ('note', description(data)), ('translation', data['qprot'])])
            else:
                # feature.qualifiers = {'locus_tag':'{0}_{1}'.format(out_prefix, n),
                # 'product':data['tid'].split('::')[0],
                #                   'note':description(data)}
                feature.qualifiers = OrderedDict([('product', data['tid'].split('::')[0]),
                                                  ('note', description(data))])
            rec.features.append(feature)

        rec.features = sorted(rec.features, key=lambda feature: feature.location.start)

        for feature in rec.features:
            feature.qualifiers = OrderedDict([('locus_tag', f'DET_{n + 1}')] + list(feature.qualifiers.items()))
            n += 1

        out_file = os.path.join(out_dir, f'{str(rec.id).replace(".", "-")}.gbk')
        with open(out_file, 'w') as out_f:
            SeqIO.write([rec], out_f, 'genbank')


def main(args):
    print(f"Version diamDetector: {version()}\n")
    wk_dir = args.wkdir
    # Check working directory
    if glob.glob(wk_dir) is None:
        print(f"\nDirectory {wk_dir} not found!\n")
        exit(1)
    else:
        print(f"\nDirectory {wk_dir} found")

    sample_file = args.sample_file

    # Assign default sample file
    if not sample_file:
        sample_file = os.path.join(os.path.dirname(wk_dir), 'sample.tsv')

    # Check sample file
    if not os.path.exists(sample_file):
        print(f"\nSample file {sample_file} not found!\n")
        exit(1)
    else:
        print(f"\nSample file: {sample_file}")

    # Load sample name and species if provided (not required)
    sample_dic = load_sample_name(sample_file)

    # Search the assembly files corresponding to the sample names in the working directory
    # (e.g. ./working_directory/*/sample_name.fasta)
    print("\nSearch the assembly files:")
    file_list = []
    for sample_name in sample_dic.keys():
        try:
            print(f"Search for the assembly file: {os.path.join(wk_dir, f'{sample_name}.fasta')}")
            sample_file = os.path.join(wk_dir, f'{sample_name}.fasta')
            file_list.append(sample_file)
            print(f'Assembly for {sample_name}: {sample_file}')
        except IndexError:
            print(f'Assembly for {sample_name}: file not found!')

    # Check the location of databases:
    print("\nSearch for the databases:")
    check_db = 0
    cds_process = True
    cds_target_file = os.path.abspath(args.cds_target_file)
    if os.path.exists(cds_target_file):
        if os.stat(cds_target_file).st_size == 0:
            print(f"  CDS database {cds_target_file} found but EMPTY !!!!!!!!")
            cds_process = False
        else:
            check_db += 1
            print(f"  CDS database {cds_target_file} found")
    else:
        print("  No CDS database")
    dna_target_file = os.path.abspath(args.dna_target_file)
    dna_process = True
    if os.path.exists(dna_target_file):
        if os.stat(dna_target_file).st_size == 0:
            print(f"  DNA database {dna_target_file} found but EMPTY !!!!!!!!")
            dna_process = False
        else:
            check_db += 1
            print(f"  DNA database {dna_target_file} found")
    else:
        print("  No DNA database")
    if check_db == 0:
        print("NO database provided\n")
        exit(1)
    else:
        print(f"{check_db} database provided\n")

    database_split = os.path.basename(cds_target_file).split(".")[0].split("_")

    database = ""
    if "dbVIR" in dna_target_file:
        database = "VIR"
    elif "dbARM" in dna_target_file:
        database = "ARM"

    # Print version and subset database
    print(f"\nList {database}-DB Subset: {database_split[2]}\n")
    print(f"\nVersion {database}-DB: {database_split[1]}\n")

    print("\nFeature detection parameters:")
    pass_pid = float(args.perc_id)
    print(f"  Identity threshold: {pass_pid}%")
    pass_pcv = float(args.perc_cv)
    print(f"  Coverage threshold: {pass_pcv}%")
    pass_overlap = int(args.overlap)
    print(f"  Maximum feature overlap: {pass_overlap}-bases\n")

    # Number of threads
    threads = int(args.threads)
    # Force the overwrite
    force = args.Force

    # Start the detection of features in databases
    for n, query_file in enumerate(file_list):
        # Extract sample names from file name
        sample_name = os.path.splitext(os.path.basename(query_file))[0]
        # Set bam file
        bam_file = os.path.splitext(query_file)[0] + '.bam'
        # Set taxonomy for taxonomy-based filtering
        taxonomy = sample_dic[sample_name]
        title = f"~~~~  {n + 1}/{len(file_list)} sample: {sample_name} species: {taxonomy}  ~~~~"
        print("\n\n")
        print("~" * len(title))
        print(title)
        print("~" * len(title))
        # Set taxonomy filtering level
        taxonomy_filter_detect = args.taxonomy_filter
        if taxonomy_filter_detect not in ['strict', 'lax', 'none']:
            taxonomy_filter_detect = 'lax'
        if taxonomy == 'No taxonomy provided':
            taxonomy_filter_detect = 'none'
            print('No taxonomy provided: none taxonomy filtering will be performed')

        out_diamond_file = os.path.join(os.path.dirname(query_file), f'diam_output_{sample_name}_{database}.csv')
        # Launch CDS detection
        if os.path.exists(cds_target_file) and cds_process:
            dmnd_db = make_dmnd_database(cds_target_file, force)
            dmnd_result_file = run_diam(dmnd_db, query_file, pass_pid, pass_pcv, threads, force, out_diamond_file)
            dmnd_results = load_dmnd_result(dmnd_result_file, cds_target_file)
        else:
            print('No CDS to search')
            dmnd_results = []

        out_blastn_file = os.path.join(os.path.dirname(query_file), f'blastn_output_{sample_name}_{database}.csv')
        # Launch DNA detection
        if os.path.exists(dna_target_file) and dna_process:
            blastn_db = make_blastn_database(dna_target_file, force)
            blastn_result_file = run_blastn(blastn_db, query_file, pass_pid, force, 0.0001, 8, out_blastn_file)
            blastn_results = load_blastn_result(blastn_result_file, dna_target_file, pass_pid, pass_pcv)
        else:
            print('No DNA sequence to search')
            blastn_results = []

        # Filter the DNA and CDS features on the bases of taxonomy and overlaps
        print(f'\nNumber of detected features: {len(dmnd_results) + len(blastn_results)}')
        if taxonomy_filter_detect == 'strict' and taxonomy != '':
            dmnd_results = taxonomy_filter_detect(dmnd_results, taxonomy)
            blastn_results = taxonomy_filter_detect(blastn_results, taxonomy)
            print(f'Number of detected features after stric taxonomic filtering:  '
                  f'{len(dmnd_results) + len(blastn_results)}')

        if taxonomy_filter_detect == 'lax':
            dmnd_results = overlap_filter(dmnd_results, taxonomy, pass_overlap)
            print(f'\nNumber of detected features after lax taxonomic and overlap filtering: {len(dmnd_results)} \n')
            blastn_results = overlap_filter(blastn_results, taxonomy, pass_overlap)
            print(f'\nNumber of detected features after lax taxonomic and overlap filtering: {len(blastn_results)}')
            print(
                f'\n######## Number of detected features after lax taxonomic and overlap filtering:'
                f' {len(dmnd_results) + len(blastn_results)} ########')
        else:
            dmnd_results = overlap_filter(dmnd_results, '', pass_overlap)
            blastn_results = overlap_filter(blastn_results, '', pass_overlap)
            print(f'Number of detected features after overlap filtering: {len(dmnd_results) + len(blastn_results)}')
        print('')

        # Set the prefix of the output
        if args.outPrefix == '':
            out_prefix = os.path.splitext(os.path.basename(cds_target_file))[0]
        else:
            out_prefix = args.outPrefix

        query_dic = load_fasta(query_file)
        if len(dmnd_results) > 0:
            # Global alignment of CDS and mutation extraction if CDS features detected
            dmnd_results = cds_global_alignment(dmnd_results, query_dic, wk_dir, pass_pid, pass_pcv)
            if os.path.exists(bam_file):
                # Extraction quality of bases and sequencing depth if bam detected
                dmnd_results = cds_extract_quality_and_depth(bam_file, query_file, dmnd_results, out_prefix, force)
            # Show the detected CDS features
            show_cds_result(dmnd_results)

        if len(blastn_results) > 0:
            # Global alignment of DNA and mutation extraction if DNA features detected
            target_dic = load_fasta(dna_target_file)
            blastn_results = dna_global_alignemnt(blastn_results, query_dic, target_dic, pass_pid, pass_pcv)
            if os.path.exists(bam_file):
                # Extraction quality of bases and sequencing depth if bam detected
                blastn_results = dna_extract_quality_and_depth(bam_file, query_file, blastn_results, out_prefix, force)
            # Show the detected DNA features
            show_dna_result(blastn_results)

        # Merge DNA and CDS feature result files
        merged_results = dmnd_results + blastn_results

        if len(merged_results) > 0:
            # Set the directory to store depth and quality data
            mut_dir = os.path.join(os.path.dirname(bam_file), 'Mutations_depth_quality')
            if not os.path.exists(mut_dir):
                os.mkdir(mut_dir)
            mut_prefix = os.path.join(mut_dir, os.path.splitext(os.path.basename(query_file))[0].split('_')[0])
            id_prefix = os.path.splitext(query_file)[0]
            # Write the results as tsv and html files
            write_csv_html(merged_results, mut_prefix, id_prefix)
            # Set the directory to store the sequences
            out_dir = os.path.join(os.path.dirname(query_file), 'Detected_sequences')
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            print("\nStart to write gbk\n")
            # Write the sequences as gbk and fasta files
            write_gbk(merged_results, query_dic, out_dir, out_prefix)
            print("\nStart to write fasta\n")
            write_fasta(merged_results, out_dir, out_prefix)

            if n >= 2:
                # remove file
                if os.path.exists(out_diamond_file):
                    os.remove(out_diamond_file)
                if os.path.exists(out_blastn_file):
                    os.remove(out_blastn_file)
                shutil.rmtree(out_dir)
                shutil.rmtree(mut_dir)
        else:
            print('\nNo results!\n')


def version():
    return "1.0"


def run():
    global usage
    usage = "diamDetector.py [-cdsDB cds database 'armDB_6_GN_cds.fas'] [-dnaDB dna database 'armDB_6_GN_dna.fas']" \
            " [-wd work directory with the assembly] [-sf sample file] "

    parser = argparse.ArgumentParser(
        prog='diamDetector',
        usage=usage,
        description='DiamDetector: pipeline CNR Resistance for diamond Detection - Version ' + version(),
    )

    parser.add_argument('-sf', '--sample_file', dest="sample_file",
                        help='Tab-separated file containing the names of sample and the corresponding taxonomy '
                             '(species/genus) [$wkdir/../sample.csv]')
    parser.add_argument('-wd', '--wkdir', dest="wkdir", help='Directory containing the assembly files')
    parser.add_argument('-cdsDB', '--cdsDBFile', dest="cds_target_file",
                        default="/usr/local/readmapper-v0.1/dbARM/armDB_6_GN_cds.fas",
                        help='CDS Database file in fasta format')
    parser.add_argument('-dnaDB', '--dnaDBFile', dest="dna_target_file",
                        default="/usr/local/readmapper-v0.1/dbARM/armDB_6_GN_dna.fas",
                        help='DNA Database file in fasta format')
    parser.add_argument('-id', '--perc_id', dest="perc_id", default="70",
                        help="Minimum identity percentage with the target [70]")
    parser.add_argument('-cv', '--perc_cv', dest="perc_cv", default="70",
                        help="Minimum target coverage percentage [70]")
    parser.add_argument('-ov', '--overlap', dest="overlap", default="20",
                        help="Maximun overlap between detected features as base number [20]")
    parser.add_argument('-tf', '--taxonFilter', dest="taxonomy_filter", default='lax',
                        help='taxonomy filter type: strict, lax or none')
    parser.add_argument('-th', '--threads', dest="threads", default='8', help='Number of threads [8]')
    parser.add_argument('-F', '--Force', dest="Force", action="store_true", default=False,
                        help="Overwrite the detected files")
    parser.add_argument('-o', '--outPrefix', dest="outPrefix", default='', help="Outprefix [<database_name>]")
    parser.add_argument('-V', '--version', action='version', version='diamDetector-' + version(),
                        help="Prints version number")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    main(args)


if __name__ == '__main__':
    run()
