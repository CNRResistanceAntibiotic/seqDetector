#!/usr/bin/python3
import os
import subprocess

from seqdetector import bam2data
from seqdetector.misc_functions import load_fasta, read_bam_count, read_bam_stat, reverse_complement


def make_blastn_database(fasta_file, force=False, threads=8):
    blastn_db = os.path.abspath(fasta_file)
    if os.path.exists(blastn_db + '.nsp') and os.path.exists(blastn_db + '.nin') and os.path.exists(blastn_db + '.nhr') \
            and not force:
        print('\nDatabase {0} already exists'.format(blastn_db))
    else:
        cmd = '$(which makeblastdb) -in {0} -dbtype nucl -out {1}'.format(fasta_file, blastn_db)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
        print("\n{0}\n{1}".format(cmd, process.decode("utf-8")))
    return blastn_db


def run_blastn(blastn_db, query_file, pass_pid=70, force=True, evalue=0.0001, threads=8, out_file=""):

    if not force and os.path.exists(out_file):
        print('\nResult file {0} already exists'.format(out_file))

    else:
        fmt = '\"6 qseqid frames sallseqid slen qstart qend sstart send length pident nident ppos positive mismatch ' \
              'gapopen gaps qseq sseq\"'
        cmd = '$(which blastn) -out {0} -outfmt {1} -query {2} -db {3} -num_threads {4} -perc_identity {5} -evalue {6}' \
            .format(out_file, fmt, query_file, blastn_db, threads, pass_pid, evalue)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
        print("\n{0}\n{1}".format(cmd, process.decode("utf-8")))
    return out_file


def load_blastn_result(result_file, target_file, pass_pid=70, pass_pcv=70):
    target_dic = load_fasta(target_file)
    blastn_results = []
    header = 'qid strand tid tlen qstart qend tstart tend ' \
             'alen pid nid ppos npos mismatch gapopen gap ' \
             'qseq tseq'.split(' ')
    with open(result_file) as inf:
        for line in inf:
            data = dict(zip(header, line.strip().split('\t')))
            q_id = data['qid']
            t_id = data['tid']
            t_des = {}
            """
            print(target_dic[t_id].id)
            for item in target_dic[t_id].description.split(';'):
                try:
                    key, value = item.split(':')
                    if key == 'dna_snp':
                        key = 'known_dna_snp'
                    t_des[key] = value
                except Exception as e:
                    print(e)
                    t_des = {}
            print("mama")
            """
            item_dic = {'seqtype': 'dna', 'func': 'divers', 'mech': 'divers', 'known_dna_snp': ''}
            for item in item_dic.keys():
                if item in t_des.keys():
                    data[item] = t_des[item]
                else:
                    data[item] = item_dic[item]

            q_frame = int(data['strand'].split('/')[0])
            t_frame = int(data['strand'].split('/')[1])

            if q_frame > 0 and t_frame > 0:
                data['strand'] = 1
            elif q_frame < 0 and t_frame > 0:
                data['strand'] = -1
                data['qstart'], data['qend'] = data['qend'], data['qstart']
            elif q_frame > 0 and t_frame < 0:
                data['strand'] = -1
                data['tstart'], data['tend'] = data['tend'], data['tstart']
                data['tseq'] = reverse_complement(data['tseq'])
                data['qseq'] = reverse_complement(data['qseq'])
            elif q_frame < 0 and t_frame < 0:
                data['strand'] = 1
                data['tstart'], data['tend'] = data['tend'], data['tstart']
                data['qstart'], data['qend'] = data['qend'], data['qstart']
                data['tseq'] = reverse_complement(data['tseq'])
                data['qseq'] = reverse_complement(data['qseq'])

            for item in ['tlen', 'qstart', 'qend', 'tstart', 'tend', 'alen', 'nid', 'npos', 'mismatch', 'gapopen',
                         'gap']:
                data[item] = int(data[item])
            for item in ['pid', 'ppos']:
                data[item] = round(float(data[item]), 2)

            data['pcv'] = round(100 * len(data['qseq'].replace('-', '')) / float(len(target_dic[t_id].seq)), 2)

            if data['pcv'] >= pass_pcv and data['pid'] >= pass_pid:
                blastn_results.append(data)
    return blastn_results


def dna_extract_quality_and_depth(bam_file, fas_file, blastn_results, out_prefix, force):
    for data in blastn_results:
        feature_name = data['tid'].replace(':', '_').replace('(', '').replace(')', '').replace('\'', 'pr')
        ctg = data['qid']
        strand = data['strand']
        start = data['qstart']
        end = data['qend']
        position = '{0}:{1}-{2}'.format(ctg, start, end)
        out_dir = os.path.join(os.path.dirname(bam_file), '{}_depth_quality'.format(out_prefix))
        stat_outfile = '{0}_{1}_{2}_{3}-{4}_stats.csv'.format(
            os.path.join(out_dir, os.path.splitext(os.path.basename(bam_file))[0].split('_')[0]),
            feature_name, ctg, start, end)

        data_outfile = '{0}_{1}_{2}_{3}-{4}_count.csv'.format(
            os.path.join(out_dir, os.path.splitext(os.path.basename(bam_file))[0].split('_')[0]),
            feature_name, ctg, start, end)

        feature_name = feature_name + "_{0}_{1}-{2}".format(ctg, start, end)

        if not (os.path.exists(data_outfile) and not os.path.exists(stat_outfile)) or force:

            bam2data.main(bam_file=bam_file, fasta_ref=fas_file, position=position, output_dir=out_dir,
                          feature_name=feature_name, force=True, mapping_qual=0, base_qual=0, stats=True, data=True)

        result = read_bam_stat(stat_outfile)
        data['mean_qual'] = round(float(result['Ref_quali_mean']), 1)
        data['min_qual'] = int(float(result['Ref_quali_min']))
        data['max_qual'] = int(float(result['Ref_quali_max']))
        data['mean_depth'] = round(float(result['Ref_depth_mean']), 1)
        data['min_depth'] = int(float(result['Ref_depth_min']))
        data['max_depth'] = int(float(result['Ref_depth_max']))

        result, warning = read_bam_count(data_outfile)
        data['warning'] = warning
        for item in ['dna_snp', 'dna_sub']:
            if data[item]:
                for snp in data[item]:
                    dna_pos = snp['q_dna_pos']
                    t_base = snp['t_base']
                    q_base = snp['q_base']

                    if strand > 0 and q_base != 'd':
                        dna_pos = start + (dna_pos - 1)
                        if str(dna_pos) in result[ctg]:
                            d = result[ctg][str(dna_pos)]
                            ref = d['reference'].upper()
                            if ref in ['A', 'T', 'C', 'G']:
                                ref_depth = '{0}/{1}'.format(d['{0}_depth'.format(ref)], d['total_depth'])
                                ref_qual = str(round(float(d['{0}_quality'.format(ref)]), 1))
                            else:
                                continue

                            txt = 'p:{0}; f:{1}; b:{2}; d:{3}; q:{4}'.format(dna_pos, strand, ref, ref_depth, ref_qual)
                        else:
                            txt = ''

                    elif strand < 0 and q_base != 'd':
                        dna_pos = end - (dna_pos + 1)
                        # depth known at "dna_pos" position
                        if str(dna_pos) in result[ctg]:
                            d = result[ctg][str(dna_pos)]
                            ref = d['reference'].upper()
                            if ref in ['A', 'T', 'C', 'G']:
                                ref_depth = '{0}/{1}'.format(d['{}_depth'.format(ref)], d['total_depth'])
                                ref_qual = str(round(float(d['{}_quality'.format(ref)]), 1))
                            else:
                                continue
                            txt = 'p:{0}; f:{1}; b:{2}; d:{3}; q:{4}'.format(dna_pos, strand, ref, ref_depth, ref_qual)
                        # depth unknown at "dna_pos" position
                        else:
                            txt = ''
                    else:
                        txt = ''
                    snp['dna_data'] = txt

    return blastn_results


def cds_extract_quality_and_depth(bam_file, fas_file, dmnd_results, out_prefix, force):
    for data in dmnd_results:
        feature_name = data['tid'].replace(':', '_').replace('(', '').replace(')', '').replace('\'', 'pr')
        ctg = data['qid']
        strand = data['strand']
        start = data['qstart']
        end = data['qend']
        position = '{0}:{1}-{2}'.format(ctg, start, end)
        out_dir = os.path.join(os.path.dirname(bam_file), '{0}_depth_quality'.format(out_prefix))

        feature_name = feature_name + "_{0}_{1}-{2}".format(ctg, start, end)

        bam2data.main(bam_file=bam_file, fasta_ref=fas_file, position=position, output_dir=out_dir,
                      feature_name=feature_name, force=True, mapping_qual=0, base_qual=0, stats=True, data=True)

        stat_outfile = '{0}_{1}_stats.csv'.format(
            os.path.join(out_dir, os.path.splitext(os.path.basename(bam_file))[0].split('_')[0]),
            feature_name, ctg, start, end)

        result = read_bam_stat(stat_outfile)
        data['mean_qual'] = round(float(result['Ref_quali_mean']), 1)
        data['min_qual'] = int(float(result['Ref_quali_min']))
        data['max_qual'] = int(float(result['Ref_quali_max']))
        data['mean_depth'] = round(float(result['Ref_depth_mean']), 1)
        data['min_depth'] = int(float(result['Ref_depth_min']))
        data['max_depth'] = int(float(result['Ref_depth_max']))

        data_outfile = '{0}_{1}_count.csv'.format(
            os.path.join(out_dir, os.path.splitext(os.path.basename(bam_file))[0].split('_')[0]),
            feature_name, ctg, start, end)
        # if os.path.exists(data_outfile) == False or force == True:

        result, warning = read_bam_count(data_outfile)

        data['warning'] = warning
        for item in ['prot_snp', 'prot_sub']:
            if data[item]:
                for snp in data[item]:
                    # print data['tid'], snp, strand
                    prot_pos = snp['q_prot_pos']
                    t_aa = snp['t_aa']
                    q_aa = snp['q_aa']

                    if q_aa == 'd':
                        size = len(t_aa.replace('[', '').replace(']', ''))
                    else:
                        size = len(q_aa.replace('[', '').replace(']', ''))
                    if strand > 0 and q_aa != 'd':
                        dna_pos = start + ((prot_pos - 1) * 3)
                        dna_pos = range(dna_pos, dna_pos + (size * 3))
                        bases, quals, depths = [], [], []
                        base, qual, depth = 0, 0, 0
                        # print data['tid']
                        for n, pos in enumerate(dna_pos):
                            # x = result[ctg].keys()
                            # x.sort()
                            # print x
                            try:
                                d = result[ctg][str(pos)]
                                ref = d['reference'].upper()
                                if ref in ['A', 'T', 'C', 'G']:
                                    ref_depth = '{0}/{1}'.format(d['{0}_depth'.format(ref)], d['total_depth'])
                                    ref_qual = str(round(float(d['{0}_quality'.format(ref)]), 1))

                                else:
                                    continue

                            except KeyError:
                                print('Position {0} in {1} for feature {2} with reference depth of 0!'
                                      .format(pos, data['qid'], data['tid']))
                                ref = 'N'
                                ref_depth = '{0}/{1}'.format('0', '0')
                                ref_qual = '0'

                            if n % 3 == 0:
                                if base != 0:
                                    bases.append(base)
                                    quals.append(qual)
                                    depths.append(depth)
                                base, qual, depth = 'p:{0}; f:{1}; b:{2}'.format(pos, strand, ref), ref_qual, ref_depth
                            else:
                                if base == 0:
                                    base, qual, depth = 'p:{0}; f:{1}; b:{2}'.format(pos, strand,
                                                                                     ref), ref_qual, ref_depth
                                else:
                                    base = base + ref
                                    depth = depth + ',' + ref_depth
                                    qual = qual + ',' + ref_qual

                        if depth == 0:
                            continue

                        bases.append(base)
                        quals.append(qual)
                        depths.append(depth)
                        txt = ''
                        for i in range(0, len(bases)):
                            txt = txt + '{0}; d:{1}; q:{2} | '.format(bases[i], depths[i], quals[i])
                        txt = txt[:-3]

                    elif strand < 0 and q_aa != 'd':
                        dna_pos = end - ((prot_pos - 1) * 3) + 1
                        dna_pos = range(dna_pos - (size * 3), dna_pos)
                        bases, quals, depths = [], [], []
                        base, qual, depth = 0, 0, 0
                        for n, pos in enumerate(dna_pos):
                            try:
                                d = result[ctg][str(pos)]
                                ref = d['reference'].upper()
                                if ref in ['A', 'T', 'C', 'G']:
                                    ref_depth = '{0}/{1}'.format(d['{0}_depth'.format(ref)], d['total_depth'])
                                    ref_qual = str(round(float(d['{0}_quality'.format(ref)]), 1))

                                else:
                                    continue

                            except KeyError:
                                print('Position {0} in {1} for feature {2} with reference depth of 0!'
                                      .format(pos, data['qid'], data['tid']))
                                ref = 'N'
                                ref_depth = '{0}/{1}'.format('0', '0')
                                ref_qual = '0'

                            if n % 3 == 0:
                                if base != 0:
                                    bases.append(base)
                                    quals.append(qual)
                                    depths.append(depth)
                                base, qual, depth = 'p:{0}; f:{1}; b:{2}'.format(pos, strand, ref), ref_qual, ref_depth
                            else:
                                if base == 0:
                                    base, qual, depth = 'p:{0}; f:{1}; b:{2}'.format(pos, strand,
                                                                                     ref), ref_qual, ref_depth
                                else:
                                    base = base + ref
                                    depth = depth + ',' + ref_depth
                                    qual = qual + ',' + ref_qual

                        if depth == 0:
                            continue

                        bases.append(base)
                        quals.append(qual)
                        depths.append(depth)
                        txt = ''
                        for i in range(0, len(bases)):
                            txt = txt + '{0}; d:{1}; q:{2} | '.format(bases[i], depths[i], quals[i])
                        txt = txt[:-3]
                    else:
                        txt = ''
                    snp['dna_data'] = txt
    return dmnd_results
