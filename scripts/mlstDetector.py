#!/usr/bin/python3
import os, glob, re, itertools, argparse
import sys
import numpy as np
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
import pandas as pd
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature

from seqdetector.blast_functions import make_blastn_database, run_blastn, load_blastn_result,\
    dna_extract_quality_and_depth
from seqdetector.misc_functions import load_sample_name, load_fasta, dna_global_alignemnt


def load_set_file(set_file, sep='\t'):
    set_dic = {}
    with open(set_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line != '':
                key, data = line.split(sep)
                set_dic[key.lower()] = {}
                data_list = data.split(',')
                for data in data_list:
                    key2, data2 = data.split(':')
                    set_dic[key.lower()][key2] = data2.split('|')
    return set_dic


def overlap_filter(results, pass_overlap=50):
    filtered_results = []
    ctgs = list(set([d['qid'] for d in results]))
    ctgs.sort()
    for ctg in ctgs:
        subset_results = []
        for d in results:
            if d['qid'] == ctg:
                subset_results.append(d)
        print(ctg, len(subset_results), 'features -> kept:', end=' ')
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
                    score1 = (data1['nid'] - data1['gap']) / float(data1['tlen'])
                    score2 = (data2['nid'] - data2['gap']) / float(data2['tlen'])
                    if score1 >= score2:
                        del_list.append(subset_results.index(data2))
                    else:
                        del_list.append(subset_results.index(data1))

        del_list = list(set(del_list))
        del_list.sort()
        del_list.reverse()
        for item in del_list:
            del subset_results[item]
        print(len(subset_results))
        filtered_results = filtered_results + subset_results

    return filtered_results


def view_dna_result(blastn_results):
    for data in blastn_results:
        print('\nTarget: {0}\ttarget length: {1}\tstart: {2}\tend: {3}'.format(data['tid'], data['tlen'],
                                                                               data['tstart'], data['tend']))
        print('Query: {0}\tquery start: {1}\tquery end: {2}\tstrand: {3}'.format(data['qid'], data['qstart'], data['qend'],
                                                                             data['strand']))
        print('Perc_pos: {0}\tPerc_id: {1}\ttarget_cov: {2}'.format(data['ppos'], data['pid'], data['pcv']))
        print('Blastn:')
        # if data['strand'] > 0:
        print('Detected: {0}'.format(data['qseq']))
        print('DBase   : {0}'.format(data['tseq']))
        if 'warning' in data:
            print('Sequence warning: {0}'.format(data['warning']))

        if 'mean_qual' in data:
            print('Max base quality: {0}\tMean base quality: {1}\tMin base quality: {2}'
                  .format(data['max_qual'], data['mean_qual'], data['min_qual']))
            print('Max base depth: {0}\tMean base depth: {1}\tMin base depth: {2}'
                  .format(data['max_depth'], data['mean_depth'], data['min_depth']))
        if 'dna_snp' in data:
            snps = []
            for d in data['dna_snp']:
                # print d['dna_data']
                if 'dna_data' in d:
                    snps.append('{0}{1}{2} [{3}]'.format(d['t_base'], d['t_dna_pos'], d['q_base'], d['dna_data']))
                else:
                    snps.append('{0}{1}{2}'.format(d['t_base'], d['t_dna_pos'], d['q_base']))
            print('{0} dna snp(s): {1}'.format(len(data['dna_snp']), ', '.join(snps)))
        if 'dna_sub' in data:
            subs = []
            for d in data['dna_sub']:
                # print d['dna_data']
                if 'dna_data' in d:
                    subs.append('{0}{1}{2} [{3}]'.format(d['t_base'], d['t_dna_pos'], d['q_base'], d['dna_data']))
                else:
                    subs.append('{0}{1}{2}'.format(d['t_base'], d['t_dna_pos'], d['q_base']))
            print('{0} dna sub(s): {1}'.format(len(data['dna_sub']), ', '.join(subs)))
        if 'known_dna_snp' in data:
            print('DBse snp: {0}'.format(', '.join(data['known_dna_snp'].split('|'))))
        print('')


def dna_data_to_dict(id, d, dna_data, item):
    pos = dna_data['p']
    strand = dna_data['f']
    base = dna_data['b']
    depth = dna_data['d']
    qual_base = dna_data['q']
    type_mut = dna_data['t']
    mut = ""

    if item in ['prot_snp', 'prot_sub']:
        mut = OrderedDict([('Feature', id), ('Mutation', '{0}{1}{2}'.format(d['t_aa'], d['t_prot_pos'], d['q_aa'])),
                           ('Mutation type', type_mut), ('DNA position', pos), ('DNA strand', strand), ('Codon', base),
                           ('Sequencing depth', depth), ('Base quality', qual_base)])
    elif item in ['dna_snp', 'dna_sub']:
        mut = OrderedDict([('Feature', id), ('Mutation', '{0}{1}{2}'.format(d['t_base'], d['t_dna_pos'], d['q_base'])),
                           ('Mutation type', type_mut), ('DNA position', pos), ('DNA strand', strand),
                           ('Codon/base', base),
                           ('Sequencing depth', depth), ('Base quality', qual_base)])
    return mut


def description(data):
    snps = []
    if 'prot_snp' in data:
        for d in data['prot_snp']:
            snps.append('{0}{1}{2}'.format(d['t_aa'], d['t_prot_pos'], d['q_aa']))
        snps = '|'.join(snps)
    elif 'dna_snp' in data:
        for d in data['dna_snp']:
            snps.append('{0}{1}{2}'.format(d['t_base'], d['t_dna_pos'], d['q_base']))
        snps = '|'.join(snps)

    subs = []
    if 'prot_sub' in data:
        for d in data['prot_sub']:
            subs.append('{0}{1}{2}'.format(d['t_aa'], d['t_prot_pos'], d['q_aa']))
        subs = '|'.join(subs)
    elif 'dna_sub' in data:
        for d in data['dna_sub']:
            subs.append('{0}{1}{2}'.format(d['t_base'], d['t_dna_pos'], d['q_base']))
        subs = '|'.join(subs)

    try:
        des = 'function:{0}, mechanism:{1}, reference_sequence: {2}, perc_identity:{3}, perc_coverage:{4}, ' \
              'min_depth:{5}, mean_depth:{6}, max_depth:{7}, ' \
              'min_quality:{8}, mean_quality:{9}, max_quality:{10}, ' \
              'known_sub:{11}'.format(data['func'], data['mech'], data['tid'].split('::')[-1], data['pid'], data['pcv'],
                                    data['min_depth'], data['mean_depth'], data['max_depth'],
                                    data['min_qual'], data['mean_qual'], data['max_qual'],
                                    snps)
    except KeyError:
        des = 'function:{0}, mechanism:{1}, reference_sequence: {2}, perc_identity:{3}, perc_coverage:{4}, ' \
              'known_sub:{5}'.format(data['func'], data['mech'], data['tid'].split('::')[-1], \
                                    data['pid'], data['pcv'], snps)
    return des


def write_fasta(results, out_dir, out_prefix):
    aa_outfile = os.path.join(out_dir, '{0}.faa'.format(out_prefix))
    dna_outfile = os.path.join(out_dir, '{0}.fna'.format(out_prefix))
    aa_records = []
    dna_records = []
    for data in results:
        id = '{0}__{1}__{2}__f{3}__{4}'.format(data['qid'], data['qstart'], data['qend'], data['strand'], data['tid'])
        des = description(data)

        if 'qprot' in data:
            q_prot_seq = str(data['qprot']).replace('-', '')
            aa_rec = SeqRecord(Seq(q_prot_seq), id=id, description=des)
            aa_records.append(aa_rec)

            if data['strand'] > 0:
                q_dna_seq = Seq(data['qseq'].replace('-', ''))
            else:
                q_dna_seq = Seq(data['qseq'].replace('-', '')).reverse_complement()
            dna_rec = SeqRecord(q_dna_seq, id=id, description=des)
            dna_records.append(dna_rec)
        else:
            dna_rec = SeqRecord(Seq(str(data['qseq']).replace('-', '')), id=id, description=des)
            dna_records.append(dna_rec)

    with open(aa_outfile, 'w') as aa_out_f, open(dna_outfile, 'w') as dna_out_f:
        if aa_records:
            aa_records = sorted(dna_records, key=lambda x: (int(x.id.split('_')[1]) * 1E8 + int(x.id.split('__')[1])))
            SeqIO.write(aa_records, aa_out_f, 'fasta')
        if dna_records:
            dna_records = sorted(dna_records, key=lambda x: (int(x.id.split('_')[1]) * 1E8 + int(x.id.split('__')[1])))
            SeqIO.write(dna_records, dna_out_f, 'fasta')


def write_gbk(results, query_dic, out_dir, out_prefix):
    rec_dic = {}
    for data in results:
        key = data['qid']
        if key not in rec_dic:
            rec_dic[key] = [data]
        else:
            rec_dic[key].append(data)
    # get keys in a list
    *keys, = rec_dic
    keys.sort()
    n = 0
    for key in keys:
        records = rec_dic[key]
        rec = SeqRecord(Seq(str(query_dic[key].seq), IUPACAmbiguousDNA()), id=key, name=key, description='')
        # source = ''
        for data in records:

            feature = SeqFeature(FeatureLocation(data['qstart'] - 1, data['qend'], strand=data['strand']),
                                 type='misc_feature', qualifiers={})

            if 'qprot' in data:
                # feature.qualifiers = {'locus_tag':'{}_{}'.format(outprefix, n), 'product':data['tid'].split('::')[0],
                #                  'note':description(data), 'translation':data['qprot']}
                feature.qualifiers = OrderedDict([('product', data['tid'].split('::')[0]),
                                                  ('note', description(data)), ('translation', data['qprot'])])
            else:
                # feature.qualifiers = {'locus_tag':'{}_{}'.format(outprefix, n), 'product':data['tid'].split('::')[0],
                #                   'note':description(data)}
                feature.qualifiers = OrderedDict([('product', data['tid'].split('::')[0]),
                                                  ('note', description(data))])
            rec.features.append(feature)

        rec.features = sorted(rec.features, key=lambda feature: feature.location.start)

        for feature in rec.features:
            feature.qualifiers = OrderedDict([('locus_tag', 'DET_{0}'.format(n + 1))] + list(feature.qualifiers.items()))
            n += 1

        out_file = os.path.join(out_dir, '{0}_{1}_mlst.gbk'.format(out_prefix, rec.id))
        with open(out_file, 'w') as out_f:
            SeqIO.write([rec], out_f, 'genbank')


def read_mlst_scheme(mlst_scheme_file, sep='\t', mlst_size=8):
    mlst_dic = {}
    mlst_present_list = []
    found_clonal_complex = False
    found_species = False
    with open(mlst_scheme_file) as in_f:
        for n, line in enumerate(in_f):
            if n == 0:
                mlst_present_list = line.strip().split(sep)[1:mlst_size + 1]
                #  remove alternative information
                if "clonal_complex" in mlst_present_list:
                    mlst_present_list.remove("clonal_complex")
                    found_clonal_complex = True
                if "species" in mlst_present_list:
                    mlst_present_list.remove("species")
                    found_species = True
            else:
                line = line.strip().split(sep)
                mlst_name = line[0]
                # treat mlst list with mlst schema name len
                list_mlst = line[1:len(mlst_present_list)+2]
                mlst_barcode = ' '.join(list_mlst)
                mlst_dic[mlst_barcode] = mlst_name

    return mlst_dic, mlst_present_list


def identify_mlst_profile(mlst_dic, mlst_list, blastn_results, id_prefix, out_prefix):
    mlst_barcode = []
    for item in mlst_list:
        pid = 0
        pcv = 0
        barcode = '0'
        found = 0
        for data in blastn_results:
            tid = data['tid']
            motif = re.compile('(^[A-Za-z0-9_]+)[-_.]([0-9]+$)')
            match = motif.match(tid)
            if match:
                gene, allele = match.groups()
                if gene.lower() == item.lower():
                    if data['pid'] == 100 and data['pcv'] == 100:
                        mlst_barcode.append(allele)
                        found = 1
                        break
                    else:
                        if data['pcv'] >= pcv and data['pid'] > pid:
                            pid = data['pid']
                            pcv = data['pcv']
                            barcode = allele
                        elif data['pcv'] > pcv:
                            pid = data['pid']
                            pcv = data['pcv']
                            barcode = allele
        if found == 0:
            mlst_barcode.append(barcode + '?')

    if ' '.join(mlst_barcode) in mlst_dic:
        ST = mlst_dic[' '.join(mlst_barcode)]
    else:
        ST = '?'

    zipped = zip(mlst_list, mlst_barcode)
    zip_list = list(map(list, zipped))

    result = OrderedDict([('MLST_name', out_prefix), ('ST', ST)] + zip_list)
    df = pd.DataFrame(result, index=[out_prefix])
    print(df)
    df.to_csv('{0}_{1}.csv'.format(id_prefix, out_prefix), sep='\t', index=False)
    df.to_html('{0}_{1}.html'.format(id_prefix, out_prefix))
    return blastn_results


def write_csv_html(results, mut_prefix, id_prefix, out_prefix, pass_alarm_qual=20, pass_alarm_depth=30,
                   pass_alarm_frac=0.9):
    header = ['Function', 'DBase name', '% ident', '% cov', 'Sequence Warning',
              'Min depth', 'Mean depth', 'Max depth',
              'Min qual', 'Mean qual', 'Max qual',
              'SNP', 'SUB', 'Known prot SNP', 'Known DNA SNP',
              'DBase start', 'DBase end', 'DBase length',
              'Query name', 'Query start', 'Query end', 'Query strand', 'Query length',
              'Query DNA seq', 'Query prot seq', 'DBase dna seq', 'DBase prot seq']
    records = []
    for data in results:
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
        subs = []
        if 'dna_sub' in data:
            for d in data['dna_sub']:
                alarm = []
                depth_alarm = ''
                qual_alarm = ''
                frac_alarm = ''
                if 'dna_data' in d:
                    if d['dna_data'] != '':
                        for dna_data in d['dna_data'].split(' | '):
                            dna_data = dict(x.split(':') for x in dna_data.strip().split('; '))
                            dna_data['t'] = 'Unknown'
                            mut = dna_data_to_dict(data['tid'], d, dna_data, 'dna_sub')
                            mutations.append(mut)

                            if depth_alarm == '':
                                depths = mut['Sequencing depth']
                                for depth in depths.split(','):
                                    depth_ref, depth_total = depth.split('/')
                                    try:
                                        if int(depth_ref) < pass_alarm_depth:
                                            depth_alarm = 'depth<{0}'.format(pass_alarm_depth)
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
                                        frac = int(depth_ref) / float(depth_total)
                                        if frac < pass_alarm_frac:
                                            frac_alarm = 'frac<{0}'.format(pass_alarm_frac)
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
                                            qual_alarm = 'qual<{0}'.format(pass_alarm_qual)
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

                if alarm == '':
                    subs.append('{0}{1}{2}'.format(d['t_base'], d['t_dna_pos'], d['q_base']))
                else:
                    subs.append('{0}{1}{2} #{3}#'.format(d['t_base'], d['t_dna_pos'], d['q_base'], alarm))

                if mutations:
                    df = pd.DataFrame.from_dict(mutations)
                    df.to_html('{0}_{1}_{2}_{3}-{4}_mut.html'.format(mut_prefix,
                                                                     data['tid'].replace(':', '_').replace('(', '').replace(
                                                                    ')', '').replace('\'', 'pr'), data['qid'],
                                                                     data['qstart'], data['qend']), index=False)
                    df.to_csv('{0}_{1}_{2}_{3}-{4}_mut.csv'.format(mut_prefix,
                                                                   data['tid'].replace(':', '_').replace('(', '').replace(
                                                                  ')', '').replace('\'', 'pr'), data['qid'],
                                                                   data['qstart'], data['qend']), sep='\t', index=False)

        subs = ', '.join(subs)

        values = values + ['', subs, '', '']
        values = values + [data['tstart'], data['tend'], data['tlen']]
        values = values + [data['qid'], data['qstart'], data['qend'], data['strand'], data['qend'] - data['qstart'] + 1]

        if data['strand'] > 0:
            q_dna_seq = str(data['qseq'])
        else:
            q_dna_seq = str(Seq(data['qseq']).reverse_complement())

        if 'tseq' in data:
            t_dna_seq = str(data['tseq'])
        else:
            t_dna_seq = ''

        values = values + [q_dna_seq, '', t_dna_seq, '']
        d = dict(zip(header, values))
        records.append(d)
    df = pd.DataFrame.from_dict(records)
    df.sort_values(['Function', 'DBase name'], ascending=[True, True], inplace=True)
    df[header].to_csv('{0}_{1}_details.csv'.format(id_prefix, out_prefix), sep='\t', index=False)
    df[header].to_html('{0}_{1}_details.html'.format(id_prefix, out_prefix), index=False)


def main(args):

    # Check working directory
    wk_dir = args.wk_dir
    if glob.glob(wk_dir) is None:
        print("\nDirectory {0} not found!\n".format(wk_dir))
        exit(1)
    else:
        print("\nDirectory {0} found".format(wk_dir))

    # Prefix
    out_prefix = args.out_prefix

    # Assign sample file
    sample_file = args.sample_file
    if sample_file is None:
        sample_file = os.path.join(os.path.dirname(wk_dir), 'sample.csv')

    # Check sample file
    if not os.path.exists(sample_file):
        print("\nSample file {0} not found!\n".format(sample_file))
        exit(1)
    else:
        print("\nSample file: {0}".format(sample_file))

    # Load sample name and species if provided (not required)
    sample_dic = load_sample_name(sample_file)

    # Search the assembly files corresponding to the sample names in the working directory
    # (e.g. ./working_directory/*/sample_name.fasta)
    print("\nSearch the assembly files:")
    file_list = []
    for sample_name in sample_dic.keys():
        try:
            sample_file = glob.glob(os.path.join(wk_dir, sample_name + '.fasta'))[0]
            file_list.append(sample_file)
            print('Sample {0}: {1}'.format(sample_name, sample_file))
        except IndexError:
            print('Sample {0}: file not found!'.format(sample_name))

    # Check the location of MLST file:
    set_file = args.set_file
    set_dic = {}
    db_dir = ""
    if os.path.exists(set_file):
        print("\nSetting file: {0}".format(set_file))
        set_dic = load_set_file(set_file)
        db_dir = os.path.dirname(set_file)
    else:
        print('\nSetting file {0} not found!\n'.format(set_file))
        exit(1)

    print("\nFeature detection parameters:")
    pass_pid = float(args.perc_cv)
    print("  Identity threshold: {0}%".format(pass_pid))
    pass_pcv = float(args.perc_id)
    print("  Coverage threshold: {0}%".format(pass_pcv))
    pass_overlap = int(args.overlap)
    print("  Maximum feature overlap: {0}-bases\n".format(pass_overlap))

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
        taxonomy = sample_dic[sample_name].lower()
        title = "~~~~  {0}/{1} sample: {2} species: {3}  ~~~~".format(n + 1, len(file_list), sample_name, taxonomy)
        print("\n\n")
        print("~" * len(title))
        print(title)
        print("~" * len(title))

        # Set MLST scheme:
        mlst_db = ""
        if 'mlst' in set_dic[taxonomy]:
            mlst_db = set_dic[taxonomy]['mlst']

        else:
            print("No schema MLST for {0}", taxonomy)
            exit(1)

        mlst_schema = []
        if "&" in mlst_db[1]:
            for name_mlst in mlst_db[1].split("&"):
                mlst_schema.append( name_mlst)
        else:
            mlst_schema = [mlst_db[1]]

        for schema in mlst_schema:
            print("Current schema: {0}".format(schema))
            mlst_name = '{0}_{1}'.format(mlst_db[0], schema)

            mlst_dir = os.path.join(db_dir, "dbMLST", mlst_name, 'pubmlst_download')
            print("MLST directory selected : {0}".format(mlst_dir))
            mlst_scheme_file = os.path.join(mlst_dir, 'profile.txt')
            dna_target_file = os.path.join(wk_dir, 'profile_{0}.fasta'.format(schema))

            if not os.path.exists(mlst_scheme_file):
                print("\nNo MLST scheme defined for {0} in {1}\n".format(taxonomy, set_file))
                dna_target_file, mlst_scheme_file = 0, 0
                exit(1)
            if not os.path.exists(dna_target_file):
                print("\nNo fasta file for {0}\n".format(taxonomy))
                cmd = 'cat {0}/*.tfa > {1}/profile_{2}.fasta'.format(mlst_dir, wk_dir, schema)
                print("Fasta file created ! ")
                print(cmd)
                os.system(cmd)

            # Make blastn database, launch blastn and load the results
            blastn_results = []

            out_blastn_file = os.path.join(os.path.dirname(query_file), 'blastn_output_{0}.csv'.format(schema))

            if os.path.exists(dna_target_file):
                blastn_db = make_blastn_database(dna_target_file, force)
                blastn_result_file = run_blastn(blastn_db, query_file, pass_pid, force, 0.0001, 8, out_blastn_file)

                blastn_results = load_blastn_result(blastn_result_file, dna_target_file, pass_pid, pass_pcv)
                print('\nNumber of detected features: {0}'.format(len(blastn_results)))

                for element in blastn_results:
                    if "uidA" in element:
                        print(element)

                # Filter the results for overlaps
                blastn_results = overlap_filter(blastn_results, pass_overlap)
                print('Number of detected features after overlap filtering: {0}'.format(len(blastn_results)))
            else:
                print("\nMLST database file {0} not found!\n".format(dna_target_file))
                exit(1)

            query_dic = {}
            if len(blastn_results) > 0:
                query_dic = load_fasta(query_file)
                # Set the prefix of the output
                print('')
                if args.out_prefix == '':
                    out_prefix = mlst_name
                else:
                    out_prefix = args.out_prefix
                # Global alignement of DNA and mutation extraction if DNA features detected
                target_dic = load_fasta(dna_target_file)
                blastn_results = dna_global_alignemnt(blastn_results, query_dic, target_dic, pass_pid, pass_pcv)
                if os.path.exists(bam_file):
                    # Extaction quality of bases and sequencing depth if bam detected
                    blastn_results = dna_extract_quality_and_depth(bam_file, query_file, blastn_results, out_prefix, force)
                # Show the detected DNA features
                view_dna_result(blastn_results)

            # Set the directory to store depth and quality data
            mut_dir = os.path.join(os.path.dirname(bam_file), 'Mutations_depth_quality_{0}'.format(schema))
            if not os.path.exists(mut_dir):
                os.mkdir(mut_dir)
            mut_prefix = os.path.join(mut_dir, os.path.splitext(os.path.basename(query_file))[0].split('_')[0])
            mlst_id_prefix = os.path.splitext(query_file)[0]

            if os.path.exists(mlst_scheme_file):
                # Load MLST profils
                print('MLST name: {0}'.format(mlst_name))
                mlst_dic, mlst_present_list = read_mlst_scheme(mlst_scheme_file)

                # Identify ST from MLST profils
                print('Number of MLST profiles: {0} for genes: {1}\n'.format(len(mlst_dic.keys()), ', '.join(mlst_present_list)))
                identify_mlst_profile(mlst_dic, mlst_present_list, blastn_results, mlst_id_prefix, out_prefix)

            if len(blastn_results) > 0:
                write_csv_html(blastn_results, mut_prefix, mlst_id_prefix, out_prefix)
                out_dir = os.path.join(os.path.dirname(query_file), 'Detected_sequences_{0}'.format(schema))
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)
                write_gbk(blastn_results, query_dic, out_dir, out_prefix)
                write_fasta(blastn_results, out_dir, out_prefix)
            else:
                print('\nNo results!\n')


def version():
    return "1.0"


def run():
    global usage

    usage = "mlstDetector.py [-st setting database ] [-wd work directory with the assembly] [-sf sample file] "

    parser = argparse.ArgumentParser(
        prog='mlstDetector',
        usage=usage,
        description='MlstDetector: pipeline CNR Resistance for MLST Detection - Version ' + version(),
    )
    parser.add_argument('-sf', '--sampleFile', dest="sample_file",
                        help='Tab-separated file containing the names of sample and the corresponding taxonomy (species/genus) [$wkdir/../sample.csv]')
    parser.add_argument('-wd', '--wkdir', dest="wk_dir", help='Directory containing the assembly files')
    parser.add_argument('-st', '--setFile', dest="set_file", default="/usr/local/readmapper-v0.1/setting.txt",
                        help='Setting file')
    parser.add_argument('-id', '--perc_id', dest="perc_id", default="98",
                        help="Minimum identity percentage with the target [98]")
    parser.add_argument('-cv', '--perc_cv', dest="perc_cv", default="98",
                        help="Minimum target coverage percentage [98]")
    parser.add_argument('-ov', '--overlap', dest="overlap", default="20",
                        help="Maximun overlap between detected features as base number [20]")
    parser.add_argument('-th', '--threads', dest="threads", default='8', help='Number of threads [8]')
    parser.add_argument('-F', '--Force', dest="Force", action="store_true", default=False,
                        help="Overwrite the detected files")
    parser.add_argument('-o', '--outPrefix', dest="out_prefix", default='', help="Outprefix [<database_name>]")
    parser.add_argument('-V', '--version', action='version', version='diamDetector-' + version(),
                        help="Prints version number")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)


    main(args)


if __name__ == '__main__':
    run()
