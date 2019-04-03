#!/usr/bin/python3
import csv
import os
import re
import subprocess
from collections import OrderedDict
from shutil import rmtree

from Bio import SeqIO
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_fasta(fasta_file, fmt='fasta'):
    fasta_dic = {}
    with open(fasta_file) as fasta:
        for rec in SeqIO.parse(fasta, fmt):
            rec.description = rec.description.replace(rec.id, '').strip()
            fasta_dic[rec.id] = rec
    return fasta_dic


def load_sample_name(sample_file):
    sample_dic = OrderedDict()
    with open(sample_file) as inp_f:
        for line in inp_f:
            line = line.strip().split('\t')
            try:
                sample_dic[line[0]] = line[1]
            except IndexError:
                sample_dic[line[0]] = 'No taxonomy provided'
    return sample_dic


def evaluate_dna_alg(seq, ref_seq):
    nid = 0
    gap = 0
    gap_tag = 'off'
    open_gap = 0
    a_len = 0
    for i in range(0, len(ref_seq)):
        base, ref = seq[i], ref_seq[i]
        a_len += 1
        if base == '-' or ref == '-':
            gap += 1
            gap_tag = 'on'
            if gap == 'off':
                open_gap += 1
        elif base == ref:
            nid += 1
            gap_tag = 'off'
        else:
            gap_tag = 'off'
    pid = round(100 * nid / float(a_len), 2)
    pcv = round(100 * (len(seq.replace('-', ''))) / float(len(ref_seq.replace('-', ''))), 2)
    return {'nid': nid, 'pid': pid, 'pcv': pcv, 'pos': nid, 'ppos': pid, 'gap': gap, 'opengap': open_gap, 'alen': a_len}


def read_bam_stat(filename):
    with open(filename) as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        for row in reader:
            if row["ID"] == "Overall":
                return row


def make_dmnd_database(fasta_file, force=False, threads=8):
    dmnd_db = os.path.splitext(fasta_file)[0]
    if os.path.exists(dmnd_db + ".dmnd") and not force:
        print('\nDatabase {0} already exists'.format(dmnd_db))
    else:
        process = subprocess.Popen('$(which diamond) makedb -p {0} --in {1} -d {2}'
                                   .format(threads, fasta_file, dmnd_db),
                                   shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
        print("\n{0}\n".format(process.decode("utf-8")))

    return dmnd_db


def cds_global_alignment(dmnd_results, query_dic, wk_dir):
    for data in dmnd_results:
        q_seq = query_dic[data['qid']]
        # print '\n', data['tid'], data['strand'], data['tlen'], data['tstart'], data['tend']
        if data['strand'] > 0:
            m = data['qend']
            while str(q_seq[data['qend']:m + 3].seq)[-3:].upper() not in ['TAA', 'TGA', 'TAG'] and m + 3 < len(q_seq):
                m = m + 3
            if m + 3 < len(q_seq):
                data['qend'] = m + 3
            data['qseq'] = str(q_seq[data['qstart'] - 1:data['qend']].seq)
        else:
            m = data['qstart'] - 1
            while str(q_seq[m - 3:data['qstart']].seq[:3]).upper() not in ['TTA', 'TCA', 'CTA'] and m - 3 >= 0:
                m = m - 3
            if m - 3 >= 0:
                data['qstart'] = m - 2
            data['qseq'] = str(q_seq[data['qstart'] - 1:data['qend']].seq)
        if data['strand'] > 0:
            if str(data['qseq'][0:3]) not in ['ATG', 'GTG', 'TTG', 'ATT', 'CTG'] or data['tstart'] > 1:
                pos = []
                scores = []
                m = data['qstart'] + 17
                while m - 3 >= 0:
                    m = m - 3
                    seq = str(q_seq[m:data['qend']].seq).upper()
                    l_seq = len(seq)
                    codon = seq[:3]
                    if codon in ['TAA', 'TGA']:
                        break
                    else:
                        if codon in ['ATG', 'GTG', 'TTG', 'ATT', 'CTG']:
                            score = abs((3 * data['tlen']) - l_seq)
                            if not scores:
                                pos.append(m)
                                scores.append(score)
                            elif score < min(scores):
                                pos.append(m)
                                scores.append(score)
                            else:
                                break
                if scores:
                    m = pos[scores.index(min(scores))]
                    data['qstart'] = m + 1
                    # print data['qseq']
                    data['qseq'] = str(q_seq[data['qstart'] - 1:data['qend']].seq)
                    # print data['qseq']
        else:
            if str(data['qseq'][-3:]) not in ['CAT', 'CAC', 'CAA'] or data['tstart'] > 1:
                pos = []
                scores = []
                m = data['qend']
                while m + 3 < len(q_seq.seq):
                    m = m + 3
                    seq = str(q_seq[data['qstart'] - 1:m].seq).upper()
                    l_seq = len(seq)
                    codon = seq[-3:]
                    if codon in ['TTA', 'TCA']:
                        break
                    else:
                        if codon in ['CAT', 'CAC', 'CAA']:
                            score = abs((3 * data['tlen']) - l_seq)
                            if not scores:
                                pos.append(m)
                                scores.append(score)
                            elif score < min(scores):
                                pos.append(m)
                                scores.append(score)
                            else:
                                break
                if scores:
                    m = pos[scores.index(min(scores))]
                    data['qend'] = m
                    # print data['qseq']
                    data['qseq'] = str(q_seq[data['qstart'] - 1:data['qend']].seq)
                    # print data['qseq']
        data = extract_substitutions(data, wk_dir)
    return dmnd_results


def blastn_to_global_alignment(blastn_results, query_dic, target_dic):
    for data in blastn_results:
        tseq = data['tseq']
        qseq = data['qseq']
        full_t_seq = target_dic[data['tid']].seq
        if str(full_t_seq) != str(tseq).replace('-', ''):
            full_qseq = query_dic[data['qid']].seq
            if data['strand'] > 0:
                q_start = data['qstart'] - (data['tstart'] - 1)
                q_end = data['qend'] + (data['tlen'] - data['tend'])
                if q_start < 1:
                    data['qseq'] = '-' * (abs(q_start) + 1) + full_qseq[0:data['qstart'] - 1] + data['qseq']
                    data['qstart'] = 1
                else:
                    data['qseq'] = full_qseq[q_start - 1:data['qstart'] - 1] + data['qseq']
                    data['qstart'] = q_start

                if q_end > len(full_qseq):
                    data['qseq'] = data['qseq'] + full_qseq[data['qend']:] + '-' * (q_end - len(full_qseq))
                    data['qend'] = len(full_qseq)
                else:
                    data['qseq'] = data['qseq'] + full_qseq[data['qend']:q_end]
                    data['qend'] = q_end
            else:
                q_start = data['qstart'] - (data['tlen'] - data['tend'])
                q_end = data['qend'] + data['tstart'] - 1
                if q_start < 1:
                    data['qseq'] = data['qseq'] + reverse_complement(full_qseq[0:data['qstart'] - 1]) + '-' * (
                            abs(q_start) + 1)
                    data['qstart'] = 1
                else:
                    data['qseq'] = data['qseq'] + reverse_complement(full_qseq[q_start - 1:data['qstart'] - 1])
                    data['qstart'] = q_start

                if q_end > len(full_qseq):
                    data['qseq'] = '-' * (q_end - len(full_qseq)) + reverse_complement(full_qseq[data['qend']:]) + data[
                        'qseq']
                    data['qend'] = len(full_qseq)
                else:
                    data['qseq'] = reverse_complement(full_qseq[data['qend']:q_end]) + data['qseq']
                    data['qend'] = q_end
            data['tseq'] = full_t_seq[0:data['tstart'] - 1] + data['tseq'] + full_t_seq[data['tend']:]
            data['tstart'] = 1
            data['tend'] = len(full_t_seq)
            d = evaluate_dna_alg(str(data['qseq']), str(data['tseq']))
            data['pid'] = d['pid']
            data['nid'] = d['nid']
            data['pos'] = d['pos']
            data['ppos'] = d['ppos']
            data['alen'] = d['alen']
            data['gap'] = d['gap']
            data['opengap'] = d['opengap']
            data['pcv'] = d['pcv']
        data = extract_mutations(data)
    return blastn_results


def extract_substitutions(data, wk_dir):
    # translate query sequence
    q_seq = data['qseq']
    t_seq = data['fulltseq']
    strand = data['strand']
    known_snps = data['known_prot_snp']
    if strand > 0:
        try:
            q_seq = Seq(q_seq).translate(table='Bacterial', cds=True)
        except Exception as e:
            print(e)
            if len([x for x in list(set(q_seq)) if x.upper() not in ['A', 'T', 'C', 'G']]) > 0:
                try:
                    q_seq = Seq(q_seq, IUPACAmbiguousDNA()).translate(table='Bacterial', cds=True)
                except Exception as e:
                    print(e)
                    q_seq = Seq(q_seq, IUPACAmbiguousDNA()).translate(table='Bacterial', cds=False)
            else:
                q_seq = Seq(q_seq, IUPACAmbiguousDNA()).translate(table='Bacterial', cds=False)
    else:
        try:
            q_seq = Seq(q_seq).reverse_complement().translate(table='Bacterial', cds=True)
        except Exception as e:
            print(e)
            if len([x for x in list(set(q_seq)) if x.upper() not in ['A', 'T', 'C', 'G']]) > 0:
                try:
                    q_seq = Seq(q_seq, IUPACAmbiguousDNA()).reverse_complement().translate(table='Bacterial', cds=True)
                except Exception as e:
                    print(e)
                    q_seq = Seq(q_seq, IUPACAmbiguousDNA()).reverse_complement().translate(table='Bacterial', cds=False)
            else:
                q_seq = Seq(q_seq).reverse_complement().translate(table='Bacterial', cds=False)

    # Convert target sequence as Seq object
    try:
        t_seq = Seq(t_seq)
    except Exception as e:
        print(e)
        t_seq = Seq(t_seq, ExtendedIUPACProtein)

    # Generate alignment with clustalo
    records = []
    for ID, seq in [('qseq', q_seq), ('tseq', t_seq)]:
        rec = SeqRecord(seq, id=ID, description='')
        records.append(rec)

    tmp_clustalo_path = os.path.join(wk_dir, "tmp_clustalo")

    os.makedirs(tmp_clustalo_path)

    in_file = os.path.join(tmp_clustalo_path, 'tmp.fasta')
    out_file = os.path.join(tmp_clustalo_path,'tmp.clu')
    with open(in_file, 'w') as in_f:
        SeqIO.write(records, in_f, 'fasta')
    cmd = "$(which clustalo) -i {0} -o {1} --outfmt=fa --force".format(in_file, out_file)

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    print("\n{0}\n{1}\n".format(cmd, process.decode("utf-8")))


    # parse alignment
    records = []
    with open(out_file) as inf_f:
        for rec in SeqIO.parse(inf_f, 'fasta'):
            records.append(rec)

    # remove clustalo work dir
    rmtree(tmp_clustalo_path)

    q_seq = records[0]
    t_seq = records[1]
    gaps = q_seq.seq.count('-') + t_seq.seq.count('-')
    open_gaps = 0
    for seq in [q_seq, t_seq]:
        seq = str(seq.seq)
        while seq.count('--') > 0:
            seq = seq.replace('--', '-')
        open_gaps = open_gaps + seq.count('-')
    dif_indexes = [i for i in range(0, len(t_seq.seq)) if t_seq.seq[i] != q_seq.seq[i]]
    dif_indexes.sort()
    mismatch = len(dif_indexes)
    nid = len(t_seq.seq) - mismatch
    pid = round(100 * nid / float(len(t_seq.seq)), 2)
    pcv = round(100 * len(str(q_seq.seq).replace('-', '')) / float(len(str(t_seq.seq).replace('-', ''))), 2)

    # extract mutations from alignment
    unclassified_subs = []
    known_snps = known_snps.split(',')
    known_pos = list(set([int(x[1:len(x) - 1]) for x in known_snps if re.match('[A-Zi][0-9]+[A-Zd]', x)]))
    indel = 0
    for i in dif_indexes:
        t_prot_pos = i + 1 - t_seq.seq[0:i].count('-')
        q_prot_pos = i + 1 - q_seq.seq[0:i].count('-')
        a_prot_pos = i + 1
        t_aa = t_seq.seq[i]
        q_aa = q_seq.seq[i]
        if t_aa == '-':
            t_aa = 'i'
            indel = 1
        if q_aa == '-':
            q_aa = 'd'
            indel = 1
        unclassified_subs.append(
            {'t_prot_pos': t_prot_pos, 'q_prot_pos': q_prot_pos, 'a_prot_pos': a_prot_pos, 't_aa': t_aa, 'q_aa': q_aa})

    # format indels
    if indel == 1:
        pre_pos, pre_query, pre_target = 0, '', ''
        index_dels = []
        for n, d in enumerate(unclassified_subs):
            target = d['t_aa']
            pos = d['a_prot_pos']
            query = d['q_aa']

            if query == 'd' and pre_query == 'd' and pos == pre_pos + 1:
                unclassified_subs[open_deletion_index]['t_aa'] = unclassified_subs[open_deletion_index]['t_aa'] + \
                                                                 unclassified_subs[n]['t_aa']
                index_dels.append(n)
            elif query == 'd':
                open_deletion_index = n

            if target == 'i' and pre_target == 'i' and pos == pre_pos + 1:
                unclassified_subs[open_insertion_index]['q_aa'] = unclassified_subs[open_insertion_index]['q_aa'] + \
                                                                  unclassified_subs[n]['q_aa']
                index_dels.append(n)
            elif target == 'i':
                open_insertion_index = n

            pre_target = target
            pre_pos = pos
            pre_query = query

        index_dels.reverse()
        for i in index_dels:
            del unclassified_subs[i]

    # classify the mutation as known (snps) or unknown (subs)
    subs, snps = [], []
    if known_snps != ['']:
        for m in unclassified_subs:
            if m['t_prot_pos'] in known_pos:
                item = '{0}{1}{2}'.format(m['t_aa'], m['t_prot_pos'], m['q_aa'])
                if item in known_snps:
                    snps.append(m)
                else:
                    m['q_aa'] = '[{0}]'.format(m['q_aa'])
                    snps.append(m)
            else:
                subs.append(m)
    else:
        subs = unclassified_subs

    # update data before to return
    data['pid'] = pid
    data['pcv'] = pcv
    data['nid'] = nid
    data['mismatch'] = mismatch
    data['gaps'] = gaps
    data['opengaps'] = open_gaps
    data['qprot'] = q_seq.seq
    data['tprot'] = t_seq.seq
    data['prot_sub'] = subs
    data['prot_snp'] = snps
    return data


def reverse_complement(dna_seq):
    cpl_dic = {'-': '-', 'N': 'N',
               'A': 'T', 'T': 'A',
               'C': 'G', 'G': 'C',
               'S': 'S', 'W': 'W',
               'D': 'H', 'H': 'D',
               'K': 'M', 'M': 'K',
               'Y': 'R', 'R': 'Y',
               'B': 'V', 'V': 'B'}
    dna_seq = list(str(dna_seq).upper())
    size = len(dna_seq)
    rev_cplt = []
    while len(rev_cplt) < size:
        rev_cplt.append(cpl_dic[dna_seq[-1]])
        dna_seq = dna_seq[:-1]
    return ''.join(rev_cplt)


def dna_global_alignemnt(blastn_results, query_dic, target_dic):
    for data in blastn_results:
        tseq = data['tseq']
        qseq = data['qseq']
        full_tseq = target_dic[data['tid']].seq
        if str(full_tseq) != str(tseq).replace('-', ''):
            full_qseq = query_dic[data['qid']].seq
            if data['strand'] > 0:
                qstart = data['qstart'] - (data['tstart'] - 1)
                qend = data['qend'] + (data['tlen'] - data['tend'])
                if qstart < 1:
                    data['qseq'] = '-' * (abs(qstart) + 1) + full_qseq[0:data['qstart'] - 1] + data['qseq']
                    data['qstart'] = 1
                else:
                    data['qseq'] = full_qseq[qstart - 1:data['qstart'] - 1] + data['qseq']
                    data['qstart'] = qstart

                if qend > len(full_qseq):
                    data['qseq'] = data['qseq'] + full_qseq[data['qend']:] + '-' * (qend - len(full_qseq))
                    data['qend'] = len(full_qseq)
                else:
                    data['qseq'] = data['qseq'] + full_qseq[data['qend']:qend]
                    data['qend'] = qend
            else:
                qstart = data['qstart'] - (data['tlen'] - data['tend'])
                qend = data['qend'] + data['tstart'] - 1
                if qstart < 1:
                    data['qseq'] = data['qseq'] + reverse_complement(full_qseq[0:data['qstart'] - 1]) + '-' * (
                                abs(qstart) + 1)
                    data['qstart'] = 1
                else:
                    data['qseq'] = data['qseq'] + reverse_complement(full_qseq[qstart - 1:data['qstart'] - 1])
                    data['qstart'] = qstart

                if qend > len(full_qseq):
                    data['qseq'] = '-' * (qend - len(full_qseq)) + reverse_complement(full_qseq[data['qend']:]) + data[
                        'qseq']
                    data['qend'] = len(full_qseq)
                else:
                    data['qseq'] = reverse_complement(full_qseq[data['qend']:qend]) + data['qseq']
                    data['qend'] = qend
            data['tseq'] = full_tseq[0:data['tstart'] - 1] + data['tseq'] + full_tseq[data['tend']:]
            data['tstart'] = 1
            data['tend'] = len(full_tseq)
            d = evaluate_dna_alg(str(data['qseq']), str(data['tseq']))
            data['pid'] = d['pid']
            data['nid'] = d['nid']
            data['pos'] = d['pos']
            data['ppos'] = d['ppos']
            data['alen'] = d['alen']
            data['gap'] = d['gap']
            data['opengap'] = d['opengap']
            data['pcv'] = d['pcv']
        data = extract_mutations(data)
    return blastn_results


def extract_mutations(data):
    qseq = data['qseq']
    tseq = data['tseq']
    strand = data['strand']
    known_snps = data['known_dna_snp']

    dif_indexes = [i for i in range(0, len(tseq)) if str(tseq)[i] != str(qseq)[i]]
    dif_indexes.sort()
    # extract mutations from alignement and depth from bam
    muts = []
    known_snps = known_snps.split(',')
    known_pos = list(set([int(x[1:len(x) - 1]) for x in known_snps if re.match('[A-Zi][0-9]+[A-Zd]', x)]))
    indel = 0
    for i in dif_indexes:
        t_dna_pos = i + 1 - str(tseq)[0:i].count('-')
        q_dna_pos = i + 1 - str(qseq)[0:i].count('-')
        a_dna_pos = i + 1
        t_base = str(tseq)[i]
        q_base = str(qseq)[i]
        if t_base == '-':
            t_base = 'i'
            indel = 1
        if q_base == '-':
            q_base = 'd'
            indel = 1
        muts.append({'t_dna_pos': t_dna_pos, 'q_dna_pos': q_dna_pos, 'a_dna_pos': a_dna_pos, 't_base': t_base,
                     'q_base': q_base})

    # format indels
    if indel == 1:
        pre_pos, pre_query, pre_target = 0, '', ''
        index_dels = []
        open_insertion_index = ""
        open_deletion_index = ""
        for n, d in enumerate(muts):
            target = d['t_base']
            pos = d['a_dna_pos']
            query = d['q_base']

            if query == 'd' and pre_query == 'd' and pos == pre_pos + 1:
                muts[open_deletion_index]['t_base'] = muts[open_deletion_index]['t_base'] + muts[n]['t_base']
                index_dels.append(n)
            elif query == 'd':
                open_deletion_index = n

            if target == 'i' and pre_target == 'i' and pos == pre_pos + 1:
                muts[open_insertion_index]['q_base'] = muts[open_insertion_index]['q_base'] + muts[n]['q_base']
                index_dels.append(n)
            elif target == 'i':
                open_insertion_index = n

            pre_target = target
            pre_pos = pos
            pre_query = query

        index_dels.reverse()
        for i in index_dels:
            del muts[i]

    # classify the mutation as known (snps) or unknown (subs)
    subs, snps = [], []
    if known_snps != ['']:
        for m in muts:
            if m['t_dna_pos'] in known_pos:
                item = '{0}{1}{2}'.format(m['t_base'], m['t_dna_pos'], m['q_base'])
                if item in known_snps:
                    snps.append(m)
                else:
                    m['q_base'] = '[{0}]'.format(m['q_base'])
                    snps.append(m)
            else:
                subs.append(m)
    else:
        subs = muts

    # update data before to return
    data['dna_sub'] = subs
    data['dna_snp'] = snps
    return data


def read_bam_count(filename, depth_pass=20, qual_pass=20, fraction_pass=0.8):
    with open(filename) as tsv_file:

        reader = csv.DictReader(tsv_file, delimiter='\t')
        result = {}
        alarm = []
        # for n, line in enumerate(inf_f):
        for row in reader:
            # line = line.strip().split('\t')
            # if n == 0:
            #    header = line
            # else:

            ctg = row['ID']
            position = row['position']
            ref = row['reference'].upper()
            total_depth = int(row['total_depth'])
            cmt = []
            if ref in ['A', 'T', 'C', 'G']:
                base_depth = int(row['{0}_depth'.format(ref)])
                allele_fraction = round(base_depth / float(total_depth), 2)
                base_qual = float(row['{0}_quality'.format(ref)])
                if base_depth < depth_pass:
                    cmt.append('D{0}={1}'.format(position, base_depth))
                if allele_fraction < fraction_pass:
                    cmt.append('F{0}={1}'.format(position, allele_fraction))
                if base_qual < qual_pass:
                    cmt.append('Q{0}={1}'.format(position, base_qual))
            else:
                cmt.append('{0}{1}'.format(ref, position))

            cmt = ','.join(cmt)
            if cmt != '':
                alarm.append(cmt)

            if ctg not in result:
                result[ctg] = {position: row}
            else:
                result[ctg][position] = row
    return result, ';'.join(alarm)
