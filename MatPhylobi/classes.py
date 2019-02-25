"""

"""
import datetime
import sys
import os
import csv
import itertools
import logging
from collections import OrderedDict
from operator import itemgetter
from Bio import Entrez, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from .utils import check, makeconfig, readconfig, dameraulevenshtein, clean_up

logger = logging.getLogger('MatPhylobi')
logging.basicConfig(stream=sys.stdout, level=logging.INFO)


class Sequence:
    def __init__(self, marker, accession, strand='', segment=None):
        self.accession = accession
        self.marker = marker
        self.seq = ''
        self.length = 0
        self.species = self.marker.organism
        self.voucher = ''
        self.tax = ''
        if segment and len(set(segment)) == 1:
            self.segment = segment[0]
        else:
            self.segment = None
        self.author = ''
        self.num = 0
        self.sumid = 0
        self.sumqc = 0
        self.strand = strand
        self.isolate = ''
        self.Ncount = 0

    def fasta(self, acc=False):
        if not self.seq:
            self.assign_info(subsp=self.marker.run.subsp)
        if acc:
            return ">" + self.accession + "\n" + str(self.seq) + "\n"
        return ">" + self.species + "\n" + str(self.seq) + "\n"

    def download_gb_record(self):
        record = ''
        n = 0

        while not record:
            try:
                handle = Entrez.efetch(db="nucleotide", id=self.accession, retmode='text', rettype="gb")
                record = SeqIO.read(handle, "gb")
                handle.close()
            except:
                n += 1
                if n > 10:
                    logger.critical("I cannot find this acc number: " + str(self.accession) + ". Please check your "
                                                                                              "internet onnection.")
                    clean_up(self.marker.run.debug, self.marker.run.folder, self.marker.run.opt, 1)
        return record

    def assign_info(self, its=False, subsp=False):

        def assign_specimen_specifics(features):
            for feat in features:
                if feat.type.lower() == 'source':
                    if "specimen_voucher" in feat.qualifiers.keys():
                        self.voucher = str(feat.qualifiers["specimen_voucher"][0])
                    if "isolate" in feat.qualifiers.keys():
                        self.isolate = str(feat.qualifiers["isolate"])
                    if "db_xref" in feat.qualifiers.keys():
                        for dbx in feat.qualifiers["db_xref"]:
                            if dbx.startswith('taxon:'):
                                self.tax = dbx[6:]

        def check_synonyms():
            if self.tax in self.marker.run.synonyms:
                if self.marker.run.synonyms[self.tax].isnumeric():
                    self.tax = self.marker.run.synonyms[self.tax]
                    handle = Entrez.efetch(db="taxonomy", id=self.tax, retmode="xml")
                    tax_record = Entrez.read(handle)
                    self.species = "_".join(tax_record[0]["ScientificName"].split()).lower()
                else:
                    self.species = self.marker.run.synonyms[self.tax]
                    self.tax = self.marker.run.synonyms[self.tax]

        def prepare_species(annotations, subsp):
            if self.species == self.marker.organism:
                self.species = "_".join(annotations["organism"].split()).lower()
            if not subsp:
                self.species = "_".join(self.species.split('_')[:2])
            if "_sp." in self.species or "_aff." in self.species or self.species.startswith('x_'):
                return False
            self.species = self.species.strip("'")
            self.species = self.species[0].upper() + self.species[1:]
            self.tax = "_".join(self.tax.split())
            return True

        def get_its_segment(annotations):
            for ref in annotations["references"]:
                self.author = "_".join(ref.authors.split()).lower()
                if "segment" in annotations.keys():
                    self.segment = str(annotations["segment"][:1])

        def check_refseq(annotations):
            if "references" in annotations.keys():
                if its and self.marker.name == 'ITS':
                    get_its_segment(annotations)
                return True
            else:
                # probably a RefSeq sequence
                return False

        def prepare_sequence(seq):
            if self.strand == 'minus':
                self.seq = seq.reverse_complement()
                logger.debug('reversed ' + self.accession)
            else:
                self.seq = seq
            self.Ncount = self.seq.count('N') + self.seq.count('n')
            self.length = len(seq)

        record = self.download_gb_record()

        prepare_sequence(record.seq)

        assign_specimen_specifics(record.features)
        check_synonyms()
        prepared = prepare_species(record.annotations, subsp)

        not_refseq = check_refseq(record.annotations)

        return prepared & not_refseq

    def join(self, other):
        if self.segment == "1" and other.segment == "2":
            self.accession = '(' + self.accession + ',' + other.accession + ')'
            self.seq = self.seq + 100 * 'N' + other.seq
        else:
            self.accession = '(' + other.accession + ',' + self.accession + ')'
            self.seq = other.seq + 100 * 'N' + self.seq
        self.length = len(self.seq)
        self.Ncount = self.Ncount + other.Ncount + 100
        self.segment = ''


class Marker:
    def __init__(self, name, sequences, organism, qc, id, parent=None, run=None):
        # if parent is None:
        #    self.children = []
        self.children = []
        self.parent = parent
        self.name = name
        self.organism = organism
        self.query_cover = qc
        self.identity = id
        self.sequences = self._get_sequences(sequences)

        self.run = run

    def _get_sequences(self, sequences):
        seqs = []
        if isinstance(sequences[0], str):
            for sequence in sequences:
                if sequence.endswith('R'):
                    strd = 'minus'
                    sequence = sequence[:-1]
                else:
                    strd = 'plus'
                seqs.append(Sequence(self, sequence, strand=strd))
        elif isinstance(sequences[0], Sequence):
            seqs = sequences
        return seqs

    def to_fasta(self, acc=False):
        logger.info("Making fasta file for blasting...")
        if isinstance(self.organism, list):
            org = self.organism[0]
        else:
            org = self.organism
        fasta_file = open(org + self.name + '_acc.fasta', 'w')
        for sequence in self.sequences:
            fasta_file.write(sequence.fasta(acc))
        fasta_file.close()

    def _run_blast_command(self, blastn_cline, type):
        stderr = ''
        try:
            logger.info("Starting BLAST...")
            logger.info("BLAST command: " + str(blastn_cline))
            stdout, stderr = blastn_cline()
            if 'search aborted by Entrez' in stderr:
                logger.warning('Entrez aborted blastn search - probably empty query.')
        except:
            logger.critical(type +
                            ' blastn could not be properly executed. Check if it is installed or if you are connected '
                            'to internet.')
            if type == 'Standalone':
                logger.critical('If you used standalone version and do not have it configured, you may wish '
                                'to use -b, which will execute the remote version of blast databases. Be aware that it is '
                                'slower. Consider downloading and configuring the standalone database, as described in the '
                                'instruction manual.')
            logger.critical(stderr)
            clean_up(self.run.debug, self.run.folder, self.run.opt, 1)

    def _blast_search_using_query(self, file, query, out):
        blastn_cline = NcbiblastnCommandline(query=file, db="nt", task="blastn", remote=True, entrez_query=query,
                                             out=file.replace('_acc.fasta', '.txt'), outfmt=out, evalue=0.0001,
                                             max_target_seqs=50000, perc_identity=self.identity)
        self._run_blast_command(blastn_cline, 'Remote')

    def _blast_search_using_idlist(self, file, out):
        blastn_cline = NcbiblastnCommandline(query=file, db="nt", task="blastn", seqidlist='acc.txt',
                                             out=file.replace('_acc.fasta', '.txt'), outfmt=out, evalue=0.0001,
                                             max_target_seqs=50000, perc_identity=self.identity)
        self._run_blast_command(blastn_cline, 'Standalone')

    def blast_search(self, qr):
        if self.name == 'ITS' and self.run.its:
            out = '"6 sallseqid qcovs sstrand qlen qstart qend slen sstart send qaccver"'
        else:
            out = '"6 sallseqid qcovs sstrand"'

        if isinstance(self.organism, list):
            org = self.organism[0]
        else:
            org = self.organism

        file = org + self.name + "_acc.fasta"

        if qr:
            self._blast_search_using_query(file, qr, out)
        else:
            self._blast_search_using_idlist(file, out)

    def parse_blast(self):
        def open_blast_result():
            if isinstance(self.organism, list):
                org = self.organism[0]
            else:
                org = self.organism
            result = open(org + self.name + ".txt")
            result = result.readlines()
            result = sorted(result)
            result = list(result for result, _ in itertools.groupby(result))
            logger.info("Parsing BLAST results...")
            return result

        def get_ITS_record(rec, records, seqs):
            for s in self.sequences:
                if s.accession.split('.')[0] == rec[9].split('.')[0]:
                    if not s.segment:
                        qlen = int(rec[3])
                        slen = int(rec[6])
                        if abs(int(rec[5]) - int(rec[4])) > int(qlen) * 4 / 5:
                            if seqs in records.keys():
                                records[seqs][1] += '12'
                            else:
                                records[seqs] = [rec[2].lower(), '12']
                        else:
                            if int(rec[4]) <= int(qlen) / 3:  # its1
                                if slen - min(int(rec[7]), int(rec[8])) >= (int(qlen)) * 2 / 3:
                                    if seqs in records.keys():
                                        records[seqs][1] += '12'
                                    else:
                                        records[seqs] = [rec[2].lower(), '12']
                                else:
                                    if seqs in records.keys():
                                        records[seqs][1] += '1'
                                    else:
                                        records[seqs] = [rec[2].lower(), '1']
                            elif int(rec[4]) <= int(qlen) * 2 / 3:  # 5.8s
                                if min(int(rec[7]), int(rec[8])) >= int(qlen) / 3:
                                    if seqs in records.keys():
                                        records[seqs][1] += '1'
                                    else:
                                        records[seqs] = [rec[2].lower(), '1']
                                if slen - min(int(rec[7]), int(rec[8])) >= int(qlen) / 3:
                                    if seqs in records.keys():
                                        records[seqs][1] += '2'
                                    else:
                                        records[seqs] = [rec[2].lower(), '2']
                            else:  # its2
                                if min(int(rec[7]), int(rec[8])) >= int(qlen) * 2 / 3:
                                    if seqs in records.keys():
                                        records[seqs][1] += '12'
                                    else:
                                        records[seqs] = [rec[2].lower(), '12']
                                else:
                                    if seqs in records.keys():
                                        records[seqs][1] += '2'
                                    else:
                                        records[seqs] = [rec[2].lower(), '2']
                    else:
                        if seqs in records.keys():
                            records[seqs][1] += s.segment
                        else:
                            records[seqs] = [rec[2].lower(), s.segment]
            return records

        records = {}
        for line in open_blast_result():

            rec = line.split()
            seqs = tuple([r.split('|')[-2] for r in rec[0].split(';')])

            if self.run.remote:
                for s in seqs:
                    if s in self.run.blacklist:
                        continue

            if int(rec[1]) >= self.query_cover:
                if self.name == 'ITS' and self.run.its:
                    records = get_ITS_record(rec, records, seqs)

                else:
                    records[seqs] = [rec[2].lower(), None]

        self.sequences = []
        for record in records:
            for r in record:
                self.sequences.append(Sequence(self, r, records[record][0], segment=records[record][1]))

    def makecsv(self, its=None, dist=None):
        def remove_bad(sequences, bad):
            for b in bad:
                sequences.remove(b)
            return sequences

        def save_debug_files(rows):
            check(self.name + '.csv', 'r')
            file = open(self.name + '.csv', 'w', newline='')
            writer = csv.writer(file)
            writer.writerows(rows)
            file.close()

        rowlist = []
        bad = []
        itsjoinlist = []

        for s in self.sequences:
            r = s.assign_info(its, self.run.subsp)
            if r:
                if its and s.author:
                    itsjoinlist.append(s)
                else:
                    rowlist.append([s.species, s.tax, s.voucher, s.accession, s.length, str(s.seq)])
            else:
                bad.append(s)

        self.sequences = remove_bad(self.sequences, bad)

        if its:
            logger.info('ITS joiner is ON.')
            joinedlist = self.joinits(itsjoinlist, its, dist)
            rowlist = rowlist + joinedlist

        rows = sorted(rowlist, key=itemgetter(0, 1, 2))
        if self.run.debug:
            save_debug_files(rows)

        return rows

    def joinits(self, itslist, its, dist):
        joinrows = []
        added = set()
        for i in range(1, its + 1):
            itslist, added, join = self.joinpart(itslist, i, 'normal', dist)
            joinrows += join
            if len(added) != (2 * len(itslist)):
                itslist -= added
                itslist, added, join = self.joinpart(itslist, i, 'lownospace', dist)
                joinrows += join
                if len(added) != (2 * len(itslist)):
                    itslist -= added
                    itslist, added, join = self.joinpart(itslist, i, 'damerau', dist)
                    joinrows += join
        if len(added) != (2 * len(itslist)):
            for ss in itslist:
                if ss not in added:
                    if ss.segment == '1':
                        ss.accession = '(' + ss.accession + ',-)'
                    elif ss.segment == '2':
                        ss.accession = '(-,' + ss.accession + ')'
                    joinrows.append([ss.species, ss.tax, ss.voucher, ss.accession, ss.length, str(ss.seq)])
        return joinrows

    def joinpart(self, itslist, num, type, dist):
        def normal(num, sequence, sequence2, added, joinrows, itslist):
            if num < 3 and (
                        (sequence.voucher != '' and sequence.voucher == sequence2.voucher) or (
                                    sequence.voucher == '' and sequence.tax == sequence2.tax)):
                if num < 2 and sequence.author == sequence2.author:
                    if sequence.segment != sequence2.segment:
                        added.append(sequence)
                        added.append(sequence2)
                        sequence.join(sequence2)
                        joinrows.append([sequence.species, sequence.tax, sequence.voucher,
                                         sequence.accession,
                                         sequence.length, str(sequence.seq)])
                        self.sequences.remove(sequence2)
                        itslist.remove(sequence2)
                if num == 2:
                    if sequence.segment != sequence2.segment:
                        added.append(sequence)
                        added.append(sequence2)
                        sequence.join(sequence2)
                        joinrows.append([sequence.species, sequence.tax, sequence.voucher,
                                         sequence.accession,
                                         sequence.length, str(sequence.seq)])
                        self.sequences.remove(sequence2)
                        itslist.remove(sequence2)
            return added, joinrows, itslist

        def lownospace(num, sequence, sequence2, added, joinrows, itslist):
            if num < 3 and ((sequence.voucher != '' and "".join(
                    sequence.voucher.split('_')).lower() == "".join(
                sequence2.voucher.split('_')).lower()) or (
                            sequence.voucher == '' and sequence.tax == sequence2.tax)):
                if num < 2 and sequence.author == sequence2.author:
                    if sequence.segment != sequence2.segment:
                        added.append(sequence)
                        added.append(sequence2)
                        sequence.join(sequence2)
                        joinrows.append([sequence.species, sequence.tax, sequence.voucher,
                                         sequence.accession,
                                         sequence.length, str(sequence.seq)])
                        self.sequences.remove(sequence2)
                        itslist.remove(sequence2)
                if num == 2:
                    if sequence.segment != sequence2.segment:
                        added.append(sequence)
                        added.append(sequence2)
                        sequence.join(sequence2)
                        joinrows.append([sequence.species, sequence.tax, sequence.voucher,
                                         sequence.accession,
                                         sequence.length, str(sequence.seq)])
                        self.sequences.remove(sequence2)
                        itslist.remove(sequence2)
            return added, joinrows, itslist

        def damerau(num, sequence, sequence2, added, joinrows, itslist):
            if num < 3 and (
                        (sequence.voucher != '' and
                                 dameraulevenshtein(sequence.voucher.lower(),
                                                    sequence2.voucher.lower()) <= dist) or (
                                    sequence.voucher == '' and sequence.tax == sequence2.tax)):
                if num < 2 and sequence.author == sequence2.author:
                    if sequence.segment != sequence2.segment:
                        added.append(sequence)
                        added.append(sequence2)
                        sequence.join(sequence2)
                        joinrows.append([sequence.species, sequence.tax, sequence.voucher,
                                         sequence.accession,
                                         sequence.length, str(sequence.seq)])
                        self.sequences.remove(sequence2)
                        itslist.remove(sequence2)
                if num == 2:
                    if sequence.segment != sequence2.segment:
                        added.append(sequence)
                        added.append(sequence2)
                        sequence.join(sequence2)
                        joinrows.append([sequence.species, sequence.tax, sequence.voucher,
                                         sequence.accession,
                                         sequence.length, str(sequence.seq)])
                        self.sequences.remove(sequence2)
                        itslist.remove(sequence2)
            return added, joinrows, itslist

        def bruteforce(sequence, sequence2, added, joinrows, itslist):
            if sequence.segment != sequence2.segment:
                added.append(sequence)
                added.append(sequence2)
                sequence.join(sequence2)
                joinrows.append(
                    [sequence.species, sequence.tax, sequence.voucher, sequence.accession,
                     sequence.length,
                     str(sequence.seq)])
                self.sequences.remove(sequence2)
                itslist.remove(sequence2)
            return added, joinrows, itslist

        added = []
        joinrows = []
        for sequence in itslist.copy():
            if sequence not in added:
                for sequence2 in itslist.copy():
                    if sequence2 not in added:
                        if sequence.segment and sequence2.segment:
                            if sequence.species == sequence2.species:
                                if type == 'normal':
                                    added, joinrows, itslist = normal(num, sequence, sequence2, added, joinrows,
                                                                      itslist)
                                elif type == 'lownospace':
                                    added, joinrows, itslist = lownospace(num, sequence, sequence2, added, joinrows,
                                                                          itslist)

                                elif type == 'damerau':
                                    added, joinrows, itslist = damerau(num, sequence, sequence2, added, joinrows,
                                                                       itslist)

                                if num == 3:
                                    added, joinrows, itslist = bruteforce(sequence, sequence2, added, joinrows, itslist)

        return set(itslist), set(added), joinrows

    def breakdown(self):
        logger.info('Breaking down the markers by species')
        breakdic = {}
        for s in self.sequences:
            if s.species in breakdic.keys():
                breakdic[s.species].append(s)
            else:
                breakdic[s.species] = [s]
        for species in sorted(breakdic.keys()):
            self.children.append(
                Marker(self.name, breakdic[species], species, self.query_cover, self.identity, parent=self,
                       run=self.run))
        self.sequences = []

    def _blast_all_vs_all(self):
        self.to_fasta(acc=True)
        if isinstance(self.organism, list):
            org = self.organism[0]
        else:
            org = self.organism
        blastn_cline = NcbiblastnCommandline(query=org + self.name + '_acc.fasta',
                                             subject=org + self.name + '_acc.fasta', out='blast.txt',
                                             outfmt='"6 sseqid qseqid pident qcovs qcovhsp slen sseq"')
        stdout, stderr = blastn_cline()
        blast = open('blast.txt', 'r')
        lines = blast.readlines()
        lines = sorted(lines, key=itemgetter(0, 1))
        lines = list(lines for lines, _ in itertools.groupby(lines))
        return [line.split() for line in lines]

    def chooser(self):
        def choose_from_2():
            if (self.sequences[0].length - self.sequences[0].Ncount) > (
                        self.sequences[1].length - self.sequences[1].Ncount):
                # difference but taking into account number of Ns
                self.sequences[0].num = 2
                return self.sequences[0]
            else:
                self.sequences[1].num = 2
                return self.sequences[1]

        def choose_from_many(itsprioritize=False):

            def get_sumid_sumqc(sequences):
                sumids = 0
                sumqcs = 0

                for s in sequences:
                    for i in range(len(lines)):
                        if s.accession == lines[i][0]:
                            s.sumid += float(lines[i][2]) * (float(lines[i][4]) / float(lines[i][3]))
                            s.sumqc += float(lines[i][4])
                    sumqcs += s.sumqc
                    sumids += s.sumid
                return sequences, sumqcs, sumids

            def get_maxs(sequences):
                idmaxid = []
                idmaxqc = []

                for s in sequences:

                    if s.sumid >= 0.95 * sumids / len(self.sequences):
                        idmaxid.append(s)
                    if s.sumqc >= 0.95 * sumqcs / len(self.sequences):
                        idmaxqc.append(s)
                return set(idmaxid), set(idmaxqc)

            def choose_when_good_id(seqs):
                maxqc = 0
                maxid = 0
                maxlen = 0
                idmaxid2 = []
                idmaxqc2 = []
                for s in seqs:
                    if s.sumqc > maxqc:
                        idmaxqc2 = [s]
                        maxqc = s.sumqc
                    elif s.sumqc == maxqc:
                        idmaxqc2.append(s)
                if len(idmaxqc2) > 1:
                    for s in idmaxqc2:
                        if s.sumid > maxid:
                            idmaxid2 = [s]
                            maxid = s.sumid
                        elif s.sumid == maxid:
                            idmaxid2.append(s)
                    if len(idmaxid2) > 1:
                        idmaxlen = []
                        for s in idmaxid2:
                            if (s.length - s.Ncount) >= maxlen:
                                idmaxlen.append(s)
                                maxlen = (s.length - s.Ncount)
                                idmaxlen[0].num = len(self.sequences)
                                return idmaxlen[0]
                    else:
                        idmaxid2[0].num = len(self.sequences)
                        return idmaxid2[0]
                else:
                    idmaxqc2[0].num = len(self.sequences)
                    return idmaxqc2[0]

            def choose_when_good_qc(idmaxqc):
                if len(idmaxqc) > 1:
                    maxid = 0
                    idmaxid2 = []
                    for s in idmaxqc:
                        if s.sumid > maxid:
                            idmaxid2 = [s]
                            maxid = s.sumid
                        elif s.sumid == maxid:
                            idmaxid2.append(s)
                    if len(idmaxid2) > 1:
                        maxlen = 0
                        idmaxlen = []
                        for s in idmaxid2:
                            if (s.length - s.Ncount) >= maxlen:
                                idmaxlen.append(s)
                                maxlen = (s.length - s.Ncount)
                                idmaxlen[0].num = len(self.sequences)
                                return idmaxlen[0]
                    else:
                        idmaxid2[0].num = len(self.sequences)
                        return idmaxid2[0]
                else:
                    next(iter(idmaxqc)).num = len(self.sequences)
                    return next(iter(idmaxqc))

            lines = self._blast_all_vs_all()

            sequences = []

            if itsprioritize:  # prioritize joined its - run just on them
                for seq in self.sequences:
                    if not seq.segment:
                        # it means it is joined
                        sequences.append(seq)
                if not sequences:
                    return None
            else:
                sequences = self.sequences

            sequences, sumqcs, sumids = get_sumid_sumqc(sequences)
            idmaxid, idmaxqc = get_maxs(sequences)

            if idmaxid & idmaxqc:
                return choose_when_good_id(idmaxid & idmaxqc)
            else:
                return choose_when_good_qc(idmaxqc)

        if self.parent:
            if len(self.sequences) == 1:
                self.sequences[0].num = 1
                return self.sequences[0]
            elif len(self.sequences) == 2:
                return choose_from_2()
            else:
                if self.name == 'ITS' and self.run.its:
                    seq = choose_from_many(True)
                    if seq:
                        return seq
                return choose_from_many()
        else:
            for kid in self.children:
                self.sequences.append(kid.chooser())
            self.children = []


class Run:
    def __init__(self, opt='n', folder='.', debug=False):
        self.opt = opt
        self.folder = folder
        self.outfolder = None
        self.moddate = None
        self.remote = None
        self.markers = []
        self.today = str(datetime.datetime.now())[:10].replace('-', '/')
        self.blacklist = []
        self.synonyms = {}
        self.debug = debug
        self.organisms = None
        self.len_threshold = None
        self.its = None
        self.distance = None
        self.subsp = None
        self.excluded = None
        self.readconf()  # reading from config

    def readconf(self):
        def get_qc_or_id(arg, num_args, name):
            if arg:
                if len(arg) == 1:
                    argument = [arg[0] for _ in range(num_args)]
                elif len(arg) == num_args:
                    argument = arg
                else:
                    logger.critical('Wrong length of ' + name + ' parameters list!')
                    clean_up(self.debug, self.folder, self.opt, 1)
            else:
                argument = [0 for _ in range(num_args)]
            return argument

        startargs = readconfig(folder=self.folder)
        genes, seqs, self.organisms, self.len_threshold, self.its, qc, idy, self.distance, self.subsp, self.excluded, \
        self.remote, blacklist, synonyms, self.moddate, fin, upd, date = startargs

        query_cover = get_qc_or_id(qc, len(genes), 'query_cover')
        identity = get_qc_or_id(idy, len(genes), 'identity')

        for i in range(len(genes)):
            self.markers.append(
                Marker(genes[i], seqs[i], self.organisms, query_cover[i], identity[i], run=self))

        if self.opt == 'n':
            self.outfolder = self.folder + '/normal/'
        elif self.opt == 'u':
            self.outfolder = self.folder + '/update/'
            self.moddate = startargs[16]
            makeconfig(*startargs[:13], mdate=self.moddate, update='True', date=self.today, folder=self.folder)
        else:
            logger.critical('Wrong type of simulation')
            clean_up(self.debug, self.folder, self.opt, 1)

        check(self.outfolder, 'm')

        self.read_blacklist(blacklist)
        self.read_synonyms(synonyms)

    def read_blacklist(self, blacklist):
        if blacklist != '' and check(blacklist, 'e2'):
            with open(blacklist, 'r') as blist:
                blacklist = blist.readlines()
                blacklist = [line.split('#')[0].rstrip() for line in blacklist]
                black_fetch = Entrez.efetch(db="nucleotide", id=blacklist, rettype="acc", retmode="text", retmax=10000)
                self.blacklist = [acc.strip() for acc in black_fetch]
            black_fetch.close()
        else:
            self.blacklist = []

    def read_synonyms(self, synonyms):
        if synonyms != '' and check(synonyms, 'e2'):
            with open(synonyms, 'r') as slist:
                for line in slist:
                    line = line.split('#')[0]
                    syns = [x.strip(' \t\n') for x in line.split(',')]
                    for sy in syns[1:]:
                        if not self.subsp:
                            main_syn = "_".join(syns[0].split()[:2])
                        else:
                            main_syn = "_".join(syns[0].split())
                        self.synonyms[sy] = main_syn
            logger.info('Synonyms list: ' + str(self.synonyms))
        else:
            self.synonyms = []

    def track_synonyms(self, o):
        os = []
        for k, v in self.synonyms.items():
            if v == o:
                handle = Entrez.efetch(db="taxonomy", id=k, retmode="xml")
                tax_record = Entrez.read(handle)
                o_species = "_".join(tax_record[0]["ScientificName"].split()).lower()
                os.append(o_species)
        return os + [o]

    def set_term(self, o):
        if not o:
            orgs = "((" + ("[organism])OR(".join(self.organisms)) + "[organism]))"
            mdat = self.moddate
        else:
            os = self.track_synonyms(o)
            orgs = "((" + ("[organism])OR(".join(os)) + "[organism]))"
            mdat = '1000/01/01'

        if self.excluded:
            orgs += "NOT((" + ("[organism])OR(".join(self.excluded)) + "[organism]))"

        return '"(' + orgs + ")AND(" + mdat + ":" + self.today + "[mdat])AND(" + self.len_threshold + '[slen])"'

    def limit(self, o=None):

        def search_entrez(term):
            term = term.replace('"', '')
            handle = Entrez.esearch(db="nucleotide",
                                    term=term, retmax="100000")
            record = Entrez.read(handle)
            logger.info(record[
                            "Count"] + ' sequences total for ' + term.split('AND')[0] + 'and specified modification '
                                                                                        'date and sequence length '
                                                                                        'thresholds.')
            return record

        def report_no_new():
            if self.opt == 'u':
                upd = 'True'
            else:
                upd = 'False'

            startargs = readconfig(folder=self.folder)
            makeconfig(*startargs[:13], mdate=self.moddate, finished='True', update=upd, date=self.today,
                       folder=self.folder)

            logger.info('No new sequences found.')
            clean_up(self.debug, self.folder, self.opt, 0)

        def retrieve_additional_found(term):
            term = term.replace('"', '')
            records = record["IdList"]
            if int(record["Count"]) > 100000:
                if int(record["Count"]) % 100000:
                    rets = (int(record["Count"]) // 100000) + 1
                else:
                    rets = (int(record["Count"]) // 100000)
                for i in range(1, rets):
                    handle = Entrez.esearch(db="nucleotide",
                                            term=term, retmax=100000, retstart=i * 100000)
                    rec = Entrez.read(handle)
                    records += rec["IdList"]

            return records

        def translate_to_accs(count, records):
            accs = []
            if count % 10000:
                rets = (count // 10000) + 1
            else:
                rets = count // 10000
            for i in range(rets):
                acc_fetch = Entrez.efetch(db="nucleotide", id=records, rettype="acc", retmode="text", retmax="10000",
                                          retstart=str(i * 10000))
                accs += [acc.strip() for acc in acc_fetch]
                acc_fetch.close()
            return accs

        logger.info('Limiting search and excluding from blacklist...')

        prevdir = os.getcwd()
        check(self.outfolder + 'wf/', 'm')
        os.chdir(self.outfolder + 'wf/')

        term = self.set_term(o)
        record = search_entrez(term)

        if record["Count"] == '0':
            report_no_new()
        records = retrieve_additional_found(term)

        accs = translate_to_accs(int(record["Count"]), records)
        accver_file = open(self.outfolder + 'wf/acc.txt', 'w')
        for r in accs:
            if r not in self.blacklist:
                accver_file.write(r + '\n')
        accver_file.close()

        os.chdir(prevdir)

    def writefiles(self):
        def prepare_content():
            dict = {}
            headers = ['organism/marker']
            for i in range(len(self.markers)):
                self.markers[i].to_fasta()
                for sequence in self.markers[i].sequences:
                    if sequence.species in dict.keys():
                        dict[sequence.species][i] = (sequence.accession, sequence.num)
                    else:
                        dict[sequence.species] = [('-', 0) for _ in range(len(self.markers))]
                        dict[sequence.species][i] = (sequence.accession, sequence.num)
                headers.append(self.markers[i].name)
            return OrderedDict(sorted(dict.items())), headers

        dict, headers = prepare_content()

        matrixs = [open('accver_matrix.csv', 'w', newline=''), open('num_matrix.csv', 'w', newline='')]
        for j in range(len(matrixs)):
            writer = csv.writer(matrixs[j], delimiter=";")
            writer.writerow(headers)
            for orgs, pairs in dict.items():
                row = [orgs]
                for pair in pairs:
                    row.append(pair[j])
                writer.writerow(row)

    def update_matrices(self):
        updatedorgs = []
        logger.info('Updating the final matrixes.')
        matrixs = ['accver', 'num']
        for mat in matrixs:
            normal = open(self.folder + '/normal/' + mat + '_matrix.csv')
            update = open(self.folder + '/update/' + mat + '_matrix.csv')
            normalr = list(csv.reader(normal, delimiter=';'))
            updater = list(csv.reader(update, delimiter=';'))
            headers = updater[0]
            for nrow in normalr:
                for urow in updater:
                    if nrow[0] == urow[0]:
                        for i in range(1, len(nrow)):
                            if (urow[i] == '-' or urow[i] == '0') and (nrow[i] != urow[i]):
                                urow[i] = nrow[i]
                        updatedorgs.append(nrow[0])
            for nrow in normalr:
                if nrow[0] not in updatedorgs:
                    updater.append(nrow)
            updater = sorted(updater[1:], key=itemgetter(0))
            final = open(self.folder + '/final_' + mat + '_matrix.csv', 'w', newline='')
            writer = csv.writer(final, delimiter=';')
            writer.writerow(headers)
            writer.writerows(updater)

    def update_marker_files(self):
        for marker in self.markers:
            final = {}
            if isinstance(marker.organism, list):
                org = marker.organism[0]
            else:
                org = marker.organism
            oldfasta = open(self.folder + '/normal/' + org + marker.name + '_acc.fasta')
            newfasta = open(self.folder + '/update/' + org + marker.name + '_acc.fasta')
            old = oldfasta.readlines()
            new = newfasta.readlines()
            for i in range(len(old)):
                if old[i].startswith('>'):
                    final[old[i]] = old[i + 1]
            for j in range(len(new)):
                if new[j].startswith('>'):
                    final[new[j]] = new[j + 1]
            finalfasta = open(self.folder + '/' + org + marker.name + 'final.fasta', 'w')
            final = OrderedDict(sorted(final.items()))
            for orgs in final.keys():
                finalfasta.write(orgs)
                finalfasta.write(final[orgs])

    def start(self):
        def prepare():
            check(self.outfolder + 'wf/', 'm')
            os.chdir(self.outfolder + 'wf/')

            if self.remote:
                orgs = "((" + ("[organism])OR(".join(self.organisms)) + "[organism]))"
                if self.excluded:
                    orgs += "NOT((" + ("[organism])OR(".join(self.excluded)) + "[organism]))"
                return '"(' + orgs + ")AND(" + self.moddate + ":" + self.today + "[mdat])AND(" + self.len_threshold + \
                       '[slen])"'
            else:
                self.limit()
                return ''

        def run_main_analysis(qr):
            for marker in self.markers:
                marker.to_fasta(acc=True)
                marker.blast_search(qr)
                marker.parse_blast()
                if self.its and marker.name == 'ITS':
                    marker.makecsv(self.its, self.distance)
                else:
                    marker.makecsv()
                marker.breakdown()
                if self.opt == 'u':
                    for kidmarker in marker.children:
                        if self.remote:
                            qr = self.set_term(kidmarker.organism)
                        else:
                            self.limit(o=kidmarker.organism)
                            qr = ''
                        kidmarker.to_fasta(acc=True)
                        kidmarker.blast_search(qr)
                        kidmarker.parse_blast()
                        if self.its and marker.name == 'ITS':
                            kidmarker.makecsv(self.its, self.distance)
                        else:
                            kidmarker.makecsv()
                marker.chooser()

        def save_results():
            os.chdir(self.outfolder)
            self.writefiles()
            os.chdir(self.folder)

            if self.opt == 'u':
                self.update_matrices()
                self.update_marker_files()

            startargs = readconfig(folder=self.folder)
            makeconfig(*startargs[:13], mdate=self.moddate, finished='True', date=self.today, folder=self.folder)

        qr = prepare()
        run_main_analysis(qr)
        save_results()
