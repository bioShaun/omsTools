'''
Created on Nov 2, 2010

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import subprocess
import shutil

GTF_EMPTY_FIELD = '.'
GTF_ATTR_SEP = ';'
GTF_ATTR_TAGVALUE_SEP = ' '

# Map GENCODE gene type to overarching categories
GENCODE_CATEGORY_MAP = {
    'IG_C_gene': 'protein_coding',
    'IG_D_gene': 'protein_coding',
    'IG_J_gene': 'protein_coding',
    'IG_V_gene': 'protein_coding',
    'IG_LV_gene': 'protein_coding',
    'TR_C_gene': 'protein_coding',
    'TR_J_gene': 'protein_coding',
    'TR_V_gene': 'protein_coding',
    'TR_D_gene': 'protein_coding',
    'TEC': 'protein_coding',
    'nonsense_mediated_decay': 'protein_coding',
    'non_stop_decay': 'protein_coding',
    'retained_intron': 'protein_coding',
    'protein_coding': 'protein_coding',
    'ambiguous_orf': 'protein_coding',
    'Mt_rRNA': 'ncRNA',
    'Mt_tRNA': 'ncRNA',
    'miRNA': 'ncRNA',
    'misc_RNA': 'ncRNA',
    'rRNA': 'ncRNA',
    'snRNA': 'ncRNA',
    'snoRNA': 'ncRNA',
    'ribozyme': 'ncRNA',
    'sRNA': 'ncRNA',
    'scaRNA': 'ncRNA',
    'scRNA': 'ncRNA',
    'non_coding': 'ncRNA',
    'known_ncrna': 'ncRNA',
    '3prime_overlapping_ncrna': 'ncRNA',
    'vaultRNA': 'ncRNA',
    'processed_transcript': 'lncRNA',
    'lincRNA': 'lncRNA',
    'macro_lncRNA': 'lncRNA',
    'sense_intronic': 'lncRNA',
    'sense_overlapping': 'lncRNA',
    'antisense': 'lncRNA',
    'antisense_RNA': 'lncRNA',
    'bidirectional_promoter_lncRNA': 'lncRNA',
    'IG_pseudogene': 'pseudogene',
    'IG_C_pseudogene': 'pseudogene',
    'IG_J_pseudogene': 'pseudogene',
    'IG_V_pseudogene': 'pseudogene',
    'TR_V_pseudogene': 'pseudogene',
    'TR_J_pseudogene': 'pseudogene',
    'Mt_tRNA_pseudogene': 'pseudogene',
    'tRNA_pseudogene': 'pseudogene',
    'snoRNA_pseudogene': 'pseudogene',
    'snRNA_pseudogene': 'pseudogene',
    'scRNA_pseudogene': 'pseudogene',
    'rRNA_pseudogene': 'pseudogene',
    'misc_RNA_pseudogene': 'pseudogene',
    'miRNA_pseudogene': 'pseudogene',
    'pseudogene': 'pseudogene',
    'processed_pseudogene': 'pseudogene',
    'polymorphic_pseudogene': 'pseudogene',
    'retrotransposed': 'pseudogene',
    'transcribed_processed_pseudogene': 'pseudogene',
    'transcribed_unprocessed_pseudogene': 'pseudogene',
    'transcribed_unitary_pseudogene': 'pseudogene',
    'translated_processed_pseudogene': 'pseudogene',
    'translated_unprocessed_pseudogene': 'pseudogene',
    'unitary_pseudogene': 'pseudogene',
    'unprocessed_pseudogene': 'pseudogene',
    'XH': 'lncRNA',
    'XT': 'lncRNA',
    'SD': 'lncRNA',
    'SII': 'lncRNA',
    'SOI': 'lncRNA',
    'SPI': 'lncRNA',
    'SU': 'lncRNA',
    'XIE': 'lncRNA',
    'XII': 'lncRNA',
    'XOE': 'lncRNA',
    'XOI': 'lncRNA',
    'XPE': 'lncRNA',
    'XPI': 'lncRNA',
}

# gene_type classify priority

GENE_TYPE_PRIORITY = ('protein_coding', 'pseudogene', 'TUCP', 'lncRNA',
                      'ncRNA')


class GTFError(Exception):
    pass


def get_tr_type(tr_type):
    return GENCODE_CATEGORY_MAP.get(tr_type, tr_type)


def get_gene_type(tr_types):
    ''' label gene type according to GENE_TYPE_PRIORITY '''
    for each_type in tr_types:
        each_type = GENCODE_CATEGORY_MAP.get(each_type, each_type)
        if each_type in GENE_TYPE_PRIORITY:
            return each_type
    return 'unknown'


def sort_gtf(filename, output_file, tmp_dir=None):
    args = ["sort"]
    if tmp_dir is not None:
        args.extend(["-T", tmp_dir])
    args.extend(["-k1,1", "-k4,4n", "-k3,3r", filename])
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    return subprocess.call(args, stdout=open(output_file, "w"), env=myenv)


def merge_sort_gtf_files(gtf_files, output_file, tmp_dir=None):
    tmp_file = os.path.splitext(output_file)[0] + ".unsorted.gtf"
    outfh = open(tmp_file, "w")
    for filename in gtf_files:
        shutil.copyfileobj(open(filename), outfh)
    outfh.close()
    sort_gtf(tmp_file, output_file, tmp_dir)
    os.remove(tmp_file)


def parse_assembly_loci(line_iter):
    '''
    requires that GTF file has been sorted and formatted such that a
    single 'transcript' feature appears before individual 'exon'
    features such that transcript boundaries can be ascertained. this
    greatly simplifies parsing. using this function without appropriately
    formatted GTF files will result in undefined behavior
    '''

    def window_overlap(a, b):
        if a[0] != b[0]:
            return False
        return (a[1] <= b[2]) and (b[1] <= a[2])

    def get_intervals(line_iter):
        for line in line_iter:
            if line.startswith("#"):
                continue
            # read the essential part of the GTF line
            line = line.rstrip()
            fields = line.split('\t', 5)
            if len(fields) < 5:
                continue
            seqid = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            yield seqid, start, end, line

    try:
        interval_iter = get_intervals(line_iter)
        # initialize window
        seqid, start, end, line = interval_iter.next()
        window = [line]
        window_range = (seqid, start, end)
        # separate into loci
        for seqid, start, end, line in interval_iter:
            # check if next transcript is outside current window
            interval = (seqid, start, end)
            if not window_overlap(interval, window_range):
                # yield current window
                yield window
                # reset window
                window = [line]
                window_range = (seqid, start, end)
            else:
                # add transcript to window
                window.append(line)
                newstart = (start
                            if start < window_range[1] else window_range[1])
                newend = (end if end > window_range[2] else window_range[2])
                window_range = (seqid, newstart, newend)
    except StopIteration:
        pass
    # yield last window
    if len(window) > 0:
        yield window


def parse_loci(line_iter):
    '''
    return a dict of transcripts GTFFeature grouped by gene_id
    '''
    gene_dict = {}
    for eachline in line_iter:
        # skip annotation
        if eachline.startswith('#'):
            continue
        f = GTFFeature.from_string(eachline)
        gene_id = f.attrs['gene_id']
        gene_dict.setdefault(gene_id, []).append(eachline)
    return gene_dict


class GTFFeature(object):
    '''
    1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    5. end - The ending position of the feature (inclusive).
    6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
    7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    8. phase - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.

    chr1    Cufflinks       transcript      136546  137059  1000    .       .       gene_id "VCAP_SHEZH2.657699"; transcript_id "VCAP_SHEZH2.657699.1"; FPKM "100.7219943204"; frac "1.000000"; conf_lo "80.649925"; conf_hi "120.794064"; cov "2.198209";
    '''
    __slots__ = ('seqid', 'source', 'feature_type', 'start', 'end', 'score',
                 'strand', 'phase', 'attrs')

    def __str__(self):
        line = [
            self.seqid,
            self.source,
            self.feature_type,
            # convert to 1-based intervals
            str(self.start + 1),
            str(self.end),
            str(self.score),
            str(self.strand),
            self.phase
        ]
        attr_str = ' '.join(
            '%s "%s";' % (k, v) for (k, v) in self.attrs.iteritems())
        line.append(attr_str)
        return '\t'.join(line)

    @staticmethod
    def from_string(line, attr_defs=None):
        f = GTFFeature()
        # read the GTF line
        fields = line.strip().split('\t')
        f.seqid = fields[0]
        f.source = fields[1]
        f.feature_type = fields[2]
        # convert from 1-based (inclusive) to 0-based (exclusive) intervals
        f.start = int(fields[3]) - 1
        f.end = int(fields[4])
        f.score = 0 if (fields[5] == '.') else float(fields[5])
        strand = fields[6]
        if not (strand == '+' or strand == '-'):
            strand = GTF_EMPTY_FIELD
        f.strand = strand
        f.phase = fields[7]
        attrs = {}
        if fields[8] != GTF_EMPTY_FIELD:
            attr_strings = fields[8].split(GTF_ATTR_SEP)
            for a in attr_strings:
                a = a.strip()
                if len(a) == 0:
                    continue
                tag, value = a.split(GTF_ATTR_TAGVALUE_SEP, 1)
                # remove quotes
                value = value.strip('"')
                #value = value.split('"')[1]
                # apply parsing function
                # if (attr_defs != None) and (tag in attr_defs) and (attr_defs[tag] != None):
                #     value = attr_defs[tag](value)
                if (attr_defs is not None) and (tag not in attr_defs):
                    continue
                attrs[tag] = value
        f.attrs = attrs
        return f

    @staticmethod
    def parse(line_iter, attr_defs=None):
        for line in line_iter:
            # read the GTF line
            if not line:
                continue
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            yield GTFFeature.from_string(line, attr_defs)
