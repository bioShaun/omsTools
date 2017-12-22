'''
Created on Nov 30, 2011

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
import collections
from gtf import parse_loci, GTFFeature, GTFError

# attributes
TRANSCRIPT_ID = "transcript_id"
GENE_ID = "gene_id"

# strand constants
POS_STRAND = 0
NEG_STRAND = 1
NO_STRAND = 2

# strand format conversion
_strand_int_to_str = ["+", "-", "."]
_strand_str_to_int = {"+": POS_STRAND, "-": NEG_STRAND, ".": NO_STRAND}


def strand_str_to_int(strand):
    return _strand_str_to_int[strand]


def strand_int_to_str(strand):
    return _strand_int_to_str[strand]


def merge_strand(strand1, strand2):
    if strand1 == strand2:
        return strand1
    elif strand1 == NO_STRAND:
        return strand2
    elif strand2 == NO_STRAND:
        return strand1
    else:
        assert False


def cmp_strand(a, b):
    if (a == NO_STRAND) or (b == NO_STRAND):
        return True
    return a == b


def interval_overlap(a, b):
    return (a.start < b.end) and (b.start < a.end)


def interval_overlap_threshold(a, b, d=0):
    return (a.start < b.end + d) and (b.start < a.end + d)


class Exon(object):
    __slots__ = ('start', 'end')

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __str__(self):
        return ("%s-%s" % (self.start, self.end))

    def __repr__(self):
        return ("<%s(start=%d end=%d)>" % (self.__class__.__name__, self.start,
                                           self.end))

    def __eq__(self, other):
        return (self.start == other.start) and (self.end == other.end)

    def __ne__(self, other):
        return (self.start != other.start) or (self.end != other.end)

    def __lt__(self, other):
        return self.start < other.start

    def __le__(self, other):
        return self.start <= other.start

    def __gr__(self, other):
        return self.start > other.start

    def __ge__(self, other):
        return self.start >= other.start

    def __hash__(self):
        return (self.start << 16) | (self.end)

    def is_overlapping(self, other):
        return interval_overlap(self, other)


class Transcript(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'score', 'exons', 'attrs')

    def __init__(self):
        self.chrom = None
        self.start = -1
        self.end = -1
        self.strand = NO_STRAND
        self.score = 0.0
        self.exons = None
        self.attrs = {}

    def __str__(self):
        return ("<%s(chrom='%s', start='%d', end='%d', strand='%s', "
                "score='%s' exons='%s', attrs='%s'" %
                (self.__class__.__name__, self.chrom, self.start, self.end,
                 strand_int_to_str(self.strand), str(self.score), self.exons,
                 self.attrs))

    @property
    def length(self):
        return sum((e.end - e.start) for e in self.exons)

    def iterintrons(self):
        e1 = self.exons[0]
        for j in xrange(1, len(self.exons)):
            e2 = self.exons[j]
            yield e1.end, e2.start
            e1 = e2

    def introns(self):
        return list(self.iterintrons())

    # add function to return tss of transcript
    def tss(self, left, right):
        tss_start = self.start + left
        tss_end = self.start + right
        tss_start = tss_start if tss_start > 0 else 0
        tss_end = tss_end if tss_end < self.end else self.end
        return tss_start, tss_end

    def to_feature_bed6(self, tss_left=-2200, tss_right=500):
        bed6 = []
        tss_start, tss_end = self.tss(tss_left, tss_right)

        def add_feature_pos(start, end, feature):
            bed_line = [
                self.chrom,
                str(start),
                str(end),
                str(self.attrs["transcript_id"]), '0',
                strand_int_to_str(self.strand),
                str(self.attrs["gene_id"]), feature
            ]
            return '\t'.join(bed_line)

        bed6.append(add_feature_pos(tss_start, tss_end, 'tss'))
        for exon in self.exons:
            bed6.append(add_feature_pos(exon.start, exon.end, 'exon'))
        for intron in self.iterintrons():
            bed6.append(add_feature_pos(intron[0], intron[1], 'intron'))
        for each_line in bed6:
            yield each_line

    def to_bed12(self):
        block_sizes = []
        block_starts = []
        tr_start = self.start
        for e in self.exons:
            block_starts.append(e.start - tr_start)
            block_sizes.append(e.end - e.start)
        # write
        s = '\t'.join([
            self.chrom,
            str(self.start),
            str(self.end),
            str(self.attrs["transcript_id"]), '0',
            strand_int_to_str(self.strand),
            str(self.start),
            str(self.start), '0',
            str(len(self.exons)), ','.join(map(str, block_sizes)) + ',',
            ','.join(map(str, block_starts)) + ','
        ])
        return s

    def to_gtf_features(self, source=None, score=1000):
        if source is None:
            source = 'assemblyline'
        # transcript feature
        f = GTFFeature()
        f.seqid = self.chrom
        f.source = source
        f.feature_type = 'transcript'
        f.start = self.start
        f.end = self.end
        f.score = score
        f.strand = strand_int_to_str(self.strand)
        f.phase = '.'
        f.attrs = self.attrs
        features = [f]
        # exon features
        for i, e in enumerate(self.exons):
            f = GTFFeature()
            f.seqid = self.chrom
            f.source = source
            f.feature_type = 'exon'
            f.start = e.start
            f.end = e.end
            f.score = score
            f.strand = strand_int_to_str(self.strand)
            f.phase = '.'
            f.attrs = self.attrs.copy()
            f.attrs["exon_number"] = i
            features.append(f)
        return features


def to_formatted_gtf(lines, gtf_file, attr_defs=None):
    transcripts = collections.OrderedDict()
    for line in lines:
        feature = GTFFeature.from_string(line, attr_defs)
        # skip gene annotation in gtf files
        if "transcript_id" not in feature.attrs:
            continue
        t_id = feature.attrs["transcript_id"]
        # extract exon information
        if t_id not in transcripts:
            if feature.feature_type == "exon":
                t = Transcript()
                t.chrom = feature.seqid
                t.start = feature.start
                t.end = feature.end
                t.strand = strand_str_to_int(feature.strand)
                t.exons = [Exon(feature.start, feature.end)]
                t.attrs = dict()
                for each_attr in feature.attrs:
                    if 'exon' not in each_attr.lower():
                        t.attrs.update({each_attr: feature.attrs[each_attr]})
                transcripts[t_id] = t
        else:
            t.start = t.start if t.start <= feature.start else feature.start
            t.end = t.end if t.end >= feature.end else feature.end
            t.exons.append(Exon(feature.start, feature.end))
    with open(gtf_file, 'w') as gtf_output:
        for each_tr in transcripts:
            each_tr_obj = transcripts[each_tr]
            for each_feature in each_tr_obj.to_gtf_features():
                gtf_output.write('{gtf_line}\n'.format(
                    gtf_line=str(each_feature)))


def transcripts_from_gtf_lines(lines, attr_defs=None):
    transcripts = collections.OrderedDict()
    for line in lines:
        feature = GTFFeature.from_string(line, attr_defs)
        # skip gene annotation in gtf files
        if "transcript_id" not in feature.attrs:
            continue
        t_id = feature.attrs["transcript_id"]
        if t_id not in transcripts:
            if feature.feature_type != "transcript":
                raise GTFError(
                    "Feature type '%s' found before 'transcript' record: %s" %
                    (feature.feature_type, str(feature)))
            t = Transcript()
            t.chrom = feature.seqid
            t.start = feature.start
            t.end = feature.end
            # convert from string strand notation ("+", "-", ".")
            # to integer (0, 1)
            t.strand = strand_str_to_int(feature.strand)
            t.exons = []
            t.attrs = feature.attrs
            transcripts[t_id] = t
        else:
            t = transcripts[t_id]
        if feature.feature_type == "exon":
            t.exons.append(Exon(feature.start, feature.end))
    # sort transcript exons by genomic position
    for t in transcripts.itervalues():
        t.exons.sort()
    return transcripts.values()


def parse_gtf(fileh, attr_defs=None):
    gene_line_dict = parse_loci(fileh)
    for each_gene in gene_line_dict:
        locus_features = gene_line_dict[each_gene]
        yield each_gene, transcripts_from_gtf_lines(locus_features, attr_defs)
