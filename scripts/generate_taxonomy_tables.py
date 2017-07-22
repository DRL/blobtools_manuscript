#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: analyse_blobtools_tables.py                 -t <FILE> [-r <FLOAT>] -d <DIR> [-m <FILE>] -b <FILE> --taxrank <STR> [-h|--help]

    Options:
        -h --help                                   show this
        -m, --mask_taxonomy <FILE>                  Mask true taxonomy of contigs in list (contig_id,mask)
        -t, --taxonomy <FILE>                       Taxonomy based on read mappings
        -r, --min_read_fraction <FLOAT>             Minimum fraction of reads to set true taxonomy based on read mappings [default: 1.0]
        -b, --blastorder <FILE>                     Order of tables in output (file = Name)
        -d, --dir <DIR>                             Directory containing blobDB tables (*.blobDB.table.txt)
        --taxrank <STR>                             Taxonomic rank at which to parse blobDB tables (order, phylum)

"""
from __future__ import division
import sys
from docopt import docopt
import os


def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")


class ContigObj():
    def __init__(self, contig_id, CELEG_reads, ECOLI_reads, HSAPI_reads, PAERU_reads, ref_id, true_taxonomy, read_fraction):
        self.contig_id = contig_id
        self.ref_id = ref_id
        self.true_taxonomy = true_taxonomy
        self.CELEG_reads = CELEG_reads
        self.ECOLI_reads = ECOLI_reads
        self.HSAPI_reads = HSAPI_reads
        self.PAERU_reads = PAERU_reads
        self.read_fraction = read_fraction
        self.length = 0
        self.gc = 0.0
        self.N = 0
        self.cov_A = 0.0
        self.cov_B = 0.0
        self.cov_sum = 0.0
        self.taxonomy_by_table = {}

    def add_table(self, table, length, gc, N, cov_A, cov_B, cov_sum, table_taxonomy):
        if self.length == 0:
            self.length = length
            self.gc = gc
            self.N = N
            self.cov_A = cov_A
            self.cov_B = cov_B
            self.cov_sum = cov_sum
        self.taxonomy_by_table[table] = table_taxonomy

    def get_data(self):
        return [self.contig_id, self.ref_id, self.true_taxonomy, self.CELEG_reads, self.ECOLI_reads, self.HSAPI_reads, self.PAERU_reads, self.read_fraction, self.length, self.gc, self.N, self.cov_A, self.cov_B, self.cov_sum]


class Main():
    def __init__(self, taxonomy_f, table_dir, mask_taxonomy_f, blastorder_f, taxrank, min_read_fraction):
        self.taxonomy_f = taxonomy_f
        self.mask_taxonomy_f = mask_taxonomy_f
        self.blastorder_f = blastorder_f
        self.table_dir = table_dir
        self.min_read_fraction = min_read_fraction
        self.taxrank = taxrank
        self.contigObjs_by_contig_id = {}
        self.table_order = []
        self.ref_id_order = ['CELEG', 'HSAPI', 'ECOLI', 'PAERU', 'Unknown']
        self.true_taxonomy_order_by_taxrank = {'phylum': ['Nematoda', 'Chordata', 'Proteobacteria'],
                                               'order': ['Rhabditida', 'Primates', 'Pseudomonadales', 'Enterobacterales']}
        self.true_order = ['Nematoda', 'Chordata', 'Proteobacteria']
        self.table_names = {}
        self.parameters = ['taxon_present', 'BLASTn', 'BLAST_MTS1', 'BLAST_MTS10', 'BLAST_MAXHSP1', 'BLAST_CUL10', 'Diamond', 'DIAMOND_MTS1', 'DIAMOND_MAXHSP1', 'DIAMOND_MTS10']
        self.parameters_by_table_name = {}
        self.parse_blast_order_f()
        self.masked_contigs = {}
        if self.mask_taxonomy_f:
            self.masked_contigs = self.parse_mask_f()
        self.parse_taxonomy_f()
        self.analyse_tables()
        self.blobtools_contig_table_f = None
        self.write_contig_table()
        self.contig_ids_by_ref_id = {}
        self.count_by_ref_id = {}
        self.span_by_ref_id = {}
        self.total_count = 0
        self.total_span = 0
        self.allocate_contigs_to_ref_ids()
        self.write_taxonomy_table()
        self.write_true_taxonomy_metrics()
        self.calculate_efficiency()

    def calculate_efficiency(self):
        count_by_taxonomy_by_table = {}
        span_by_taxonomy_by_table = {}
        # set_of_test_taxonomies = set(self.true_taxonomy_order_by_taxrank[self.taxrank] + ['global'])
        set_of_test_taxonomies = set(self.true_taxonomy_order_by_taxrank[self.taxrank])
        for table in self.table_order:
            count_by_taxonomy_by_table[table] = {}
            span_by_taxonomy_by_table[table] = {}
            for true_taxonomy in self.true_taxonomy_order_by_taxrank[self.taxrank]:
                count_by_taxonomy_by_table[table][true_taxonomy] = {}
                span_by_taxonomy_by_table[table][true_taxonomy] = {}
        for contig_id, contigObj in self.contigObjs_by_contig_id.items():
            length = contigObj.length
            contig_true_taxonomy = contigObj.true_taxonomy
            # condition : contig_true_taxonomy
            # testA : contig_table_taxonomy, excluding no-hit
            # testB : contig_table_taxonomy, including no-hit
            if contig_true_taxonomy is not 'Unknown':
                # print "# ", contig_id
                for table in self.table_order:
                    # print "- ", table
                    contig_table_taxonomy = contigObj.taxonomy_by_table[table]
                    for test_taxonomy in set_of_test_taxonomies:
                        result = self.determine_result(test_taxonomy, contig_true_taxonomy, contig_table_taxonomy)
                        # print "Test:%s, Cond:%s, Test:%s => %s" % (test_taxonomy, contig_true_taxonomy, contig_table_taxonomy, result)
                        count_by_taxonomy_by_table[table][test_taxonomy][result] = count_by_taxonomy_by_table[table][test_taxonomy].get(result, 0) + 1
                        span_by_taxonomy_by_table[table][test_taxonomy][result] = span_by_taxonomy_by_table[table][test_taxonomy].get(result, 0) + length
                    # global
                    # result = self.determine_result('global', contig_true_taxonomy, contig_table_taxonomy)
                    # count_by_taxonomy_by_table[table]['global'][result] = count_by_taxonomy_by_table[table]['global'].get(result, 0) + 1
                    # span_by_taxonomy_by_table[table]['global'][result] = span_by_taxonomy_by_table[table]['global'].get(result, 0) + length
        self.write_metrics(count_by_taxonomy_by_table, 'count')
        self.write_metrics(span_by_taxonomy_by_table, 'span')

    def write_metrics(self, metric_by_taxonomy_by_table, label):
        # order_of_groups = ['global'] + self.true_taxonomy_order_by_taxrank[self.taxrank]
        order_of_groups = self.true_taxonomy_order_by_taxrank[self.taxrank]
        output_precision_recall = []
        output_raw = []
        output_f_score = []
        header_precision_recall_1 = []
        header_raw_1 = []
        header_f_score = []
        header_done = False
        for parameter in self.parameters:
            header_precision_recall_1.append(parameter)
            header_f_score.append(parameter)
            header_raw_1.append(parameter)
        for table in self.table_order:
            line_precision_recall = []
            line_raw = []
            line_f_score = []
            for parameter in self.parameters:
                line_precision_recall.append(self.parameters_by_table_name[self.table_names[table]][parameter])
                line_raw.append(self.parameters_by_table_name[self.table_names[table]][parameter])
                line_f_score.append(self.parameters_by_table_name[self.table_names[table]][parameter])
            for group in order_of_groups:
                # if not group == 'global':
                #    accuracy = (metric_by_taxonomy_by_table[table][group].get('TP', 0) + metric_by_taxonomy_by_table[table][group].get('TN', 0)) / (metric_by_taxonomy_by_table[table][group].get('TP', 0) + metric_by_taxonomy_by_table[table][group].get('TN', 0) + metric_by_taxonomy_by_table[table][group].get('FP', 0) + metric_by_taxonomy_by_table[table][group].get('FN', 0))
                #    line.append(accuracy)
                # else:
                #    line.append("N/A")
                if not header_done:
                    header_precision_recall_1.append("%s Precision" % group)
                    header_precision_recall_1.append("%s Recall" % group)
                    header_raw_1.append("%s TP" % group)
                    header_raw_1.append("%s FP" % group)
                    header_raw_1.append("%s TN" % group)
                    header_raw_1.append("%s FN" % group)
                    header_f_score.append("%s F-score" %group)
                precision = metric_by_taxonomy_by_table[table][group].get('TP', 0) / (metric_by_taxonomy_by_table[table][group].get('TP', 0) + metric_by_taxonomy_by_table[table][group].get('FP', 0))
                line_precision_recall.append(precision)
                recall = metric_by_taxonomy_by_table[table][group].get('TP', 0) / (metric_by_taxonomy_by_table[table][group].get('TP', 0) + metric_by_taxonomy_by_table[table][group].get('FN', 0))
                line_precision_recall.append(recall)
                f_score = 2 * (precision * recall)/(precision + recall)
                line_f_score.append(f_score)
                line_raw.append(metric_by_taxonomy_by_table[table][group].get('TP', 0))
                line_raw.append(metric_by_taxonomy_by_table[table][group].get('FP', 0))
                line_raw.append(metric_by_taxonomy_by_table[table][group].get('TN', 0))
                line_raw.append(metric_by_taxonomy_by_table[table][group].get('FN', 0))
            if not header_done:
                output_precision_recall.append("\t".join(header_precision_recall_1))
                output_raw.append("\t".join(header_raw_1))
                output_f_score.append("\t".join(header_f_score))
                header_done = True
            output_precision_recall.append("\t".join([str(x) for x in line_precision_recall]))
            output_raw.append("\t".join([str(x) for x in line_raw]))
            output_f_score.append("\t".join([str(x) for x in line_f_score]))
        with open("efficiency.precision_recall.%s.min_read_frac_%s.%s.tsv" % (self.taxrank, self.min_read_fraction, label), 'w') as fh:
            fh.write("\n".join(output_precision_recall) + "\n")
        with open("efficiency.raw.%s.min_read_frac_%s.%s.tsv" % (self.taxrank, self.min_read_fraction, label), 'w') as fh:
            fh.write("\n".join(output_raw) + "\n")
        with open("efficiency.f_score.%s.min_read_frac_%s.%s.tsv" % (self.taxrank, self.min_read_fraction, label), 'w') as fh:
            fh.write("\n".join(output_f_score) + "\n")

    def determine_result(self, test_taxonomy, condition, test):
        if not test_taxonomy == 'global':
            if test_taxonomy == condition:
                if test == condition:
                    return "TP"
                else:
                    return "FN"
            else:
                if test == test_taxonomy:
                    return "FP"
                else:
                    return "TN"
        else:
            if condition == test:
                return "TP"
            else:
                if test in set(self.true_taxonomy_order_by_taxrank[self.taxrank]):
                    return "FP"
                else:
                    return "FN"


    def parse_blast_order_f(self):
        for line in read_file(self.blastorder_f):
            col = line.split(" = ")
            table = col[0]
            table_name = col[1]
            self.table_order.append(table)
            self.table_names[table] = table_name
            split_names = set(table_name.split())
            self.parameters_by_table_name[table_name] = {}
            for parameter in self.parameters:
                if parameter in split_names:
                    self.parameters_by_table_name[table_name][parameter] = 1
                else:
                    self.parameters_by_table_name[table_name][parameter] = 0

    def allocate_contigs_to_ref_ids(self):
        for contig_id in sorted(self.contigObjs_by_contig_id):
            contigObj = self.contigObjs_by_contig_id[contig_id]
            if contigObj.ref_id not in self.contig_ids_by_ref_id:
                self.contig_ids_by_ref_id[contigObj.ref_id] = []
            self.contig_ids_by_ref_id[contigObj.ref_id].append(contig_id)
            if contigObj.ref_id not in self.count_by_ref_id:
                self.count_by_ref_id[contigObj.ref_id] = 0
            self.count_by_ref_id[contigObj.ref_id] += 1
            if contigObj.ref_id not in self.span_by_ref_id:
                self.span_by_ref_id[contigObj.ref_id] = 0
            self.span_by_ref_id[contigObj.ref_id] += contigObj.length

    def write_true_taxonomy_metrics(self):
        output = []
        header = "#ref_id\tcount\tperc_count\tspan\tperc_span"
        output.append(header)
        for ref_id, count in self.count_by_ref_id.items():
            self.total_count += count
        for ref_id, span in self.span_by_ref_id.items():
            self.total_span += span
        for ref_id in self.ref_id_order:
            line = []
            line.append(ref_id)
            count = self.count_by_ref_id[ref_id]
            line.append(count)
            perc_count = "{0:.6f}".format(count / self.total_count)
            line.append(perc_count)
            span = self.span_by_ref_id[ref_id]
            line.append(span)
            perc_span = "{0:.6f}".format(span / self.total_span)
            line.append(perc_span)
            output.append("\t".join([str(x) for x in line]))
        with open("true_taxonomy_metrics_by_ref_id.%s.min_read_frac_%s.tsv" % (self.taxrank, self.min_read_fraction), 'w') as fh:
            fh.write("\n".join(output) + "\n")

    def write_taxonomy_table(self):
        counts_by_table_by_taxonomy_by_ref_id = {}
        span_by_table_by_taxonomy_by_ref_id = {}
        for ref_id in self.ref_id_order:
            set_of_taxonomies = set()
            counts_by_table_by_taxonomy_by_ref_id[ref_id] = {}
            span_by_table_by_taxonomy_by_ref_id[ref_id] = {}
            for contig_id in self.contig_ids_by_ref_id[ref_id]:
                contigObj = self.contigObjs_by_contig_id[contig_id]
                for table in self.table_order:
                    table_taxonomy = contigObj.taxonomy_by_table[table]
                    set_of_taxonomies.add(table_taxonomy)
            for taxonomy in list(set_of_taxonomies):
                counts_by_table_by_taxonomy_by_ref_id[ref_id][taxonomy] = {}
                span_by_table_by_taxonomy_by_ref_id[ref_id][taxonomy] = {}
                for table in self.table_order:
                    counts_by_table_by_taxonomy_by_ref_id[ref_id][taxonomy][table] = 0
                    span_by_table_by_taxonomy_by_ref_id[ref_id][taxonomy][table] = 0
            for contig_id in self.contig_ids_by_ref_id[ref_id]:
                contigObj = self.contigObjs_by_contig_id[contig_id]
                for table in self.table_order:
                    table_taxonomy = contigObj.taxonomy_by_table[table]
                    counts_by_table_by_taxonomy_by_ref_id[ref_id][table_taxonomy][table] += 1
                    span_by_table_by_taxonomy_by_ref_id[ref_id][table_taxonomy][table] += contigObj.length
        header = []
        header.append("#blobtools_taxonomy")
        for table in self.table_order:
            # header.append("%s count" % self.table_names[table])
            # header.append("%s count fract." % self.table_names[table])
            # header.append("%s span" % self.table_names[table])
            header.append("%s span fract." % self.table_names[table])
        for ref_id in self.ref_id_order:
            output = []
            output.append("\t".join(header))
            for taxonomy in sorted(counts_by_table_by_taxonomy_by_ref_id[ref_id]):
                line = []
                line.append(taxonomy)
                for table in self.table_order:
                    # count = counts_by_table_by_taxonomy_by_ref_id[ref_id][taxonomy][table]
                    # line.append(count)
                    # line.append("{0:.3f}".format(count / self.count_by_taxonomy[taxonomy]))
                    span = span_by_table_by_taxonomy_by_ref_id[ref_id][taxonomy][table]
                    # line.append(span)
                    line.append("{0:.3f}".format(span / self.span_by_ref_id[ref_id]))
                output.append("\t".join([str(x) for x in line]))
            with open("%s.%s.min_read_frac_%s.tsv" % (ref_id, self.taxrank, self.min_read_fraction), 'w') as fh:
                fh.write("\n".join(output) + "\n")

    def write_contig_table(self):
        output = []
        header = []
        header.append("#contig_id")
        header.append("ref_id")
        header.append("true_taxonomy")
        header.append("CELEG_reads")
        header.append("ECOLI_reads")
        header.append("HSAPI_reads")
        header.append("PAERU_reads")
        header.append("read_fraction")
        header.append("length")
        header.append("gc")
        header.append("N")
        header.append("cov_A")
        header.append("cov_B")
        header.append("cov_sum")
        for table in self.table_order:
            header.append(self.table_names[table])
        output.append("\t".join(header))
        for contig_id in sorted(self.contigObjs_by_contig_id):
            contigObj = self.contigObjs_by_contig_id[contig_id]
            line = contigObj.get_data()
            for table in self.table_order:
                line.append(contigObj.taxonomy_by_table[table])
            output.append("\t".join([str(x) for x in line]))
        self.blobtools_contig_table_f = "blobtools_contig_table.%s.min_read_frac_%s.tsv" % (self.taxrank, self.min_read_fraction)
        with open(self.blobtools_contig_table_f, 'w') as fh:
            fh.write("\n".join(output) + "\n")

    def analyse_tables(self):
        for filename in os.listdir(self.table_dir):
            if filename.endswith(".table.txt"):
                table = filename
                idx_by_key = {}
                for line in read_file(os.path.join(self.table_dir, table)):
                    if line.startswith("# name"):
                        col = line.split()
                        for idx, key in enumerate(col):
                            idx_by_key[key] = idx - 1
                    if not line.startswith('#'):
                        col = line.split()
                        contig_id = col[idx_by_key['name']]
                        length = int(col[idx_by_key['length']])
                        gc = float(col[idx_by_key['GC']])
                        N = int(col[idx_by_key['N']])
                        cov_A = float(col[idx_by_key['cov0']])
                        cov_B = float(col[idx_by_key['cov1']])
                        cov_sum = float(col[idx_by_key['cov_sum']])
                        table_taxonomy = None
                        if self.taxrank == 'phylum':
                            table_taxonomy = col[idx_by_key['phylum.t.12']]
                        elif self.taxrank == 'order':
                            table_taxonomy = col[idx_by_key['order.t.16']]
                        else:
                            sys.exit("[X] - wrong taxrank %s" (self.taxrank))
                        if contig_id not in self.contigObjs_by_contig_id:  # contig did not receive reads mapped
                            self.contigObjs_by_contig_id[contig_id] = ContigObj(contig_id, 0, 0, 0, 0, 'Unknown', 'Unknown', 'N/A')
                        self.contigObjs_by_contig_id[contig_id].add_table(table, length, gc, N, cov_A, cov_B, cov_sum, table_taxonomy)

    def parse_mask_f(self):
        masked_contigs = {}
        for line in read_file(self.mask_taxonomy_f):
            col = line.split(",")
            contig_id = col[0]
            mask = col[1]
            masked_contigs[contig_id] = mask
        return masked_contigs

    def parse_taxonomy_f(self):
        true_phylum_by_ref_id = {'CELEG': 'Nematoda', 'HSAPI': 'Chordata', 'ECOLI': 'Proteobacteria', 'PAERU': 'Proteobacteria', 'AMBIGOUS': 'AMBIGOUS', 'Unknown': 'Unknown'}
        true_order_by_ref_id = {'CELEG': 'Rhabditida', 'HSAPI': 'Primates', 'ECOLI': 'Enterobacterales', 'PAERU': 'Pseudomonadales', 'AMBIGOUS': 'AMBIGOUS', 'Unknown': 'Unknown'}
        for line in read_file(self.taxonomy_f):
            if not line.startswith("#") and not line.startswith("*"):
                col = line.split()
                contig_id = col[0]
                CELEG_reads = int(col[1])
                ECOLI_reads = int(col[2])
                HSAPI_reads = int(col[3])
                PAERU_reads = int(col[4])
                ref_id = self.masked_contigs.get(contig_id, col[5])
                # if taxonomy == "HS19" or taxonomy == "HSMT":
                #    taxonomy = "HSAPI"
                true_taxonomy = None
                if self.taxrank == 'phylum':
                    true_taxonomy = true_phylum_by_ref_id[ref_id]
                elif self.taxrank == 'order':
                    true_taxonomy = true_order_by_ref_id[ref_id]
                else:
                    sys.exit("[X] - unsupported taxrank %s" (self.taxrank))
                read_fraction = float(col[6])
                if read_fraction >= self.min_read_fraction:
                    contigObj = ContigObj(contig_id, CELEG_reads, ECOLI_reads, HSAPI_reads, PAERU_reads, ref_id, true_taxonomy, read_fraction)
                else:
                    contigObj = ContigObj(contig_id, CELEG_reads, ECOLI_reads, HSAPI_reads, PAERU_reads, 'Unknown', 'Unknown', read_fraction)
                self.contigObjs_by_contig_id[contig_id] = contigObj


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    mask_taxonomy_f = args['--mask_taxonomy']
    taxonomy_f = args['--taxonomy']
    table_dir = args['--dir']
    blastorder_f = args['--blastorder']
    try:
        min_read_fraction = float(args['--min_read_fraction'])
    except ValueError:
        sys.exit("[X] --min_read_fraction must be float")
    taxrank = args['--taxrank']
    Main(taxonomy_f, table_dir, mask_taxonomy_f, blastorder_f, taxrank, min_read_fraction)
